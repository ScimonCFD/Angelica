"""
Sensitivity studies and nodal analysis (roadmap item 7.2).

Parameters are named by a target string:

  "pressure:<node>"            pressure of a PressureBoundary (Pa)
  "flow:<node>"                mass flow of a FlowBoundary (kg/s)
  "temperature:<node>"         temperature of a fixed-temperature ThermalBoundary (°C)
  "component:<id>.<field>"     any numeric field of a component, e.g.
                               "component:trunk.diameter_m",
                               "component:BPR.set_pressure_pa"

``sensitivity`` solves the case once per value (``sensitivity_2d`` on a
grid); every point keeps its ``SolveResult`` — converged or not, never an
exception from the solver.

``nodal_analysis`` splits a network at a well's wellhead:
  - inflow: the well's deliverability curve, wellhead pressure against
    rate (IPR and tubing, the curve the network solver itself uses);
  - outflow: the wellhead pressure the rest of the network needs to take a
    rate, from solving the network with the well replaced by that rate
    entering at its wellhead (with its fluid, at the tubing's wellhead
    temperature);
  - operating point: where the two meet on the stable branch of the
    inflow curve, refined by bisection on the outflow solves.
The operating point equals the well's rate in the full network solve
(``tests/test_studies.py``); that is the consistency of the split.
Not supported (ValueError): wells with gas lift (the lift gas joins at the
wellhead), a wellhead that is a boundary node, or one that also receives
other links (the split would seed their mixture with the well's fluid).
A crossing the bisection cannot refine (an outflow solve that does not
converge) gives no operating point, only a note.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from types import SimpleNamespace
from typing import Any, Callable, Sequence

from angelica.core.case import (
    FlowBoundary,
    NetworkCase,
    ThermalBoundary,
)

SolverFactory = Callable[[], object]


def set_parameter(case: NetworkCase, target: str, value: float) -> NetworkCase:
    """A copy of ``case`` with the parameter named by ``target`` set to ``value``."""
    kind, _, rest = target.partition(":")
    if kind == "pressure":
        node = int(rest)
        found = False

        def p(bcs):
            nonlocal found
            out = []
            for b in bcs:
                if b.node_id == node:
                    found = True
                    b = replace(b, pressure_pa=float(value))
                out.append(b)
            return tuple(out)
        new = replace(case, pressure_inlets=p(case.pressure_inlets), pressure_outlets=p(case.pressure_outlets))
    elif kind == "flow":
        node = int(rest)
        found = any(b.node_id == node for b in case.flow_inlets + case.flow_outlets)
        new = replace(case,
                      flow_inlets=tuple(replace(b, mass_flow_kg_per_s=float(value)) if b.node_id == node else b
                                        for b in case.flow_inlets),
                      flow_outlets=tuple(replace(b, mass_flow_kg_per_s=float(value)) if b.node_id == node else b
                                         for b in case.flow_outlets))
    elif kind == "temperature":
        node = int(rest)
        found = any(t.node_id == node and t.bc_type == "fixed_temperature" for t in case.thermal_inlets)
        new = replace(case, thermal_inlets=tuple(replace(t, temperature_c=float(value))
                                                 if t.node_id == node and t.bc_type == "fixed_temperature" else t
                                                 for t in case.thermal_inlets))
    elif kind == "component":
        import dataclasses
        cid, _, field_name = rest.rpartition(".")
        found = False
        comps = []
        for c in case.components:
            if c.component_id == cid:
                names = {f.name for f in dataclasses.fields(c)}
                if field_name not in names or field_name in ("start_node", "end_node", "component_id"):
                    raise ValueError(f"component {cid!r} has no numeric parameter {field_name!r}")
                found = True
                old = getattr(c, field_name)
                if isinstance(old, bool) or not isinstance(old, (int, float)):
                    raise ValueError(f"{cid}.{field_name} is not numeric")
                if isinstance(old, int) and float(value) != int(value):
                    raise ValueError(f"{cid}.{field_name} is an integer (got {value})")
                changes: dict[str, Any] = {field_name: int(value) if isinstance(old, int) else float(value)}
                try:
                    c = replace(c, **changes)
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"{cid}.{field_name} = {value}: {exc}") from None
            comps.append(c)
        new = replace(case, components=tuple(comps))
    else:
        raise ValueError(f"unknown target {target!r} (use pressure:, flow:, temperature: or component:)")
    if not found:
        raise ValueError(f"target {target!r} not found in the case")
    return new


@dataclass(frozen=True)
class StudyPoint:
    values: tuple[float, ...]
    result: Any                        # SolveResult, or None when the point raised
    error: str = ""                    # why there is no result (invalid value, solver exception)

    @property
    def converged(self) -> bool:
        return self.result is not None and bool(self.result.converged)

    def flow(self, component_id: str) -> float:
        """Mass flow (kg/s) through a component of the case."""
        for f in self.result.component_flows:
            if f.label.split(":")[-1] == component_id:
                return float(f.mass_flow_kg_per_s)
        raise KeyError(component_id)

    def pressure(self, node_id: int) -> float:
        return float(self.result.node_pressures_pa[node_id])


def _solve_point(case: NetworkCase, solver_factory: SolverFactory, settings, values) -> StudyPoint:
    """One study point; an invalid value or a solver exception is kept as the
    point's ``error`` instead of stopping the study."""
    try:
        c = case
        for target, value in settings:
            c = set_parameter(c, target, value)
        return StudyPoint(tuple(values), solver_factory().solve(c))  # type: ignore[attr-defined]
    except Exception as exc:                 # noqa: BLE001 — reported per point
        return StudyPoint(tuple(values), None, f"{type(exc).__name__}: {exc}")


def sensitivity(case: NetworkCase, solver_factory: SolverFactory, target: str,
                values: Sequence[float]) -> list[StudyPoint]:
    """One solve per value of the parameter ``target``.  The target is
    checked first (an unknown target raises ValueError); a value the
    component rejects, or a solver exception, is kept in that point."""
    set_parameter(case, target, values[0]) if values else None
    return [_solve_point(case, solver_factory, [(target, float(v))], (float(v),)) for v in values]


def sensitivity_2d(case: NetworkCase, solver_factory: SolverFactory, target_1: str, values_1: Sequence[float],
                   target_2: str, values_2: Sequence[float]) -> list[StudyPoint]:
    """One solve per pair of values, ``values_1`` outer."""
    if values_1 and values_2:
        set_parameter(set_parameter(case, target_1, values_1[0]), target_2, values_2[0])
    return [_solve_point(case, solver_factory, [(target_1, float(v1)), (target_2, float(v2))], (float(v1), float(v2)))
            for v1 in values_1 for v2 in values_2]


@dataclass(frozen=True)
class NodalAnalysis:
    well_id: str
    wellhead_node: int
    inflow: tuple[tuple[float, float], ...]           # (kg/s, wellhead Pa) of the well's curve
    outflow: tuple[tuple[float, float | None], ...]   # (kg/s, wellhead Pa needed; None: not converged)
    operating_point: tuple[float, float] | None       # (kg/s, Pa) on the stable branch, or None
    note: str = ""


def _well_inputs(case: NetworkCase, well, isothermal_c):
    """Reservoir pressure, temperature and the fluid of a well, as the solver
    sees them: the case fluid's correlation and gas properties with the
    reservoir node's black-oil parameters or composition; the reservoir
    temperature of its fixed-temperature boundary, the solver's isothermal
    temperature, or (no boundary) the solver's initial temperature."""
    from angelica.solvers.base import BaseSolver
    res = well.start_node
    p_res = next((b.pressure_pa for b in case.pressure_inlets if b.node_id == res), None)
    if p_res is None:
        raise ValueError(f"well {well.component_id!r}: its reservoir node {res} has no pressure boundary")
    t_bc = next((t.temperature_c for t in case.thermal_inlets
                 if t.node_id == res and t.bc_type == "fixed_temperature"), None)
    if isothermal_c is not None:
        t_res = float(isothermal_c)
    elif t_bc is not None:
        t_res = float(t_bc)
    else:
        t_res = BaseSolver._initial_temperature(case)
    fluid = case.fluid_model
    bc = next((b for b in case.inlet_fluid_bcs if b.node_id == res), None)
    if bc is not None:
        from angelica.properties.black_oil import BlackOilFluid
        params: dict[str, Any] = {"api_gravity": bc.api_gravity, "gas_gravity": bc.gas_gravity,
                  "gor_sc_m3_per_m3": bc.gor_sc_m3_per_m3, "wor_sc_m3_per_m3": bc.wor_sc_m3_per_m3}
        fluid = replace(fluid, **params) if isinstance(fluid, BlackOilFluid) else BlackOilFluid(**params)
    cbc = next((b for b in case.inlet_composition_bcs if b.node_id == res), None)
    if cbc is not None:
        from angelica.properties.compositional_fluid import CompositionalFluid
        fluid = CompositionalFluid(list(getattr(case.fluid_model, "component_names")), list(cbc.zs),
                                   eos_name=getattr(case.fluid_model, "eos_name", "PR"),
                                   hydrate_inhibitor=getattr(case.fluid_model, "hydrate_inhibitor", None))
    return float(p_res), float(t_res), fluid, bc, cbc


def _outflow_case(case: NetworkCase, well, m: float, t_wh: float, bc, cbc) -> NetworkCase:
    """The network without the well, its rate entering at the wellhead."""
    wh = well.end_node
    comps = tuple(c for c in case.components if c.component_id != well.component_id)
    used = {n for c in comps for n in (c.start_node, c.end_node)}
    res = well.start_node
    kw: dict[str, Any] = {}
    if bc is not None:      # the reservoir node keeps its fluid when other links still draw from it
        kw["inlet_fluid_bcs"] = tuple(b for b in case.inlet_fluid_bcs
                                      if b.node_id != res or res in used) + (replace(bc, node_id=wh),)
    if cbc is not None:
        kw["inlet_composition_bcs"] = tuple(b for b in case.inlet_composition_bcs
                                            if b.node_id != res or res in used) + (replace(cbc, node_id=wh),)
    return replace(
        case,
        components=comps,
        pressure_inlets=tuple(b for b in case.pressure_inlets if b.node_id in used),
        flow_inlets=tuple(case.flow_inlets) + (FlowBoundary(node_id=wh, mass_flow_kg_per_s=float(m)),),
        thermal_inlets=tuple(t for t in case.thermal_inlets if t.node_id in used and t.node_id != wh)
        + (ThermalBoundary(node_id=wh, temperature_c=float(t_wh)),),
        node_ids=(),
        **kw,
    )


def nodal_analysis(case: NetworkCase, solver_factory: SolverFactory, well_id: str,
                   rates_kg_per_s: Sequence[float] | None = None, n_rates: int = 12,
                   tolerance_rel: float = 1e-4) -> NodalAnalysis:
    """Inflow and outflow curves at a well's wellhead and their crossing."""
    from angelica.closures.well_model import well_curve, well_tubing_state
    from angelica.core.components import Well

    well = next((c for c in case.components if c.component_id == well_id), None)
    if not isinstance(well, Well):
        raise ValueError(f"{well_id!r} is not a Well of the case")
    if getattr(well, "has_gas_lift", False):
        raise ValueError(f"well {well_id!r} has gas lift: nodal analysis does not support it")
    wh = well.end_node
    boundary = {b.node_id for b in tuple(case.pressure_inlets) + tuple(case.pressure_outlets)
                + tuple(case.flow_inlets) + tuple(case.flow_outlets)}
    if wh in boundary:
        raise ValueError(f"well {well_id!r}: its wellhead node {wh} is a boundary node; nodal analysis needs "
                         "the network downstream of the wellhead")
    others_in = [c.component_id for c in case.components if c is not well and c.end_node == wh]
    if others_in:
        raise ValueError(f"well {well_id!r}: its wellhead node {wh} also receives {', '.join(others_in)}; "
                         "nodal analysis splits at a node fed by the well alone")
    probe = solver_factory()
    isothermal_c = None
    for name in ("black_oil_settings", "compositional_settings", "compressible_settings"):
        isothermal_c = getattr(getattr(probe, name, None), "isothermal_temperature_c", None)
        if isothermal_c is not None:
            break
    p_res, t_res, fluid, bc, cbc = _well_inputs(case, well, isothermal_c)
    link = SimpleNamespace(isothermal=isothermal_c is not None, zs=tuple(getattr(cbc, "zs", ()) or ()))
    curve = well_curve(well, fluid, p_res, t_res, link)
    m_per = curve.mass_per_oil_sc_m3 / 86_400.0
    inflow = tuple((q * m_per, p) for q, p in zip(curve.oil_rates_sc_m3_per_day, curve.wellhead_pressures_pa))
    peak = curve.peak_index
    m_max = inflow[-1][0]
    m_min = inflow[peak][0]
    if m_max <= 0.0:
        return NodalAnalysis(well_id, wh, inflow, (), None, "the well cannot flow at any wellhead pressure")

    def p_in(m: float) -> float:            # stable branch, linear in the tabulated curve
        pts = inflow[peak:]
        for (m0, p0), (m1, p1) in zip(pts, pts[1:]):
            if m0 <= m <= m1:
                return p0 + (p1 - p0) * (m - m0) / (m1 - m0) if m1 > m0 else p0
        return pts[-1][1] if m > pts[-1][0] else pts[0][1]

    def p_out(m: float) -> float | None:
        if isothermal_c is not None:
            t_wh = t_res
        else:
            st = well_tubing_state(well, fluid, p_res, t_res, m, link_state=link)
            t_wh = st.wellhead_temperature_c if st is not None else t_res
        try:
            r = solver_factory().solve(_outflow_case(case, well, m, t_wh, bc, cbc))  # type: ignore[attr-defined]
        except Exception:                    # noqa: BLE001 — a failed outflow point is reported as None
            return None
        return float(r.node_pressures_pa[wh]) if r.converged else None

    if rates_kg_per_s is None:
        lo = max(m_min, 1e-6)
        rates_kg_per_s = [lo + (m_max - lo) * i / (n_rates - 1) for i in range(n_rates)]
    outflow = tuple((float(m), p_out(float(m))) for m in rates_kg_per_s)

    # bracket a sign change of p_in − p_out on the stable branch
    pts = [(m, p) for m, p in outflow if p is not None and m >= m_min]
    op = None
    note = ""
    for (m0, q0), (m1, q1) in zip(pts, pts[1:]):
        d0, d1 = p_in(m0) - q0, p_in(m1) - q1
        if d0 == 0.0:
            op = (m0, q0)
            break
        if d0 * d1 < 0.0:
            a, b, da = m0, m1, d0
            refined = False
            for _ in range(60):
                mid = 0.5 * (a + b)
                q = p_out(mid)
                if q is None:
                    break
                d = p_in(mid) - q
                if da * d <= 0.0:
                    b = mid
                else:
                    a, da = mid, d
                if b - a <= tolerance_rel * max(b, 1e-9):
                    refined = True
                    break
            if refined:
                mid = 0.5 * (a + b)
                op = (mid, p_in(mid))
            else:
                note = (f"crossing between {a:.4g} and {b:.4g} kg/s, not refined: the outflow solve did not "
                        "converge inside that interval")
            break
    if op is None and not note:
        note = ("no crossing on the stable branch: the well cannot flow into this network"
                if pts and all(p_in(m) < q for m, q in pts) else "no crossing in the rates studied")
    return NodalAnalysis(well_id, wh, inflow, outflow, op, note)



@dataclass(frozen=True)
class SizingCandidate:
    diameter_m: float
    point: StudyPoint
    pressure_drop_pa: float | None        # start node − end node pressure of the pipe
    velocity_m_per_s: float | None        # mixture velocity (|volumetric flow| / area)
    erosional_ratio: float | None         # API RP 14E, two-phase pipes only
    passes: bool
    reasons: tuple[str, ...]              # why it fails (empty when it passes)


@dataclass(frozen=True)
class LineSizing:
    pipe_id: str
    candidates: tuple[SizingCandidate, ...]

    @property
    def selected(self) -> SizingCandidate | None:
        """The smallest candidate diameter that meets every criterion."""
        ok = [c for c in self.candidates if c.passes]
        return min(ok, key=lambda c: c.diameter_m) if ok else None


def size_line(case: NetworkCase, solver_factory: SolverFactory, pipe_id: str, diameters_m: Sequence[float],
              max_pressure_drop_pa: float | None = None, max_velocity_m_per_s: float | None = None,
              max_erosional_ratio: float | None = None,
              min_node_pressure: tuple[int, float] | None = None) -> LineSizing:
    """Line sizing: the network solved with each candidate diameter of one
    pipe line, checked against the criteria given; ``selected`` is the
    smallest diameter that meets them all.  A candidate whose solve does not
    converge fails ("not converged").

    ``pipe_id`` names a Pipe, or a line the GUI split into segments
    ("<id>_seg0", "<id>_seg1", …): every segment gets the diameter, and the
    pressure drop is from the first segment's start to the last one's end.
    The velocity is the mixture velocity at the pipe's mean conditions (the
    largest over the segments); in a gas line the outlet velocity is higher
    (15–30 % in the tests' 70 → 40 bar line).  The erosional-velocity ratio
    (API RP 14E) exists only for two-phase pipes: asked of a single-phase
    line, a candidate fails with "not applicable"."""
    from angelica.core.components import Pipe
    segs = [c for c in case.components if isinstance(c, Pipe)
            and (c.component_id == pipe_id or c.component_id.startswith(pipe_id + "_seg"))]
    if not segs:
        raise ValueError(f"{pipe_id!r} is not a pipe (or pipe line) of the case")
    if all(v is None for v in (max_pressure_drop_pa, max_velocity_m_per_s, max_erosional_ratio, min_node_pressure)):
        raise ValueError("give at least one sizing criterion")
    for name, v in (("max_pressure_drop_pa", max_pressure_drop_pa), ("max_velocity_m_per_s", max_velocity_m_per_s),
                    ("max_erosional_ratio", max_erosional_ratio)):
        if v is not None and not v > 0.0:
            raise ValueError(f"{name} must be positive (got {v})")
    nodes = {n for c in case.components for n in (c.start_node, c.end_node)}
    if min_node_pressure is not None and min_node_pressure[0] not in nodes:
        raise ValueError(f"node {min_node_pressure[0]} is not in the network")
    ids = [c.component_id for c in segs]
    first_start = segs[0].start_node
    last_end = segs[-1].end_node
    out = []
    for D in sorted(float(d) for d in diameters_m):
        pt = _solve_point(case, solver_factory, [(f"component:{i}.diameter_m", D) for i in ids], (D,))
        reasons: list[str] = []
        dp = vel = ero = None
        if not pt.converged:
            reasons.append("not converged" if pt.result is not None else f"error: {pt.error}")
        else:
            res = pt.result
            dp = pt.pressure(first_start) - pt.pressure(last_end)
            area = 3.141592653589793 * D * D / 4.0
            vels, eros = [], []
            for f in res.component_flows:
                if f.label.split(":")[-1] not in ids:
                    continue
                if f.multiphase is not None:
                    vels.append(abs(f.multiphase.mixture_velocity_m_per_s))
                    eros.append(f.multiphase.erosional_velocity_ratio)
                else:
                    vels.append(abs(f.volumetric_flow_m3_per_h) / 3600.0 / area)
            vel = max(vels) if vels else None
            ero = max(eros) if eros else None
            if max_pressure_drop_pa is not None and abs(dp) > max_pressure_drop_pa:
                reasons.append(f"pressure drop {abs(dp):.4g} Pa > {max_pressure_drop_pa:.4g}")
            if max_velocity_m_per_s is not None and vel is not None and vel > max_velocity_m_per_s:
                reasons.append(f"mean velocity {vel:.4g} m/s > {max_velocity_m_per_s:.4g}")
            if max_erosional_ratio is not None:
                if ero is None:
                    reasons.append("erosional ratio not applicable (single-phase pipe)")
                elif ero > max_erosional_ratio:
                    reasons.append(f"erosional ratio {ero:.3g} > {max_erosional_ratio:.3g}")
            if min_node_pressure is not None:
                node, p_min = min_node_pressure
                if pt.pressure(node) < p_min:
                    reasons.append(f"pressure at node {node} {pt.pressure(node):.4g} Pa < {p_min:.4g}")
        out.append(SizingCandidate(D, pt, dp, vel, ero, not reasons, tuple(reasons)))
    return LineSizing(pipe_id, tuple(out))
