"""
Calibration factors and line sizing (roadmap item 7.3, v1.33.0).

Validation (ANALYTICAL / consistency): the friction pressure drop is
f·(L/D)·ρv²/2, so at a fixed flow through a horizontal pipe (same Reynolds
number, same f) multiplying f by c multiplies the pressure drop by c
exactly — turbulent (Colebrook) and laminar (64/Re).  Tolerance 1e-6: the
outer and velocity loops' convergence (obtained ≤ 5.9e-9).  In Beggs-Brill the
holdup and the two-phase friction factor are multiplied as stated (exact
at the gradient level); in a network the response is monotone.  Factor 1
reproduces every existing result (the tutorial suite is unchanged).  No
field data were used: the factors are the user's to fit.
"""
from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT / "src", ROOT / "tests"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from angelica.core.case import FlowBoundary, NetworkCase, PressureBoundary  # noqa: E402
from angelica.core.components import Pipe  # noqa: E402
from angelica.properties.single_component import SingleComponentFluid  # noqa: E402
from angelica.solvers import SteadyIsothermalIncompressibleSolver  # noqa: E402


def liquid_line(mu: float, c: float, m: float) -> NetworkCase:
    return NetworkCase(
        name="line", fluid_model=SingleComponentFluid(density_kg_per_m3=850.0, viscosity_pa_s=mu),
        pressure_inlets=(), flow_inlets=(FlowBoundary(1, m),), pressure_outlets=(PressureBoundary(2, 1e5),),
        components=(Pipe(component_id="p", start_node=1, end_node=2, diameter_m=0.2, length_m=5000.0,
                         absolute_roughness_m=46e-6, friction_factor_multiplier=c),))


@pytest.mark.parametrize("mu, m", [(1e-3, 40.0), (0.5, 5.0)])     # turbulent (Re 2.5e5), laminar (Re 64)
@pytest.mark.parametrize("c", [0.8, 1.2, 1.5])
def test_friction_multiplier_scales_the_pressure_drop_exactly(mu, m, c):
    s = SteadyIsothermalIncompressibleSolver()
    base = s.solve(liquid_line(mu, 1.0, m))
    cal = s.solve(liquid_line(mu, c, m))
    dp0 = base.node_pressures_pa[1] - 1e5
    dp1 = cal.node_pressures_pa[1] - 1e5
    assert dp1 / dp0 == pytest.approx(c, rel=1e-6)


def test_beggs_brill_gradient_multipliers():
    from angelica.closures.beggs_brill import beggs_brill_gradient
    args = (2.0, 0.3, 750.0, 40.0, 2e-3, 1.2e-5, 0.1524, 46e-6, 0.05)
    base = beggs_brill_gradient(*args)
    f = beggs_brill_gradient(*args, friction_multiplier=1.3)
    h = beggs_brill_gradient(*args, holdup_multiplier=1.2)
    assert f.dP_dz_friction == pytest.approx(1.3 * base.dP_dz_friction, rel=1e-12)
    assert f.holdup == base.holdup
    assert h.holdup == pytest.approx(min(1.2 * base.holdup, 1.0), rel=1e-12)
    assert h.dP_dz_gravity > base.dP_dz_gravity                      # more liquid in an uphill pipe
    full = beggs_brill_gradient(*args, holdup_multiplier=10.0)
    assert full.holdup == 1.0                                         # kept within [0, 1]


def test_network_response_to_the_multipliers():
    """The two-loop black-oil network with Beggs-Brill (horizontal).  More
    friction, less flow, monotonically (obtained 25.71 → 25.05 → 22.35 kg/s
    at 1, 1.05, 1.3).  More holdup, MORE flow here (31.18 kg/s at 1.3): in a
    horizontal pipe the holdup enters only through Beggs-Brill's two-phase
    friction factor f_tp = f_n·e^S, S a function of y = C_L/H_L², and a
    larger H_L lowers it.  In an uphill pipe the heavier in-situ mixture
    adds gravity (test_beggs_brill_gradient_multipliers)."""
    from _two_loop_networks import build_case

    from angelica.closures import BeggsBrillCorrelation
    from angelica.closures.beggs_brill import two_phase_friction_factor
    from angelica.solvers import SteadyBlackOilSolver
    case = build_case("black_oil")
    flows = {}
    for key, kw in (("base", {}), ("f1.05", {"friction_factor_multiplier": 1.05}),
                    ("f1.3", {"friction_factor_multiplier": 1.3}), ("h1.3", {"holdup_multiplier": 1.3})):
        c = replace(case, components=tuple(replace(p, **kw) for p in case.components))
        r = SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation()).solve(c)
        assert r.converged
        flows[key] = next(f.mass_flow_kg_per_s for f in r.component_flows if f.label.endswith("export"))
    assert flows["f1.3"] < flows["f1.05"] < flows["base"]
    assert flows["h1.3"] > flows["base"]
    assert two_phase_friction_factor(0.02, 0.3, 0.6) < two_phase_friction_factor(0.02, 0.3, 0.5)   # the mechanism


def test_invalid_multipliers():
    for bad in ({"friction_factor_multiplier": 0.0}, {"holdup_multiplier": -1.0},
                {"friction_factor_multiplier": float("nan")}):
        with pytest.raises(ValueError):
            Pipe(component_id="p", start_node=1, end_node=2, diameter_m=0.2, length_m=10.0, **bad)


def test_line_sizing_selects_the_smallest_diameter_that_meets_the_criteria():
    from _two_loop_networks import build_case, solver_for

    from angelica.studies import size_line
    case = build_case("gas")
    sizing = size_line(case, lambda: solver_for("gas"), "export", [0.15, 0.2, 0.25, 0.3, 0.4],
                       max_pressure_drop_pa=3e5)
    sel = sizing.selected
    assert sel is not None and sel.passes and abs(sel.pressure_drop_pa) <= 3e5
    smaller = [c for c in sizing.candidates if c.diameter_m < sel.diameter_m]
    assert smaller and all(not c.passes for c in smaller)
    assert all(c.passes for c in sizing.candidates if c.diameter_m > sel.diameter_m)
    with pytest.raises(ValueError):
        size_line(case, lambda: solver_for("gas"), "export", [0.2])                 # no criterion
    with pytest.raises(ValueError):
        size_line(case, lambda: solver_for("gas"), "nope", [0.2], max_pressure_drop_pa=1.0)


def test_gui_scene_multipliers_reach_the_pipe(tmp_path):
    from angelica.gui.io import build_network_case_from_scene, load_scene_from_file
    scene = load_scene_from_file(ROOT / "src/angelica/tutorials/steady_black_oil/07_hot_oil_pipeline/"
                                 "07_hot_oil_pipeline.gui.json")
    comp = scene.links[0].components[0]
    scene.update_link_component_properties(scene.links[0].link_id, comp.component_id,
                                           {**comp.properties, "friction_factor_multiplier": "1.3",
                                            "holdup_multiplier": "0.9"})
    pipe = next(c for c in build_network_case_from_scene(scene).components if isinstance(c, Pipe))
    assert pipe.friction_factor_multiplier == 1.3 and pipe.holdup_multiplier == 0.9


@pytest.mark.parametrize("mode, corr", [
    ("liquid", "hazen_williams"), ("gas", "colebrook"), ("gas", "weymouth"), ("gas", "panhandle_a"),
    ("gas", "panhandle_b"), ("black_oil", "beggs_brill"), ("compositional", "colebrook")])
def test_multiplier_equals_scaling_the_friction_length(mode, corr):
    """At a fixed flow through a horizontal pipe, multiplier c equals a
    friction length L·c in every correlation and fluid mode (compressible
    and multiphase too, where ΔP is not proportional to c).  Tolerance
    1e-9: the two cases solve the same equations (review measured ≤ 3e-16)."""
    from _two_loop_networks import fluid_for, solver_for

    from angelica.closures import (
        BeggsBrillCorrelation,
        ColebrookPipeCorrelation,
        HazenWilliamsPipeCorrelation,
        PanhandleAPipeCorrelation,
        PanhandleBPipeCorrelation,
        WeymouthPipeCorrelation,
    )
    from angelica.core.case import ThermalBoundary
    correlation = {"colebrook": ColebrookPipeCorrelation, "hazen_williams": HazenWilliamsPipeCorrelation,
                   "weymouth": WeymouthPipeCorrelation, "panhandle_a": PanhandleAPipeCorrelation,
                   "panhandle_b": PanhandleBPipeCorrelation, "beggs_brill": BeggsBrillCorrelation}[corr]
    m, p_out = {"liquid": (30.0, 2e5), "gas": (10.0, 40e5), "black_oil": (10.0, 10e5),
                "compositional": (10.0, 40e5)}[mode]

    def case(c, length):
        return NetworkCase(name="p", fluid_model=fluid_for(mode), pressure_inlets=(),
                           flow_inlets=(FlowBoundary(1, m),), pressure_outlets=(PressureBoundary(2, p_out),),
                           thermal_inlets=() if mode == "liquid" else (ThermalBoundary(1, 40.0),),
                           components=(Pipe(component_id="p", start_node=1, end_node=2, diameter_m=0.2,
                                            length_m=length, absolute_roughness_m=46e-6,
                                            heat_transfer_coefficient_w_per_m2k=0.0,
                                            friction_factor_multiplier=c),))
    def solve(c):
        solver = solver_for(mode)
        getattr(solver, "_hydraulic_solver", solver).turbulent_pipe_correlation = correlation()
        return solver.solve(c).node_pressures_pa[1]
    a, b = solve(case(1.3, 4000.0)), solve(case(1.0, 5200.0))
    assert a == pytest.approx(b, rel=1e-9)


def test_line_sizing_guards():
    from _two_loop_networks import build_case, solver_for

    from angelica.studies import size_line
    case = build_case("liquid")
    s = size_line(case, lambda: solver_for("liquid"), "export", [0.1, 0.2], max_erosional_ratio=1.0)
    assert not any(c.passes for c in s.candidates)                     # single phase: not applicable
    assert "not applicable" in s.candidates[0].reasons[0]
    for kw in ({"max_pressure_drop_pa": -1.0}, {"min_node_pressure": (99, 1e5)}):
        with pytest.raises(ValueError):
            size_line(case, lambda: solver_for("liquid"), "export", [0.2], **kw)
    s = size_line(case, lambda: solver_for("liquid"), "export", [0.1, 0.3], min_node_pressure=(5, 1.5e5))
    assert [c.passes for c in s.candidates] == [c.pressure(5) >= 1.5e5 for c in (x.point for x in s.candidates)]


def test_line_sizing_a_segmented_line():
    """A line split into segments (as the GUI does) is sized as a whole."""
    from _two_loop_networks import build_case, solver_for

    from angelica.studies import size_line
    case = build_case("liquid")
    exp = next(c for c in case.components if c.component_id == "export")
    seg = (replace(exp, component_id="line_seg0", end_node=50, length_m=exp.length_m / 2),
           replace(exp, component_id="line_seg1", start_node=50, length_m=exp.length_m / 2))
    case = replace(case, components=tuple(c for c in case.components if c.component_id != "export") + seg)
    s = size_line(case, lambda: solver_for("liquid"), "line", [0.15, 0.3], max_pressure_drop_pa=1e5)
    for c in s.candidates:
        assert c.pressure_drop_pa == pytest.approx(c.point.pressure(5) - c.point.pressure(6), rel=1e-12)
