"""
Looped networks in every fluid mode, thermal and isothermal (roadmap goal,
part 1: what Angelica already does, checked with numbers).

The network (``tests/_two_loop_networks.py``) has two independent loops and
two inlets.  Fluids: water (isothermal and thermal), methane gas
(Peng-Robinson), black-oil crude with gas and water (WOR 0.4), and a
seven-component condensate with water that forms free water.

Checks, all by internal consistency (no published looped network for
these fluids was available; the Hanoi water network of tutorial 07 is the
published one — EPANET results, ``test_tutorial_suite.py``):

1. Mass balance at every junction.  Tolerance: the hydraulic solver's own
   stopping criterion ``nodal_mass_imbalance_rel_tolerance`` = 1e-3, so
   this check is close to automatic for a converged solve; it guards the
   bookkeeping of flows into the result.  Obtained: ≤ 1.07e-5.
2. Kirchhoff's loop law holds by construction in a solver with one
   pressure per node: around any loop the node-pressure differences sum to
   zero exactly.  What can be wrong is whether each link's pressure drop is
   the one its correlation gives at its flow.  So each pipe is solved
   alone — its upstream node's pressure, temperature and composition from
   the network, the network's flow imposed — and must give the network's
   downstream pressure.  Tolerance 1e-4 of the pipe's pressure drop, the
   velocity loop's tolerance (``velocity_loop_tolerance``; ΔP ∝ v² near
   2·1e-4) and the outer loops' density criterion 1e-4.  Obtained:
   ≤ 3.02e-5 (black-oil thermal; compositional isothermal 3.00e-5),
   ≤ 6e-8 for the liquids.  (v1.28.0 also summed these drops around each
   loop; the network pressures cancel in that sum, so it only repeated
   this check and was removed in review.)
3. Thermal runs: the global energy balance closes within 0.01 % (the
   roadmap's figure).  Obtained: ≤ 2.77e-4 %.
4. Isothermal runs (``isothermal_temperature_c``, gas, black-oil and
   compositional): every node at that temperature, no energy balance.
5. Hydrate margin per pipe in black-oil (gas-gravity chart) and
   compositional (vdW-P): equal to the stand-alone function at the pipe's
   colder end.
"""
from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT / "src", ROOT / "tests"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from _two_loop_networks import (  # noqa: E402
    COMPONENTS,
    LINKS,
    build_case,
    single_pipe_case,
    solver_for,
)

THERMAL = ["liquid", "liquid_thermal", "gas", "black_oil", "compositional"]
ISOTHERMAL = [("gas", 30.0), ("black_oil", 60.0), ("compositional", 60.0)]
RUNS = [(m, None) for m in THERMAL] + ISOTHERMAL
JUNCTIONS = (2, 3, 4, 5)
REL_TOL = 1e-4


@lru_cache(maxsize=None)
def solved(mode: str, iso: float | None):
    case = build_case(mode)
    return case, solver_for(mode, iso).solve(case)


def flows_of(case, result) -> dict[str, float]:
    return {c.component_id: f.mass_flow_kg_per_s for c, f in zip(case.components, result.component_flows)}


@lru_cache(maxsize=None)
def pipe_errors(mode: str, iso: float | None) -> dict[str, float]:
    """Per pipe: relative error of the stand-alone solve's downstream
    pressure against the network's, as a fraction of the pipe's drop."""
    case, res = solved(mode, iso)
    out = {}
    for cid, *_ in LINKS:
        single, up, down, p_down = single_pipe_case(mode, case, res, cid)
        alone = solver_for(mode, iso).solve(single)
        assert alone.converged, cid
        p_alone = alone.node_pressures_pa[down]
        out[cid] = (p_alone - p_down) / (res.node_pressures_pa[up] - p_down)
    return out


@pytest.mark.parametrize("mode, iso", RUNS)
def test_converged_with_mass_balance_at_every_junction(mode, iso):
    case, res = solved(mode, iso)
    assert res.converged
    flows = flows_of(case, res)
    net = dict.fromkeys(JUNCTIONS, 0.0)
    for c in case.components:
        if c.start_node in net:
            net[c.start_node] -= flows[c.component_id]
        if c.end_node in net:
            net[c.end_node] += flows[c.component_id]
    scale = max(abs(v) for v in flows.values())
    assert max(abs(v) for v in net.values()) / scale < 1e-3
    # every loop pipe carries flow (no dead branch hiding a loop)
    assert min(abs(flows[c]) for c in ("2-3", "2-4", "3-4", "3-5", "4-5")) > 0.05 * scale


@pytest.mark.parametrize("mode, iso", RUNS)
def test_each_pipe_alone_reproduces_the_network_pressure(mode, iso):
    for cid, err in pipe_errors(mode, iso).items():
        assert abs(err) < REL_TOL, f"{cid}: {err:.2e}"


@pytest.mark.parametrize("mode", [m for m in THERMAL if m != "liquid"])
def test_thermal_runs_close_the_energy_balance(mode):
    case, res = solved(mode, None)
    assert abs(res.global_energy_balance.energy_error_pct) < 0.01
    temps = res.node_temperatures_c
    assert temps[6] < max(temps[1], temps[7])        # the network loses heat to 15 °C ambient


@pytest.mark.parametrize("mode, iso", ISOTHERMAL)
def test_isothermal_runs_hold_every_node(mode, iso):
    case, res = solved(mode, iso)
    assert all(t == pytest.approx(iso, abs=1e-12) for t in res.node_temperatures_c.values())
    assert res.global_energy_balance is None
    assert any("isothermal" in w for w in res.warnings)   # the inlet temperatures differ from iso
    _, thermal = solved(mode, None)
    assert res.node_pressures_pa[5] != pytest.approx(thermal.node_pressures_pa[5], rel=1e-6)


def test_three_phase_fluids_are_three_phase():
    """Oil, gas and water present where the mode allows it."""
    case, res = solved("black_oil", None)
    holdups = [f.multiphase.liquid_holdup for f in res.component_flows if f.multiphase]
    assert holdups and all(0.0 < h < 1.0 for h in holdups)        # gas and liquid
    assert all(f.free_water_fraction > 0.0 for f in res.component_flows)   # water (WOR 0.4)
    case, res = solved("compositional", None)
    water = COMPONENTS.index("water")
    assert any(f.free_water_fraction > 0.0 for f in res.component_flows)
    assert any(f.liquid_fraction and f.liquid_fraction > 0.0 for f in res.component_flows)  # condensate
    assert all(f.zs[water] > 0.0 for f in res.component_flows)


@pytest.mark.parametrize("mode", ["black_oil", "compositional"])
def test_hydrate_margin_per_pipe(mode):
    from angelica.properties.hydrate import hydrate_formation_temperature, katz_hydrate_temperature
    case, res = solved(mode, None)
    P, T = res.node_pressures_pa, res.node_temperatures_c
    fluid = case.fluid_model
    for comp, f in zip(case.components, res.component_flows):
        assert f.hydrate_temperature_c is not None, comp.component_id
        worst = None
        for n in (comp.start_node, comp.end_node):
            if mode == "black_oil":
                pt = katz_hydrate_temperature(fluid.gas_gravity, P[n])
            else:
                pt = hydrate_formation_temperature(fluid.component_names, f.zs, P[n], fluid.eos_name)
            sub = pt.temperature_c - T[n]
            if worst is None or sub > worst[0]:
                worst = (sub, pt.temperature_c)
        assert f.hydrate_subcooling_k == pytest.approx(worst[0], abs=1e-6)
        assert f.hydrate_temperature_c == pytest.approx(worst[1], abs=1e-6)
