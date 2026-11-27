"""
Sensitivity studies and nodal analysis (roadmap item 7.2, v1.32.0).

Validation: CONSISTENCY only.  The nodal analysis splits the network at a
well's wellhead; its operating point must be the well's rate and wellhead
pressure in the full network solve.  Black-oil tutorial 12: well 1
6.6933 against 6.6932 kg/s (1.5e-5) and 4.04543 against 4.04548 MPa;
well 2 2.2189 against 2.2190 kg/s.  Tolerance 1e-3: the bisection stops at
1e-4 of the rate and both solves meet the hydraulic criterion (1e-3 per
node).  A dead well has no crossing on the stable branch.  No published
nodal-analysis example was read (Economides et al. and Beggs are not in the
repository).
"""
from __future__ import annotations

import importlib.util
import sys
from functools import lru_cache
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT / "src", ROOT / "tests"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from angelica.studies import nodal_analysis, sensitivity, sensitivity_2d, set_parameter  # noqa: E402


def tutorial(path: str):
    spec = importlib.util.spec_from_file_location(path.replace("/", "_"), ROOT / "src/angelica/tutorials" / path / "run.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@lru_cache(maxsize=None)
def t12():
    from angelica.closures import BeggsBrillCorrelation
    from angelica.solvers import SteadyBlackOilSolver
    mod = tutorial("steady_black_oil/12_wells_with_tubing")
    case = mod.build_case()

    def factory():
        return SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation())
    return mod, case, factory, factory().solve(case)


def test_set_parameter_targets():
    from _two_loop_networks import build_case
    case = build_case("liquid")
    assert set_parameter(case, "pressure:6", 2e5).pressure_outlets[0].pressure_pa == 2e5
    new = set_parameter(case, "component:2-3.diameter_m", 0.3)
    assert next(c for c in new.components if c.component_id == "2-3").diameter_m == 0.3
    assert next(c for c in case.components if c.component_id == "2-3").diameter_m != 0.3     # a copy
    for bad in ("pressure:99", "component:nope.diameter_m", "component:2-3.colour", "speed:1"):
        with pytest.raises(ValueError):
            set_parameter(case, bad, 1.0)


def test_sensitivity_one_and_two_parameters():
    from _two_loop_networks import build_case, solver_for
    case = build_case("liquid")
    pts = sensitivity(case, lambda: solver_for("liquid"), "pressure:6", [0.5e5, 1.0e5, 2.0e5, 3.0e5])
    flows = [p.flow("export") for p in pts]
    assert all(p.converged for p in pts)
    assert all(b < a for a, b in zip(flows, flows[1:]))           # higher outlet pressure, less flow
    # each point is the case solved with that value
    direct = solver_for("liquid").solve(set_parameter(case, "pressure:6", 2.0e5))
    assert pts[2].pressure(5) == pytest.approx(direct.node_pressures_pa[5], rel=1e-12)
    grid = sensitivity_2d(case, lambda: solver_for("liquid"), "pressure:6", [1e5, 2e5],
                          "component:3-5.diameter_m", [0.1, 0.2])
    assert [g.values for g in grid] == [(1e5, 0.1), (1e5, 0.2), (2e5, 0.1), (2e5, 0.2)]
    assert grid[1].flow("3-5") > grid[0].flow("3-5")             # wider pipe carries more


def test_sensitivity_keeps_non_converged_points():
    """A study point that does not converge is kept with its diagnosis."""
    from _two_loop_networks import build_case, solver_for
    pts = sensitivity(build_case("black_oil"), lambda: solver_for("black_oil", 5.0), "temperature:1", [70.0])
    assert len(pts) == 1 and pts[0].result is not None


@pytest.mark.parametrize("i", [0, 1])
def test_nodal_operating_point_is_the_network_solution(i):
    mod, case, factory, full = t12()
    well = [c for c in case.components if c.component_id.startswith("well")][i]
    ws = next(w for w in full.well_states if w.component_id == well.component_id)
    na = nodal_analysis(case, factory, well.component_id)
    assert na.operating_point is not None, na.note
    m, p = na.operating_point
    assert m == pytest.approx(ws.mass_flow_kg_per_s, rel=1e-3)
    assert p == pytest.approx(full.node_pressures_pa[well.end_node], rel=1e-3)
    # the curves: inflow decreasing on its stable branch, outflow increasing
    out = [q for _m, q in na.outflow if q is not None]
    assert all(b >= a for a, b in zip(out, out[1:]))


def test_nodal_dead_well_has_no_crossing():
    mod, case, factory, full = t12()
    dead = next(w for w in full.well_states if w.state != "producing")
    na = nodal_analysis(case, factory, dead.component_id, n_rates=4)
    assert na.operating_point is None and na.note


def test_nodal_rejects_what_it_cannot_split():
    mod, case, factory, _full = t12()
    with pytest.raises(ValueError):
        nodal_analysis(case, factory, "trunk")
    t13 = tutorial("steady_black_oil/13_gas_lift")
    with pytest.raises(ValueError, match="gas lift"):
        nodal_analysis(t13.build_case(60_000.0), factory, "well")


def test_gui_study_lists():
    pytest.importorskip("sv_ttk")
    from _two_loop_networks import build_case, solver_for

    from angelica.gui.studies import output_value, study_outputs, study_targets
    case = build_case("liquid")
    targets = dict(study_targets(case))
    assert targets["Pressure of node 6 (Pa)"] == "pressure:6"
    assert "component:2-3.diameter_m" in targets.values()
    outputs = dict(study_outputs(case))
    pt = sensitivity(case, lambda: solver_for("liquid"), "pressure:6", [1e5])[0]
    assert output_value(pt, outputs["Mass flow through export (kg/s)"]) == pytest.approx(pt.flow("export"))


def test_nodal_compositional_gas_condensate_well():
    """Compositional tutorial 09: obtained 6.2268 = 6.2268 kg/s, 6.97104
    against 6.97114 MPa (1.4e-5)."""
    from angelica.solvers import SteadyCompositionalSolver
    t09 = tutorial("steady_compositional/09_gas_condensate_wells")
    case = t09.build_case()
    full = SteadyCompositionalSolver().solve(case)
    ws = next(w for w in full.well_states if w.component_id == "well_Condensate")
    na = nodal_analysis(case, SteadyCompositionalSolver, "well_Condensate", n_rates=8)
    m, p = na.operating_point
    assert m == pytest.approx(ws.mass_flow_kg_per_s, rel=1e-3)
    assert p == pytest.approx(ws.wellhead_pressure_pa, rel=1e-3)


# Cases found in review (v1.32.0): before the fixes nodal_analysis returned
# a wrong operating point for each of them without a word.
def _check_against_network(case, factory, well_id):
    full = factory().solve(case)
    ws = next(w for w in full.well_states if w.component_id == well_id)
    na = nodal_analysis(case, factory, well_id, n_rates=8)
    assert na.operating_point is not None, na.note
    assert na.operating_point[0] == pytest.approx(ws.mass_flow_kg_per_s, rel=1e-3)
    assert na.operating_point[1] == pytest.approx(ws.wellhead_pressure_pa, rel=1e-3)


def test_nodal_keeps_the_case_pvt_correlation():
    """Vasquez-Beggs: the well's fluid inherits the case fluid's correlation
    (review: +3.9 % on rate with Standing in its place)."""
    from dataclasses import replace
    mod, case, factory, _ = t12()
    case = replace(case, fluid_model=replace(case.fluid_model, correlation="VASQUEZ_BEGGS"))
    _check_against_network(case, factory, "well_Well1")


def test_nodal_with_an_isothermal_solver():
    from angelica.closures import BeggsBrillCorrelation
    from angelica.solvers import SteadyBlackOilSolver
    from angelica.solvers.steady_black_oil import BlackOilSolverSettings
    mod, case, _f, _ = t12()

    def factory():
        return SteadyBlackOilSolver(black_oil_settings=BlackOilSolverSettings(isothermal_temperature_c=40.0),
                                    turbulent_pipe_correlation=BeggsBrillCorrelation())
    _check_against_network(case, factory, "well_Well1")


def test_nodal_with_a_shared_reservoir_node():
    """Well 2 drawing from Well 1's reservoir node: the node keeps its fluid
    in the outflow network (review: +2.7 % on rate)."""
    from dataclasses import replace
    mod, case, factory, _ = t12()
    comps = tuple(replace(c, start_node=10) if c.component_id == "well_Well2" else c for c in case.components)
    case = replace(case, components=comps,
                   pressure_inlets=tuple(b for b in case.pressure_inlets if b.node_id != 20),
                   thermal_inlets=tuple(t for t in case.thermal_inlets if t.node_id != 20),
                   inlet_fluid_bcs=tuple(b for b in case.inlet_fluid_bcs if b.node_id != 20))
    _check_against_network(case, factory, "well_Well1")


def test_nodal_refuses_a_shared_wellhead():
    from dataclasses import replace
    mod, case, factory, _ = t12()
    comps = tuple(replace(c, end_node=11) if c.component_id == "well_Well2" else c for c in case.components)
    case = replace(case, components=tuple(c for c in comps if c.component_id != "flowline_Well2"))
    with pytest.raises(ValueError, match="also receives"):
        nodal_analysis(case, factory, "well_Well1")


def test_sensitivity_keeps_points_that_raise():
    from _two_loop_networks import build_case, solver_for
    pts = sensitivity(build_case("liquid"), lambda: solver_for("liquid"), "component:2-3.diameter_m", [0.2, -0.1])
    assert pts[0].converged and not pts[1].converged and pts[1].result is None and pts[1].error
    with pytest.raises(ValueError):
        sensitivity(build_case("liquid"), lambda: solver_for("liquid"), "component:2-3.area_m2", [0.1])
