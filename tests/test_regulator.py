"""
Pressure regulator (``PressureRegulator``, roadmap item 5.4), the
non-convergence diagnosis (``SolveResult.diagnosis``) and two black-oil
fixes found on the way (v1.27.0).

Standing of the checks
----------------------
- Controlling: the controlled node is at the set point (exact: a fixed-
  pressure row) and the flow equals the flow of the network on the other
  side solved alone with that pressure imposed — an independent solve
  (obtained 1.2e-7; tolerance: the solver's nodal tolerance).
- Open: identical to the same network with a check valve of the same K
  (obtained 1e-7; nodal tolerance).
- Closed: zero flow when the outlet is at or above the set point
  (reducing) or the inlet below it (back-pressure), and with reversed
  pressures.
- Isenthalpic: the gas leaving a reducing regulator has the inlet's
  enthalpy at the outlet pressure — checked with the fluid's own h(P, T),
  tolerance cp·0.01 K, the solver's temperature tolerance.
- Sweeps: every point of an outlet-pressure and an inlet-pressure sweep
  converges and the flow is monotonic (it may vary by 1e-5 relative where
  the outlet temperature changes with the throttling: tolerance 1e-3).
- Not compared with a published regulator test or a manufacturer's curve:
  the wide-open valve is a K-factor fitting (K = 5 by default, an
  assumption) and there is no droop or lock-up.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures import BeggsBrillCorrelation  # noqa: E402
from angelica.core.case import (  # noqa: E402
    FlowBoundary,
    InletCompositionBC,
    InletFluidBC,
    NetworkCase,
    PressureBoundary,
    ThermalBoundary,
)
from angelica.core.components import CheckValve, Pipe, PressureRegulator  # noqa: E402
from angelica.properties.black_oil import BlackOilFluid, compute_pvt, temperature_from_enthalpy  # noqa: E402
from angelica.properties.compositional_fluid import CompositionalFluid  # noqa: E402
from angelica.solvers import (  # noqa: E402
    SteadyBlackOilSolver,
    SteadyCompositionalSolver,
    SteadyIsothermalIncompressibleSolver,
)

NODAL = 2e-3      # the solver's nodal tolerance (1e-3 of inflow + outflow)
CRUDE = {"api_gravity": 35.0, "gas_gravity": 0.70, "gor_sc_m3_per_m3": 120.0, "wor_sc_m3_per_m3": 0.3}
GAS = ("methane", "ethane", "propane", "n-butane", "nitrogen", "carbon dioxide")
GAS_ZS = (0.90, 0.05, 0.02, 0.01, 0.01, 0.01)


def _pipe(cid, a, b, L, D=0.1524, U=2.0, T_amb=15.0):
    return Pipe(component_id=cid, start_node=a, end_node=b, diameter_m=D, length_m=L, absolute_roughness_m=4.6e-5,
                heat_transfer_coefficient_w_per_m2k=U, ambient_temperature_c=T_amb)


def _line(p_in, p_out, device):
    """Inlet 1 → 2 km pipe → device 2→3 → 3 km pipe → outlet 4 (black-oil)."""
    return NetworkCase(
        name="regulator", fluid_model=BlackOilFluid(**CRUDE),
        pressure_inlets=(PressureBoundary(1, p_in),), pressure_outlets=(PressureBoundary(4, p_out),),
        thermal_inlets=(ThermalBoundary(1, 60.0),), inlet_fluid_bcs=(InletFluidBC(1, **CRUDE),),
        components=(_pipe("up", 1, 2, 2000.0), device, _pipe("down", 3, 4, 3000.0)))


def _reg(p_set=3.0e6, mode="downstream"):
    return PressureRegulator(component_id="reg", start_node=2, end_node=3, diameter_m=0.1524,
                             set_pressure_pa=p_set, mode=mode)


def _flows(res):
    return {c.label.split(":")[1]: c.mass_flow_kg_per_s for c in res.component_flows}


# ── reducing regulator ────────────────────────────────────────────────────────

def test_reducing_regulator_holds_its_outlet():
    res = SteadyBlackOilSolver().solve(_line(6.0e6, 1.0e6, _reg()))
    assert res.converged
    reg = res.regulators[0]
    assert reg.state == "controlling" and "holding the outlet" in reg.message
    assert reg.outlet_pressure_pa == 3.0e6
    # the flow is what the downstream pipe passes from 3 MPa, solved alone
    alone = NetworkCase(
        name="down", fluid_model=BlackOilFluid(**CRUDE), pressure_inlets=(PressureBoundary(3, 3.0e6),),
        pressure_outlets=(PressureBoundary(4, 1.0e6),),
        thermal_inlets=(ThermalBoundary(3, res.node_temperatures_c[3]),),
        inlet_fluid_bcs=(InletFluidBC(3, **CRUDE),), components=(_pipe("down", 3, 4, 3000.0),))
    ref = SteadyBlackOilSolver().solve(alone).component_flows[0].mass_flow_kg_per_s
    assert reg.mass_flow_kg_per_s == pytest.approx(ref, rel=NODAL)


@pytest.mark.parametrize("p_in", [3.5e6, 2.5e6])
def test_reducing_regulator_wide_open_is_a_check_valve(p_in):
    """Upstream too low for the set point: wide open, i.e. the network with
    a check valve of the same K."""
    a = SteadyBlackOilSolver().solve(_line(p_in, 1.0e6, _reg()))
    b = SteadyBlackOilSolver().solve(_line(p_in, 1.0e6, CheckValve(component_id="reg", start_node=2, end_node=3,
                                                                     diameter_m=0.1524, loss_coefficient=5.0)))
    assert a.converged and b.converged
    reg = a.regulators[0]
    assert reg.state == "open" and reg.outlet_pressure_pa < 3.0e6
    assert reg.mass_flow_kg_per_s == pytest.approx(_flows(b)["reg"], rel=NODAL)
    assert reg.outlet_pressure_pa == pytest.approx(b.node_pressures_pa[3], rel=NODAL)


@pytest.mark.parametrize("p_in, p_out", [(6.0e6, 3.0e6), (6.0e6, 5.0e6), (6.0e6, 7.0e6), (1.0e6, 3.0e6)])
def test_reducing_regulator_closed(p_in, p_out):
    """Outlet at or above the set point, or pressures reversed: no flow."""
    res = SteadyBlackOilSolver().solve(_line(p_in, p_out, _reg()))
    assert res.converged
    assert res.regulators[0].state == "closed"
    assert all(abs(q) < 1e-5 for q in _flows(res).values())


# ── back-pressure regulator ───────────────────────────────────────────────────

@pytest.mark.parametrize("p_in, state", [(4.1e6, "controlling"), (3.8e6, "closed"), (6.0e6, "open")])
def test_back_pressure_regulator_states(p_in, state):
    res = SteadyBlackOilSolver().solve(_line(p_in, 1.0e6, _reg(4.0e6, "upstream")))
    assert res.converged
    reg = res.regulators[0]
    assert reg.state == state
    if state == "controlling":
        assert reg.inlet_pressure_pa == 4.0e6
        fl = _flows(res)
        assert fl["up"] == pytest.approx(fl["reg"], rel=NODAL) and fl["down"] == pytest.approx(fl["reg"], rel=NODAL)
    elif state == "closed":
        assert abs(reg.mass_flow_kg_per_s) < 1e-5
    else:
        assert reg.inlet_pressure_pa > 4.0e6      # downstream too high to hold the inlet down to the set point


# ── sweeps (the standing rule: beyond the tutorial's points) ──────────────────

@pytest.mark.parametrize("mode, var, values, fixed", [
    ("downstream", "out", np.arange(0.5e6, 5.6e6, 0.5e6), 6.0e6),
    ("downstream", "in", np.arange(1.5e6, 8.1e6, 0.5e6), 1.0e6),
    ("upstream", "in", np.arange(1.5e6, 8.1e6, 0.5e6), 1.0e6),
])
def test_regulator_sweeps_converge_and_are_monotonic(mode, var, values, fixed):
    """Every point converges, crossing the closed / controlling / open
    boundaries; the flow falls with the outlet pressure and rises with the
    inlet pressure (1e-3: at fixed outlet the throttled fluid's temperature
    moves the flow by ~1e-5)."""
    p_set = 3.0e6 if mode == "downstream" else 4.0e6
    prev, states = None, set()
    for v in values:
        p_in, p_out = (v, fixed) if var == "in" else (fixed, v)
        res = SteadyBlackOilSolver().solve(_line(float(p_in), float(p_out), _reg(p_set, mode)))
        assert res.converged, (v, res.diagnosis.message)
        m = res.regulators[0].mass_flow_kg_per_s
        states.add(res.regulators[0].state)
        if prev is not None:
            if var == "out":
                assert m <= prev * (1 + 1e-3) + 1e-6
            else:
                assert m >= prev * (1 - 1e-3) - 1e-6
        prev = m
    assert len(states) >= 2


# ── compositional ─────────────────────────────────────────────────────────────

def _city_gate(p_in, p_out, p_set=20e5):
    return NetworkCase(
        name="city gate", fluid_model=CompositionalFluid(GAS, GAS_ZS),
        pressure_inlets=(PressureBoundary(1, p_in),), pressure_outlets=(PressureBoundary(4, p_out),),
        thermal_inlets=(ThermalBoundary(1, 20.0),), inlet_composition_bcs=(InletCompositionBC(1, GAS_ZS),),
        components=(_pipe("transmission", 1, 2, 20000.0, D=0.2032, T_amb=10.0),
                    PressureRegulator(component_id="reg", start_node=2, end_node=3, diameter_m=0.2032,
                                      set_pressure_pa=p_set),
                    _pipe("distribution", 3, 4, 2000.0, T_amb=10.0)))


def test_compositional_reducing_regulator_is_isenthalpic():
    """A 64 → 20 bar let-down cools the gas by ~27 K (Joule-Thomson): the
    outlet state has the inlet's enthalpy (fluid's own h(P, T))."""
    res = SteadyCompositionalSolver().solve(_city_gate(70e5, 5e5))
    assert res.converged
    reg = res.regulators[0]
    assert reg.state == "controlling" and reg.outlet_pressure_pa == 20e5
    from types import SimpleNamespace
    fluid = CompositionalFluid(GAS, GAS_ZS)
    link = SimpleNamespace(zs=GAS_ZS)
    T2, T3 = res.node_temperatures_c[2], res.node_temperatures_c[3]
    assert T3 < T2 - 20.0
    h_in = fluid.enthalpy_for_link(link, reg.inlet_pressure_pa, T2)
    h_out = fluid.enthalpy_for_link(link, reg.outlet_pressure_pa, T3)
    cp = fluid.specific_heat_at(link, reg.outlet_pressure_pa, T3)
    assert abs(h_out - h_in) < cp * 0.01


@pytest.mark.parametrize("p_in, state", [(70e5, "controlling"), (25e5, "open"), (15e5, "open")])
def test_compositional_supply_too_weak_opens_the_regulator(p_in, state):
    """With 25 bar at the inlet the 20 km line cannot deliver the demand at
    20 bar downstream: the regulator is wide open (before v1.27.0's limit on
    the controlled flow, the upstream pressure went to zero and the flash
    raised)."""
    res = SteadyCompositionalSolver().solve(_city_gate(p_in, 5e5))
    assert res.converged
    assert res.regulators[0].state == state


# ── regulator on a separator's gas line ───────────────────────────────────────

def _tutorial(path: str):
    import importlib.util
    spec = importlib.util.spec_from_file_location("tut", SRC_ROOT / "angelica" / "tutorials" / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("p_export, state", [(0.5e6, "controlling"), (2.0e6, "controlling"), (5.0e6, "open")])
def test_back_pressure_regulator_holds_a_separator(p_export, state):
    """Tutorial 14's separator with a back-pressure regulator (3 MPa) on its
    gas line: the separator stays at 3 MPa whatever the export pressure
    below it; above, the regulator is open.  The separator's liquid outlets
    take part of the feed, which the regulator's implicit flow accounts for
    (without it the hydraulics diverged)."""
    t14 = _tutorial("steady_black_oil/14_separator/run.py")
    case = t14.build_case(p_export=p_export)
    comps = []
    for k in case.components:
        if k.component_id == "gas_line":
            comps.append(PressureRegulator(component_id="bpr", start_node=t14.SEPARATOR, end_node=200,
                                           diameter_m=0.1524, set_pressure_pa=3.0e6, mode="upstream"))
            comps.append(type(k)(**{**k.__dict__, "start_node": 200}))
        else:
            comps.append(k)
    case = NetworkCase(**{**case.__dict__, "components": tuple(comps)})
    res = SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation()).solve(case)
    assert res.converged
    reg, sep = res.regulators[0], res.separators[0]
    assert reg.state == state
    if state == "controlling":
        assert sep.pressure_pa == 3.0e6
    else:
        assert sep.pressure_pa > p_export
    assert reg.mass_flow_kg_per_s == pytest.approx(sep.gas_mass_flow_kg_per_s, rel=NODAL)


# ── layout errors ─────────────────────────────────────────────────────────────

def test_regulator_layout_errors():
    with pytest.raises(ValueError, match="set_pressure_pa"):
        PressureRegulator(component_id="r", start_node=1, end_node=2, diameter_m=0.1)
    with pytest.raises(ValueError, match="mode"):
        PressureRegulator(component_id="r", start_node=1, end_node=2, diameter_m=0.1, set_pressure_pa=1e6,
                          mode="both")
    # controlled node is a boundary
    case = NetworkCase(
        name="bad", fluid_model=BlackOilFluid(**CRUDE), pressure_inlets=(PressureBoundary(1, 6e6),),
        pressure_outlets=(PressureBoundary(3, 1e6),), thermal_inlets=(ThermalBoundary(1, 60.0),),
        components=(_pipe("up", 1, 2, 1000.0),
                    PressureRegulator(component_id="r", start_node=2, end_node=3, diameter_m=0.1,
                                      set_pressure_pa=3e6)))
    with pytest.raises(ValueError, match="boundary"):
        SteadyBlackOilSolver().solve(case)
    # single-phase solvers
    with pytest.raises(ValueError, match="PressureRegulator"):
        SteadyIsothermalIncompressibleSolver().solve(_line(6e6, 1e6, _reg()))


# ── diagnosis of a solve that does not converge ───────────────────────────────

def test_diagnosis_names_the_node_and_the_flow_pattern():
    """Black-oil tutorial 14 at 0.5 MPa export: the oil line cannot pass
    the flow its LCV imposes (Beggs-Brill pattern boundary).  The diagnosis
    names the stuck node (103, the oil LCV outlet) and the oil line's
    alternating patterns."""
    t14 = _tutorial("steady_black_oil/14_separator/run.py")
    res = SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation()).solve(
        t14.build_case(p_export=0.5e6))
    assert not res.converged
    d = res.diagnosis
    assert "hydraulics" in d.failed
    assert d.worst_node_id == t14.OIL_LCV_OUT
    oil = [s for s in d.pattern_switches if s.component_id == "oil_line"]
    assert oil and set(oil[0].patterns) >= {"intermittent", "distributed"}
    assert "oil_line" in d.message


def test_diagnosis_of_a_single_pipe_at_an_impossible_flow():
    """The oil line alone at 7.12 kg/s: Beggs-Brill gives no solution
    between 6.87 and 7.55 kg/s (intermittent up to 6.86 at 0.799 MPa inlet,
    distributed from 7.56 at 0.730 MPa).  The pattern alternates within the
    hydraulic iterations; the diagnosis still names it."""
    fluid = BlackOilFluid(api_gravity=37.98, gas_gravity=0.696, gor_sc_m3_per_m3=3.134, wor_sc_m3_per_m3=0.0)
    case = NetworkCase(name="oil line", fluid_model=fluid, pressure_inlets=(),
                       flow_inlets=(FlowBoundary(1, 7.12),), pressure_outlets=(PressureBoundary(2, 0.2e6),),
                       thermal_inlets=(ThermalBoundary(1, 36.6),),
                       components=(_pipe("oil_line", 1, 2, 3000.0, D=0.1016, U=3.0),))
    res = SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation()).solve(case)
    assert not res.converged
    assert any(s.component_id == "oil_line" for s in res.diagnosis.pattern_switches)


def test_converged_solve_has_no_diagnosis():
    res = SteadyBlackOilSolver().solve(_line(6.0e6, 1.0e6, _reg()))
    assert res.converged and res.diagnosis is None


# ── black-oil fixes found on the way ──────────────────────────────────────────

def test_black_oil_enthalpy_inversion_returns_the_root():
    """Fixed in v1.27.0: plain regula falsi kept one end of the bracket and
    returned the bracket's midpoint after 60 steps — 61.64 °C instead of
    59.15 °C here, alternating between iterations (a check valve's outlet
    temperature cycled by 2.5 K).  Values from that case."""
    fluid = BlackOilFluid(**CRUDE)
    P_out = 4971526.188676454
    for h in (148821.56615580758, 148821.8255188298):
        T = fluid.temperature_from_enthalpy_for_link(None, h, P_out, 59.14294757173744)
        # the root to the resolution of h (T rounded to 1e-4 K: 0.65 J/kg)
        assert abs(fluid.enthalpy_for_link(None, P_out, T) - h) < 1.0
        assert T == pytest.approx(59.145, abs=0.01)
    # a plain monotonic function, far from the guess
    root = (-3.0 + (9.0 + 4.0 * 0.001 * 300.0) ** 0.5) / (2.0 * 0.001)     # 96.872
    assert temperature_from_enthalpy(lambda T: 3.0 * T + 0.001 * T * T, 300.0, 10.0) == pytest.approx(
        root, abs=1e-3)


def test_check_valve_case_that_cycled_now_converges():
    """The same network with a K = 5 check valve at 3.05 MPa outlet ended
    with converged=False (the cycle above); a plain fitting converged."""
    res = SteadyBlackOilSolver().solve(_line(6.0e6, 3.05e6, CheckValve(
        component_id="reg", start_node=2, end_node=3, diameter_m=0.1524, loss_coefficient=5.0)))
    assert res.converged


def test_black_oil_properties_below_their_range_do_not_raise():
    """Below −6.7 °C (a numerical floor, not the correlation's limit: see
    test_black_oil_range.py) properties are taken at −6.7 °C instead of
    raising (at −17.8 °C the Beggs-Robinson log of T in °F fails; a few °F
    above it 10^(10^z) overflows for heavy oils)."""
    args = (5e6, 35.0, 0.7, 120.0, 0.3)
    low = compute_pvt(args[0], -30.0, *args[1:])
    at_limit = compute_pvt(args[0], -6.7, *args[1:])
    assert low.viscosity_oil_pa_s == at_limit.viscosity_oil_pa_s
    compute_pvt(args[0], 3000.0, *args[1:])      # an unsettled energy iteration: no overflow


# ── two regulators in series (v1.27.1) ────────────────────────────────────────
#
# Each case is checked for consistency of every regulator's state with its
# pressures (controlling: the held node at the set point; open: the set
# point out of reach; closed: no flow and the set point already met) and for
# one flow through the series.  Tolerance: NODAL of the set point and of the
# flow.  Where two regulators have contradictory set points the solution is
# no flow; the line between them then has no pressure of its own and the
# result says so (a warning).  Off-suite, the full grids (black-oil 9 × 9
# set points, compositional 8 × 8, five layouts: 725 solves) all converged
# and were consistent.

_PAIR_MODES = {"bpr+reducer": ("upstream", "downstream"), "reducer+bpr": ("downstream", "upstream"),
               "reducer+reducer": ("downstream", "downstream")}


def _pair_case(layout, p1, p2, compositional=False, direct=False):
    m1, m2 = _PAIR_MODES[layout]
    D = 0.2032 if compositional else 0.1524

    def reg(cid, a, b, p, mode):
        return PressureRegulator(component_id=cid, start_node=a, end_node=b, diameter_m=D,
                                 set_pressure_pa=p, mode=mode)
    comps = [_pipe("up", 1, 2, 5000.0, D=D, U=3.0, T_amb=10.0 if compositional else 15.0), reg("r1", 2, 3, p1, m1)]
    comps += [reg("r2", 3, 5, p2, m2)] if direct else [
        _pipe("mid", 3, 4, 1000.0 if compositional else 500.0, D=D, U=3.0), reg("r2", 4, 5, p2, m2)]
    comps += [_pipe("down", 5, 9, 2000.0, D=D, U=3.0)]
    if compositional:
        return NetworkCase(name=layout, fluid_model=CompositionalFluid(GAS, GAS_ZS),
                           pressure_inlets=(PressureBoundary(1, 70e5),), pressure_outlets=(PressureBoundary(9, 5e5),),
                           thermal_inlets=(ThermalBoundary(1, 20.0),),
                           inlet_composition_bcs=(InletCompositionBC(1, GAS_ZS),), components=tuple(comps))
    return NetworkCase(name=layout, fluid_model=BlackOilFluid(**CRUDE),
                       pressure_inlets=(PressureBoundary(1, 8e6),), pressure_outlets=(PressureBoundary(9, 0.5e6),),
                       thermal_inlets=(ThermalBoundary(1, 60.0),), inlet_fluid_bcs=(InletFluidBC(1, **CRUDE),),
                       components=tuple(comps))


def _assert_consistent(res):
    assert res.converged, res.diagnosis.message
    flows = [c.mass_flow_kg_per_s for c in res.component_flows]
    top = max(max(flows), 1e-3)
    assert max(flows) - min(flows) <= NODAL * top + 1e-3          # one flow through the series
    if top <= 1e-3:
        # no flow: the line between two closed regulators is cut off
        assert all(g.state == "closed" or abs(g.mass_flow_kg_per_s) <= 1e-3 for g in res.regulators)
        return
    for g in res.regulators:
        tol = NODAL * g.set_pressure_pa
        held = g.outlet_pressure_pa if g.mode == "downstream" else g.inlet_pressure_pa
        if g.state == "controlling":
            assert held == pytest.approx(g.set_pressure_pa, abs=tol)
        elif g.state == "open":
            assert (g.outlet_pressure_pa <= g.set_pressure_pa + tol if g.mode == "downstream"
                    else g.inlet_pressure_pa >= g.set_pressure_pa - tol)
        else:
            assert abs(g.mass_flow_kg_per_s) <= 1e-3


@pytest.mark.parametrize("layout", list(_PAIR_MODES))
def test_black_oil_regulator_pairs_set_point_sweep(layout):
    for p1 in (1e6, 3e6, 5e6, 7e6, 7.8e6):
        for p2 in (1e6, 3e6, 5e6, 7e6, 7.8e6):
            _assert_consistent(SteadyBlackOilSolver().solve(_pair_case(layout, p1, p2)))


@pytest.mark.parametrize("layout, p1, p2, direct", [
    ("bpr+reducer", 50e5, 20e5, False),      # raised in v1.27.0 (line between them at 10¹⁰ bar)
    ("bpr+reducer", 30e5, 20e5, False),
    ("bpr+reducer", 50e5, 20e5, True),       # no pipe between them
    ("bpr+reducer", 60e5, 30e5, False),
    ("reducer+bpr", 40e5, 30e5, False),
    ("reducer+bpr", 20e5, 20e5, False),      # raised in v1.27.0 (−337 °C in a network at rest)
    ("reducer+bpr", 10e5, 60e5, False),      # contradictory: no flow
    ("reducer+reducer", 40e5, 20e5, False),
    ("reducer+reducer", 10e5, 20e5, False),  # raised in v1.27.0 (pressure 0 Pa)
    ("reducer+reducer", 20e5, 65e5, False),  # raised in v1.27.0
])
def test_compositional_regulator_pairs(layout, p1, p2, direct):
    _assert_consistent(SteadyCompositionalSolver().solve(_pair_case(layout, p1, p2, True, direct)))


def test_back_pressure_then_reducing_regulator_solution():
    """Black-oil, BPR at 5 MPa then reducer at 2 MPa, 8 → 0.5 MPa: more
    reaches the BPR at its set point than the reducer passes at its own, so
    the BPR is wide open and the reducer controls.  The flow equals the
    downstream pipe solved alone from 2 MPa at the reducer's outlet
    temperature (an independent solve; obtained 1.3e-7, tolerance the
    solver's nodal tolerance)."""
    res = SteadyBlackOilSolver().solve(_pair_case("bpr+reducer", 5e6, 2e6))
    _assert_consistent(res)
    bpr, red = res.regulators
    assert (bpr.state, red.state) == ("open", "controlling")
    assert bpr.inlet_pressure_pa > 5e6
    alone = SteadyBlackOilSolver().solve(NetworkCase(
        name="down", fluid_model=BlackOilFluid(**CRUDE), pressure_inlets=(PressureBoundary(5, 2e6),),
        pressure_outlets=(PressureBoundary(9, 0.5e6),),
        thermal_inlets=(ThermalBoundary(5, res.node_temperatures_c[5]),),
        inlet_fluid_bcs=(InletFluidBC(5, **CRUDE),), components=(_pipe("down", 5, 9, 2000.0, U=3.0),)))
    assert red.mass_flow_kg_per_s == pytest.approx(alone.component_flows[0].mass_flow_kg_per_s, rel=NODAL)


def test_contradictory_pair_reports_the_cut_off_line():
    """Reducer at 7.5 MPa then BPR at 7.8 MPa: nothing can pass.  Both
    closed, and the warning says the pressure between them is arbitrary."""
    res = SteadyBlackOilSolver().solve(_pair_case("reducer+bpr", 7.5e6, 7.8e6))
    assert res.converged
    assert [g.state for g in res.regulators] == ["closed", "closed"]
    assert any("cut off by closed regulators" in w and "3, 4" in w for w in res.warnings)
