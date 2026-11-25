"""
Black-oil correlations at the edge of and outside their data range (v1.27.1).

Standing of the checks
----------------------
- Data range of the Beggs-Robinson dead-oil viscosity: 16–58 °API,
  70–295 °F (21.1–146.1 °C) — read in secondary sources (pengtools wiki;
  IHS Harmony reference); the original paper (SPE-5434-PA) was not
  accessible.  Not a validation of the correlation: the tests check that the
  solver says when it is outside that range.
- Enthalpy below the numerical floor (−6.7 °C): continuous and increasing
  in T (consistency).  Until v1.27.1 h(8 MPa) fell by 35 kJ/kg between
  −7.2 and −6.2 °C, so a let-down below the floor heated the fluid by
  16.6 K and the inversion of h found a node at −22.7 °C.
- Laminar-turbulent transition: the friction factor jumps from 64/Re =
  0.0278 to Colebrook 0.0475 at Re = 2300 (ε/D = 3e-4); a viscous line
  near it can have no solution.  Checked: such a case ends converged=False
  and names the pipe and the two regimes.  Not fixed (a modelling choice for
  the maintainer; see the ROADMAP).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary  # noqa: E402
from angelica.core.components import Pipe, PressureRegulator  # noqa: E402
from angelica.properties.black_oil import (  # noqa: E402
    T_MIN_C,
    BlackOilFluid,
    black_oil_enthalpy,
    temperature_from_enthalpy,
)
from angelica.properties.dead_oil import BEGGS_ROBINSON_T_RANGE_C  # noqa: E402
from angelica.solvers import SteadyBlackOilSolver  # noqa: E402

CRUDE = {"api_gravity": 35.0, "gas_gravity": 0.70, "gor_sc_m3_per_m3": 100.0, "wor_sc_m3_per_m3": 0.2}
PARAMS = (35.0, 0.7, 100.0, 0.2, "STANDING", 2200.0, 0.035)


def _pipe(cid, a, b, L, T_amb):
    return Pipe(component_id=cid, start_node=a, end_node=b, diameter_m=0.1524, length_m=L,
                absolute_roughness_m=46e-6, heat_transfer_coefficient_w_per_m2k=3.0, ambient_temperature_c=T_amb)


def _line(T_amb, regulator=True, p_in=8e6, T_in=60.0, api=35.0):
    """20 km → reducer (2 MPa) → 20 km, 8 → 0.5 MPa (the review's case)."""
    crude = dict(CRUDE, api_gravity=api)
    comps = [_pipe("p1", 1, 2, 20000.0, T_amb)]
    if regulator:
        comps += [PressureRegulator(component_id="r1", start_node=2, end_node=3, diameter_m=0.1524,
                                    set_pressure_pa=2e6), _pipe("p3", 3, 9, 20000.0, T_amb)]
    else:
        comps = [_pipe("p3", 1, 9, 20000.0, T_amb)]
    return NetworkCase(name="line", fluid_model=BlackOilFluid(**crude),
                       pressure_inlets=(PressureBoundary(1, p_in),), pressure_outlets=(PressureBoundary(9, 0.5e6),),
                       thermal_inlets=(ThermalBoundary(1, T_in),), inlet_fluid_bcs=(InletFluidBC(1, **crude),),
                       components=tuple(comps))


# ── enthalpy below the numerical floor ────────────────────────────────────────

def test_enthalpy_is_continuous_and_increasing_through_the_floor():
    Ts = (-12.0, -10.0, -8.0, -7.2, -6.7, -6.2, -5.0, -3.0)
    for P in (2e6, 8e6):
        h = [black_oil_enthalpy(PARAMS, P, T) for T in Ts]
        assert all(b > a for a, b in zip(h, h[1:]))
        # no jump at the floor: the step −7.2 → −6.2 °C is ~1 K of cp
        # (obtained 2.27 kJ/kg; it was −35.4 kJ/kg at 8 MPa)
        step = black_oil_enthalpy(PARAMS, P, -6.2) - black_oil_enthalpy(PARAMS, P, -7.2)
        assert 1.5e3 < step < 3.5e3


def test_let_down_below_the_floor_continues_the_one_above_it():
    """Isenthalpic 7.97 → 2 MPa.  Below the floor the properties are frozen,
    so the temperature change equals the one at the floor (obtained +0.12 K
    at −6.2, −8 and −10 °C; it was +16.6 K below the floor and −11 K just
    above it).  Tolerance 0.02 K, ten times the inversion's resolution."""
    def let_down(T_in):
        h = black_oil_enthalpy(PARAMS, 7.97e6, T_in)
        return temperature_from_enthalpy(lambda T: black_oil_enthalpy(PARAMS, 2e6, round(T, 4)), h, T_in) - T_in
    at_floor = let_down(-6.2)
    for T_in in (-8.0, -10.0):
        assert let_down(T_in) == pytest.approx(at_floor, abs=0.02)
    assert abs(at_floor) < 0.5


@pytest.mark.parametrize("T_amb", [-6.0, -8.0])
def test_cold_regulated_line_has_physical_temperatures(T_amb):
    """v1.27.0: node 2 at −22.7 °C (ambient −6) and +16.6 K across the
    regulator (ambient −8).  Now every node lies between the ambient and the
    inlet temperature, the let-down changes T by less than 0.5 K (+0.1 K
    obtained), and the solve converges with range warnings."""
    res = SteadyBlackOilSolver().solve(_line(T_amb))
    assert res.converged
    T = res.node_temperatures_c
    assert all(T_amb - 0.5 <= t <= 60.0 for t in T.values())
    assert abs(T[3] - T[2]) < 0.5
    assert any("extrapolated below" in w for w in res.warnings)
    if T_amb < T_MIN_C:
        assert any("frozen at their numerical limits" in w for w in res.warnings)


# ── range warnings, not convergence failures ──────────────────────────────────

def test_warnings_outside_the_data_range_and_none_inside():
    lo, hi = BEGGS_ROBINSON_T_RANGE_C
    warm = SteadyBlackOilSolver().solve(_line(20.0))
    assert warm.converged and min(warm.node_temperatures_c.values()) >= lo
    assert not [w for w in warm.warnings if "dead-oil" in w]
    cold = SteadyBlackOilSolver().solve(_line(10.0))
    assert cold.converged and cold.diagnosis is None
    assert any(f"({lo} °C)" in w for w in cold.warnings)
    heavy = SteadyBlackOilSolver().solve(_line(20.0, api=12.0))
    assert any("°API" in w and "12" in w for w in heavy.warnings)


# ── laminar-turbulent transition ──────────────────────────────────────────────

def test_laminar_turbulent_cycle_is_diagnosed():
    """The downstream pipe alone, 2 → 0.5 MPa, inlet 13 °C, ambient 5 °C:
    Re cycles 2184 → 2306 → 3219 → 2334 → 2566 (f 0.029 ↔ 0.047) over the
    outer iterations.  Not converged, and the diagnosis says why."""
    res = SteadyBlackOilSolver().solve(_line(5.0, regulator=False, p_in=2e6, T_in=13.0))
    assert not res.converged
    sw = [s for s in res.diagnosis.pattern_switches if s.component_id == "p3"]
    assert sw and set(sw[0].patterns) == {"laminar", "turbulent"}


def test_ambient_sweep_of_the_regulated_line():
    """Ambient 20 → −5 °C in steps of 1 K.  Every case ends converged or
    not converged with the laminar-turbulent cycle named (6, 5 and 2 °C);
    none raises.  The flow falls with the ambient (6.24 kg/s at 20 °C, 0.010
    at −5 °C) except across the transition."""
    flows = {}
    for T_amb in range(20, -6, -1):
        res = SteadyBlackOilSolver().solve(_line(float(T_amb)))
        if res.converged:
            flows[T_amb] = res.regulators[0].mass_flow_kg_per_s
        else:
            regimes = {p for s in res.diagnosis.pattern_switches for p in s.patterns}
            assert regimes == {"laminar", "turbulent"}, res.diagnosis.message
    assert len(flows) >= 20
    assert flows[20] == pytest.approx(6.239, rel=1e-3) and flows[-5] == pytest.approx(0.010, rel=0.1)
    temps = sorted(flows)
    assert all(flows[a] < flows[b] for a, b in zip(temps, temps[1:]))
