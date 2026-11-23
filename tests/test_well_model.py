"""
Tests for the ``Well`` component (IPR + production tubing), v1.13.0.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures.well_inflow import oil_rate_sc_m3_per_day  # noqa: E402
from angelica.closures.well_model import (  # noqa: E402
    bottomhole_pressure,
    tubing_wellhead_pressure,
    well_curve,
)
from angelica.core.case import InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary  # noqa: E402
from angelica.core.components import Pipe, Well  # noqa: E402
from angelica.properties.black_oil import BlackOilFluid  # noqa: E402
from angelica.properties.single_component import SingleComponentFluid  # noqa: E402
from angelica.solvers import SteadyBlackOilSolver, SteadyIsothermalIncompressibleSolver  # noqa: E402

BAR = 1e5
RICH = {"api_gravity": 40.0, "gas_gravity": 0.70, "gor_sc_m3_per_m3": 150.0, "wor_sc_m3_per_m3": 0.1}
STAGNANT_KG_PER_S = 1e-5   # flow left in a stagnant pipe by the solver's pressure tolerance


def _well(J=6.0, D=0.0762, L=2100.0, dz=2000.0, model="vogel", **kw) -> Well:
    return Well(start_node=1, end_node=2, diameter_m=D, productivity_index_m3_per_day_per_bar=J,
                ipr_model=model, tubing_length_m=L, tubing_height_change_m=dz, **kw)


@pytest.mark.parametrize("model", ["linear", "vogel"])
def test_bottomhole_pressure_inverts_ipr(model):
    p_res, J, p_b = 250 * BAR, 6.0, 200 * BAR
    for p_wf in (240 * BAR, 180 * BAR, 90 * BAR, 20 * BAR):
        q = oil_rate_sc_m3_per_day(p_res, p_wf, J, model, p_b)
        assert bottomhole_pressure(q, p_res, J, model, p_b) == pytest.approx(p_wf, rel=1e-9)


def test_single_phase_tubing_matches_hydrostatic_plus_darcy():
    """Dead oil (no gas, no water): Beggs-Brill reduces to single-phase flow.
    Compare with an independent integration using the Haaland friction factor."""
    fluid = BlackOilFluid(api_gravity=30.0, gas_gravity=0.65, gor_sc_m3_per_m3=0.0, wor_sc_m3_per_m3=0.0)
    D, L, dz, eps, T, m_dot = 0.0762, 2100.0, 2000.0, 4.6e-5, 70.0, 8.0
    p_wf = 250 * BAR
    p_wh = tubing_wellhead_pressure(p_wf, m_dot, fluid, T, D, L, dz, eps, 400)

    area = math.pi * D * D / 4.0
    sin_t = dz / L
    P, n = p_wf, 4000
    ds = L / n
    for _ in range(n):
        pvt = fluid.pvt(P, T)
        rho, mu = pvt.density_oil_kg_per_m3, pvt.viscosity_oil_pa_s
        v = m_dot / (rho * area)
        Re = rho * v * D / mu
        f = (-1.8 * math.log10((eps / D / 3.7) ** 1.11 + 6.9 / Re)) ** -2   # Haaland
        P -= ds * (rho * 9.81 * sin_t + f * rho * v * v / (2 * D))
    assert p_wh == pytest.approx(P, rel=5e-3)


def test_curve_shape_and_stable_branch():
    fluid = BlackOilFluid(**RICH)
    curve = well_curve(_well(), fluid, 250 * BAR, 85.0)
    p = curve.wellhead_pressures_pa
    # Rich, gassy well: the J-curve rises above the static wellhead pressure first.
    assert curve.peak_index > 0 and curve.max_wellhead_pressure_pa > p[0]
    # Stable branch: wellhead pressure falls monotonically after the peak ...
    tail = p[curve.peak_index:]
    assert all(b < a for a, b in zip(tail, tail[1:]))
    # ... so the link law ṁ(P_wh) decreases with P_wh and has dṁ/dP_wh < 0.
    pts = [curve.mass_flow(x) for x in (2e6, 4e6, 6e6, 9e6)]
    assert all(b[0] < a[0] for a, b in zip(pts, pts[1:]))
    assert all(d < 0.0 for _, d in pts)
    # Above the maximum the well is dead.
    assert curve.mass_flow(curve.max_wellhead_pressure_pa + 1e5) == (0.0, 0.0)


def test_heavy_wet_well_is_dead():
    fluid = BlackOilFluid(api_gravity=25.0, gas_gravity=0.60, gor_sc_m3_per_m3=15.0, wor_sc_m3_per_m3=0.8)
    curve = well_curve(_well(J=3.0), fluid, 170 * BAR, 75.0)
    assert curve.mass_flow(1.5e6)[0] == 0.0


def _network(p_res=250 * BAR, p_sep=15 * BAR) -> NetworkCase:
    return NetworkCase(
        name="one-well",
        fluid_model=BlackOilFluid(**RICH),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=p_res),),
        pressure_outlets=(PressureBoundary(node_id=3, pressure_pa=p_sep),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=85.0, bc_type="fixed_temperature"),),
        inlet_fluid_bcs=(InletFluidBC(node_id=1, **RICH),),
        components=(_well(),
                    Pipe(start_node=2, end_node=3, diameter_m=0.1016, length_m=3000.0,
                         absolute_roughness_m=4.6e-5)),
    )


def test_network_operating_point_lies_on_the_curve():
    res = SteadyBlackOilSolver().solve(_network())
    assert res.converged
    m = res.component_flows[0].mass_flow_kg_per_s
    curve = well_curve(_well(), BlackOilFluid(**RICH), 250 * BAR, 85.0)
    assert m == pytest.approx(curve.mass_flow(res.node_pressures_pa[2])[0], rel=1e-4)
    assert m > 1.0
    assert res.global_energy_balance.energy_error_pct < 1e-3


def test_dead_well_in_network_produces_nothing():
    res = SteadyBlackOilSolver().solve(_network(p_res=120 * BAR))   # static column > P_res − P_sep
    assert res.converged
    assert res.component_flows[0].mass_flow_kg_per_s == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize("correlation", ["colebrook", "beggs_brill"])
@pytest.mark.parametrize("p_sep_over_max", [1.01, 1.4, 1.3 * 220 / 85.5])   # last: above P_res
def test_well_dead_by_back_pressure_converges(correlation, p_sep_over_max):
    """v1.23.1.  The separator is above the maximum of the deliverability
    curve, so the well cannot flow.  Before the fix the flowline, stagnant at
    ΔP = 0, carried a phantom ρ·A·1e-6 kg/s (the iteration floor stored as
    its velocity); with dense oil at the separator that exceeded the
    'network at rest' threshold and the solve never converged."""
    from angelica.closures import BeggsBrillCorrelation

    p_res = 220 * BAR
    p_max = well_curve(_well(), BlackOilFluid(**RICH), p_res, 85.0).max_wellhead_pressure_pa
    assert p_max == pytest.approx(85.5 * BAR, rel=0.01)
    p_sep = p_sep_over_max * p_max
    solver = (SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation())
              if correlation == "beggs_brill" else SteadyBlackOilSolver())
    case = _network(p_res=p_res, p_sep=p_sep)
    res = solver.solve(case)
    assert res.converged
    well, flowline = res.component_flows
    assert well.mass_flow_kg_per_s == 0.0
    # The flowline carries the laminar flow of the residual ΔP the solver
    # accepts (pressure tolerance 1e-3 Pa → ~1e-6 kg/s in 3 km of 4-in. pipe).
    assert flowline.mass_flow_kg_per_s == pytest.approx(0.0, abs=STAGNANT_KG_PER_S)
    assert res.node_pressures_pa[2] == pytest.approx(p_sep, abs=1.0)


def _two_well_network(p_weak: float, p_sep: float) -> NetworkCase:
    """A strong well (300 bar) and a weak one sharing a flowline to the separator."""
    def flowline(a, b, L):
        return Pipe(start_node=a, end_node=b, diameter_m=0.1016, length_m=L, absolute_roughness_m=4.6e-5)
    weak = Well(start_node=4, end_node=5, diameter_m=0.0762, productivity_index_m3_per_day_per_bar=6.0,
                tubing_length_m=2100.0, tubing_height_change_m=2000.0)
    return NetworkCase(
        name="two-wells",
        fluid_model=BlackOilFluid(**RICH),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=300 * BAR),
                         PressureBoundary(node_id=4, pressure_pa=p_weak)),
        pressure_outlets=(PressureBoundary(node_id=9, pressure_pa=p_sep),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=85.0, bc_type="fixed_temperature"),
                        ThermalBoundary(node_id=4, temperature_c=85.0, bc_type="fixed_temperature")),
        inlet_fluid_bcs=(InletFluidBC(node_id=1, **RICH), InletFluidBC(node_id=4, **RICH)),
        components=(_well(), weak, flowline(2, 7, 1000.0), flowline(5, 7, 1000.0), flowline(7, 9, 8000.0)),
    )


def test_weak_well_dies_against_the_manifold_while_the_other_produces():
    """The strong well holds the manifold at ~5.9 MPa, above the weak well's
    maximum wellhead pressure (5.5 MPa): the weak well is dead, its flowline
    stagnant, and all the oil comes from the strong well."""
    res = SteadyBlackOilSolver().solve(_two_well_network(p_weak=180 * BAR, p_sep=30 * BAR))
    assert res.converged
    m_strong, m_weak, _l1, l2, trunk = (cf.mass_flow_kg_per_s for cf in res.component_flows)
    p_max_weak = well_curve(_well(), BlackOilFluid(**RICH), 180 * BAR, 85.0).max_wellhead_pressure_pa
    assert res.node_pressures_pa[5] > p_max_weak
    assert m_weak == 0.0 and l2 == pytest.approx(0.0, abs=STAGNANT_KG_PER_S)
    assert m_strong > 5.0
    assert trunk == pytest.approx(m_strong, rel=1e-3)     # nodal tolerance of the solver
    assert all(math.isfinite(t) and 0.0 < t < 100.0 for t in res.node_temperatures_c.values())


def _states(res) -> dict:
    return {ws.component_id: ws for ws in res.well_states}


def test_bistable_well_left_dead_when_it_cannot_sustain_its_minimum_rate():
    """v1.23.2.  Between its static wellhead pressure (3.04 MPa) and the
    maximum of its curve (5.49 MPa at 0.94 kg/s) the weak well is bistable.
    Producing is inconsistent: at its minimum stable rate the manifold rises
    above 5.49 MPa.  Dead is consistent: with it shut in the manifold is at
    5.42 MPa, above its static pressure, so it cannot start.  v1.23.1
    reported converged=False here (and before that, an OverflowError)."""
    res = SteadyBlackOilSolver().solve(_two_well_network(p_weak=180 * BAR, p_sep=15 * BAR))
    assert res.converged
    m_strong, m_weak, _l1, l2, trunk = (cf.mass_flow_kg_per_s for cf in res.component_flows)
    weak = res.well_states[1]
    assert weak.state == "dead_minimum_rate"
    assert weak.static_wellhead_pressure_pa < weak.wellhead_pressure_pa < weak.max_wellhead_pressure_pa
    assert weak.minimum_stable_rate_kg_per_s == pytest.approx(0.94, abs=0.01)
    assert "unloading or gas lift" in weak.message
    assert m_weak == 0.0 and l2 == pytest.approx(0.0, abs=STAGNANT_KG_PER_S)
    assert m_strong > 5.0 and trunk == pytest.approx(m_strong, rel=1e-3)
    assert res.well_states[0].state == "producing"
    assert all(math.isfinite(t) and 0.0 < t < 100.0 for t in res.node_temperatures_c.values())


def test_bistable_well_produces_when_producing_is_consistent():
    """Same network, weak reservoir at 200 bar: its wellhead (6.30 MPa) lies
    in the bistable band (5.31–6.97 MPa) and producing is consistent, so the
    well produces — the result v1.23.1 already gave (7.876 / 1.946 kg/s)."""
    res = SteadyBlackOilSolver().solve(_two_well_network(p_weak=200 * BAR, p_sep=15 * BAR))
    assert res.converged
    weak = res.well_states[1]
    assert weak.state == "producing"
    assert weak.static_wellhead_pressure_pa < weak.wellhead_pressure_pa < weak.max_wellhead_pressure_pa
    assert res.component_flows[0].mass_flow_kg_per_s == pytest.approx(7.876, rel=1e-3)
    assert res.component_flows[1].mass_flow_kg_per_s == pytest.approx(1.946, rel=1e-3)


def test_well_with_no_steady_state_is_reported_unstable():
    """The case without a steady state.  Weak reservoir at 205 bar on its own
    1-in., 3 km flowline, separator at 5 bar.  Dead, its wellhead (5.25 MPa)
    is below its static pressure (5.87 MPa), so it would start; producing,
    its minimum stable rate (0.96 kg/s) through the thin line lifts the
    wellhead above its maximum (7.34 MPa), so it would die.  The solver
    says so instead of reporting either state as converged."""
    import dataclasses

    case = _two_well_network(p_weak=205 * BAR, p_sep=5 * BAR)
    comps = list(case.components)
    comps[3] = Pipe(start_node=5, end_node=7, diameter_m=0.0254, length_m=3000.0, absolute_roughness_m=4.6e-5)
    res = SteadyBlackOilSolver().solve(dataclasses.replace(case, components=tuple(comps)))
    assert not res.converged
    weak = res.well_states[1]
    assert weak.state == "unstable"
    assert weak.wellhead_pressure_pa < weak.static_wellhead_pressure_pa
    assert "no steady state" in weak.message
    assert all(math.isfinite(t) for t in res.node_temperatures_c.values())


def test_back_pressure_dead_state_is_reported():
    res = SteadyBlackOilSolver().solve(_two_well_network(p_weak=180 * BAR, p_sep=30 * BAR))
    assert res.converged
    assert res.well_states[1].state == "dead_back_pressure"
    assert res.well_states[1].message.startswith("dead: needs P_wh <")


def test_well_rejected_outside_black_oil():
    case = NetworkCase(
        name="well-water",
        fluid_model=SingleComponentFluid(density_kg_per_m3=1000.0, viscosity_pa_s=1e-3),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=3e5),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=1e5),),
        components=(_well(),),
    )
    with pytest.raises(ValueError, match="black-oil"):
        SteadyIsothermalIncompressibleSolver().solve(case)


def test_component_validation():
    with pytest.raises(ValueError):
        Well(start_node=1, end_node=2, diameter_m=0.07, productivity_index_m3_per_day_per_bar=5.0)
    with pytest.raises(ValueError):
        _well(dz=2500.0)                      # rise longer than the tubing
    assert _well(dz=float("nan")).tubing_height_change_m == 2100.0   # default: vertical


def test_gas_well_with_backpressure_ipr():
    """Back-pressure IPR in a Well: the deliverability curve's bottom-hole
    pressure inverts the IPR at the gas rate GOR·q_o (consistency)."""
    from angelica.closures.well_inflow import IPR, ipr_rate, mass_per_rate_unit
    from angelica.closures.well_model import well_tubing_state

    gas = {"api_gravity": 55.0, "gas_gravity": 0.65, "gor_sc_m3_per_m3": 3000.0, "wor_sc_m3_per_m3": 0.0}
    fluid = BlackOilFluid(**gas)
    well = Well(start_node=1, end_node=2, diameter_m=0.0762, ipr_model="backpressure", ipr_coefficient=40.0,
                ipr_exponent=0.8, tubing_length_m=2500.0, tubing_correlation="gray")
    curve = well_curve(well, fluid, 200 * BAR, 90.0)
    # at very low rates condensate loads the tubing (liquid loading), so even
    # this gas well has a small hump; the network uses the branch beyond it
    assert curve.peak_index < len(curve.oil_rates_sc_m3_per_day) // 3
    m, _ = curve.mass_flow(60 * BAR)
    assert m > 0.0
    state = well_tubing_state(well, fluid, 200 * BAR, 90.0, m)
    assert state.wellhead_pressure_pa == pytest.approx(60 * BAR, rel=0.01)
    # bottom-hole pressure from the IPR at this rate
    ipr = IPR("backpressure", coefficient=40.0, exponent=0.8)
    q_gas = m * 86_400.0 / mass_per_rate_unit(fluid, ipr)
    assert q_gas == pytest.approx(3000.0 * m * 86_400.0 / curve.mass_per_oil_sc_m3, rel=1e-9)
    assert ipr_rate(ipr, 200 * BAR, 150 * BAR) > 0.0
