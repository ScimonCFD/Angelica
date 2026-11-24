"""
Gas lift in ``Well`` (roadmap item 4.4) and standard-volume black-oil mixing.

Standing of the checks
----------------------
- Black-oil mixing: standard volumes are conserved exactly (by construction,
  checked); the mixed composition's enthalpy equals the mass average of the
  streams' within 0.01 K·cp (obtained), against 2–4 K·cp with the old
  mass-averaged parameters — the reason tutorial temperatures changed.
- Gas lift, tubing: continuity as the lift rate vanishes (1 Sm³/d changes
  P_wh by 0.003 %, obtained); first law of the two-section march
  (residual 0.01 W of 500 kW, obtained); lift mass rate = q·γ·ρ_air.
- Network: the flowline carries production + lift gas; mass and energy
  balances close; compositional flowline composition is the molar mix.
- NOT VALIDATED against published gas-lift performance data (Brown's
  gradient curves, Economides' examples not accessible).  The lift gas is
  assumed to reach the valve at the formation temperature and not to
  exchange mass with the produced fluid.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures import well_model as wm  # noqa: E402
from angelica.closures.well_fluid import BlackOilWellFluid, lift_gas_fluid  # noqa: E402
from angelica.core.case import (  # noqa: E402
    InletCompositionBC,
    InletFluidBC,
    NetworkCase,
    PressureBoundary,
    ThermalBoundary,
)
from angelica.core.components import Pipe, Well  # noqa: E402
from angelica.properties.black_oil import (  # noqa: E402
    BlackOilComposition,
    BlackOilFluid,
    GasStream,
    mix_black_oil_streams,
)
from angelica.properties.compositional_fluid import CompositionalFluid  # noqa: E402
from angelica.solvers import SteadyBlackOilSolver, SteadyCompositionalSolver  # noqa: E402

G = 9.80665
# The solver's nodal tolerance is 1e-3 of (inflow + outflow), i.e. 2e-3 of the
# flow through a node.
NODAL = 2e-3
HEAVY = {"api_gravity": 25.0, "gas_gravity": 0.60, "gor_sc_m3_per_m3": 15.0, "wor_sc_m3_per_m3": 0.8}
RICH = {"api_gravity": 40.0, "gas_gravity": 0.70, "gor_sc_m3_per_m3": 150.0, "wor_sc_m3_per_m3": 0.1}


def _well(q_gl=0.0, depth=1800.0, J=3.0, **kw) -> Well:
    lift = {"gas_lift_rate_sm3_per_day": q_gl, "gas_lift_depth_m": depth} if q_gl else {}
    return Well(component_id="w", start_node=1, end_node=2, diameter_m=0.0762,
                productivity_index_m3_per_day_per_bar=J, tubing_length_m=2100.0, tubing_height_change_m=2000.0,
                **lift, **kw)


# ── black-oil stream mixing ──────────────────────────────────────────────────

def _oil_volume(m, c):
    from angelica.properties.dead_oil import dead_oil_density_kg_per_m3
    return m / (dead_oil_density_kg_per_m3(c.api_gravity) + c.gor_sc_m3_per_m3 * c.gas_gravity * 1.2250
                + c.wor_sc_m3_per_m3 * 1025.0)


def test_black_oil_mixing_conserves_standard_volumes_and_enthalpy():
    a, b = BlackOilComposition(32.0, 0.65, 25.0, 0.3), BlackOilComposition(22.0, 0.70, 10.0, 1.5)
    ma, mb = 102.4, 65.1
    mix = mix_black_oil_streams([(ma, a), (mb, b)], a)
    qo = _oil_volume(ma, a) + _oil_volume(mb, b)
    assert _oil_volume(ma + mb, mix) == pytest.approx(qo, rel=1e-9)
    assert mix.gor_sc_m3_per_m3 * qo == pytest.approx(25.0 * _oil_volume(ma, a) + 10.0 * _oil_volume(mb, b), rel=1e-9)

    def fluid(c):
        return BlackOilFluid(api_gravity=c.api_gravity, gas_gravity=c.gas_gravity,
                             gor_sc_m3_per_m3=c.gor_sc_m3_per_m3, wor_sc_m3_per_m3=c.wor_sc_m3_per_m3)

    for P, T in ((4.76e6, 63.0), (8e6, 80.0)):
        h_avg = (ma * fluid(a).enthalpy_for_link(None, P, T) + mb * fluid(b).enthalpy_for_link(None, P, T)) / (ma + mb)
        cp = fluid(mix).specific_heat_at(None, P, T)
        assert abs(fluid(mix).enthalpy_for_link(None, P, T) - h_avg) / cp < 0.05     # K


def test_gas_stream_adds_gas_only():
    oil = BlackOilComposition(25.0, 0.60, 15.0, 0.8)
    m_oil, q_gl = 2.0, 50_000.0
    m_gl = q_gl * 0.60 * 1.2250 / 86_400.0
    mix = mix_black_oil_streams([(m_oil, oil), (m_gl, GasStream(0.60))], oil)
    q_o = _oil_volume(m_oil, oil) * 86_400.0
    assert mix.gor_sc_m3_per_m3 == pytest.approx(15.0 + q_gl / q_o, rel=1e-9)
    assert (mix.api_gravity, mix.wor_sc_m3_per_m3) == pytest.approx((25.0, 0.8))
    gas_only = mix_black_oil_streams([(m_gl, GasStream(0.65))], oil)
    assert gas_only.gor_sc_m3_per_m3 >= 1e6 and gas_only.gas_gravity == pytest.approx(0.65)


# ── tubing with gas lift ─────────────────────────────────────────────────────

def test_lift_rate_to_zero_recovers_the_natural_well():
    fluid = BlackOilFluid(**RICH)
    natural = wm.well_tubing_state(_well(), fluid, 250e5, 85.0, 3.0)
    tiny = wm.well_tubing_state(_well(1.0), fluid, 250e5, 85.0, 3.0)
    assert tiny.wellhead_pressure_pa == pytest.approx(natural.wellhead_pressure_pa, rel=1e-4)
    assert tiny.wellhead_temperature_c == pytest.approx(natural.wellhead_temperature_c, abs=0.01)


def test_two_section_march_satisfies_the_first_law():
    fluid = BlackOilFluid(**HEAVY)
    well = _well(50_000.0)
    st = wm.well_tubing_state(well, fluid, 17e6, 75.0, 2.0)
    wf = BlackOilWellFluid(fluid)
    gas = lift_gas_fluid(wf, well)
    m_gl = st.lift_mass_flow_kg_per_s
    assert m_gl == pytest.approx(50_000.0 * 0.60 * 1.2250 / 86_400.0, rel=1e-6)
    h_in = 2.0 * wf.enthalpy(17e6, 75.0) + m_gl * gas.enthalpy(st.injection_pressure_pa, st.injection_temperature_c)
    h_out = ((2.0 + m_gl) * st.wellhead_enthalpy_j_per_kg + st.formation_heat_w
             + G * (2.0 * 2000.0 + m_gl * 2000.0 * 1800.0 / 2100.0))
    assert h_in == pytest.approx(h_out, abs=1.0)                                # W
    assert 0.0 < st.injection_pressure_pa < 17e6


def test_lift_revives_a_dead_well_and_has_an_optimum():
    """The heavy, wet well is dead without lift; oil rate at 1.5 MPa wellhead
    pressure against lift rate rises, peaks and falls (friction of excess
    gas) — the shape of a gas-lift performance curve."""
    fluid = BlackOilFluid(**HEAVY)
    assert wm.well_curve(_well(), fluid, 17e6, 75.0).max_wellhead_pressure_pa <= 0.0
    rates = []
    for q in (5_000.0, 20_000.0, 50_000.0, 100_000.0, 200_000.0):
        c = wm.well_curve(_well(q), fluid, 17e6, 75.0)
        rates.append(c.mass_flow(1.5e6)[0])
    assert rates[0] < rates[1] < rates[2] < rates[3] and rates[4] < rates[3]


def test_lift_validation():
    with pytest.raises(ValueError, match="gas_lift_depth_m"):
        _well(10_000.0, depth=float("nan"))
    with pytest.raises(ValueError, match="gas_lift_depth_m"):
        _well(10_000.0, depth=2500.0)
    with pytest.raises(ValueError, match="sum to 1"):
        _well(10_000.0, gas_lift_zs=(0.5, 0.4))


# ── networks ─────────────────────────────────────────────────────────────────

def _bo_network(q_gl, p_sep=1.0e6) -> NetworkCase:
    return NetworkCase(
        name="gas lift", fluid_model=BlackOilFluid(**HEAVY),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=17e6),),
        pressure_outlets=(PressureBoundary(node_id=3, pressure_pa=p_sep),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=75.0, bc_type="fixed_temperature"),),
        inlet_fluid_bcs=(InletFluidBC(node_id=1, **HEAVY),),
        components=(_well(q_gl), Pipe(component_id="fl", start_node=2, end_node=3, diameter_m=0.1016,
                                      length_m=2000.0, absolute_roughness_m=4.6e-5,
                                      heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0)))


def test_black_oil_network_with_gas_lift():
    res = SteadyBlackOilSolver().solve(_bo_network(20_000.0))
    assert res.converged
    m_well, m_line = (cf.mass_flow_kg_per_s for cf in res.component_flows)
    m_gl = 20_000.0 * 0.60 * 1.2250 / 86_400.0
    assert res.well_states[0].state == "producing"
    assert res.well_states[0].gas_lift_kg_per_s == pytest.approx(m_gl, rel=1e-6)
    assert m_line == pytest.approx(m_well + m_gl, rel=NODAL)
    assert res.global_balance.mass_inlet_kg_per_s == pytest.approx(m_line, rel=NODAL)
    geb = res.global_energy_balance
    assert abs(geb.energy_error_kw) < m_line * 2500.0 * 0.01 / 1e3              # ṁ·cp·0.01 K


def test_dead_well_still_passes_its_lift_gas():
    res = SteadyBlackOilSolver().solve(_bo_network(20_000.0, p_sep=12.5e6))
    assert res.converged
    m_well, m_line = (cf.mass_flow_kg_per_s for cf in res.component_flows)
    assert m_well == 0.0 and res.well_states[0].state == "dead_back_pressure"
    assert m_line == pytest.approx(20_000.0 * 0.60 * 1.2250 / 86_400.0, rel=NODAL)


def test_compositional_network_with_gas_lift():
    """Volatile oil lifted with the stream's own standard-conditions gas: the
    flowline composition is the molar mix of production and lift gas."""
    names = ("methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane", "n-decane", "n-hexadecane")
    zs = (0.30, 0.06, 0.05, 0.04, 0.03, 0.04, 0.28, 0.20)
    lift_zs = (0.90, 0.06, 0.03, 0.01, 0.0, 0.0, 0.0, 0.0)
    well = Well(component_id="ow", start_node=1, end_node=2, diameter_m=0.0762, ipr_model="linear",
                productivity_index_m3_per_day_per_bar=3.0, tubing_length_m=2500.0, tubing_height_change_m=2500.0,
                gas_lift_rate_sm3_per_day=30_000.0, gas_lift_depth_m=2000.0, gas_lift_zs=lift_zs)
    fluid = CompositionalFluid(names, zs)
    case = NetworkCase(
        name="gas lift EOS", fluid_model=fluid,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=230e5),),
        pressure_outlets=(PressureBoundary(node_id=3, pressure_pa=10e5),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=90.0, bc_type="fixed_temperature"),),
        inlet_composition_bcs=(InletCompositionBC(node_id=1, zs=zs),),
        components=(well, Pipe(component_id="fl", start_node=2, end_node=3, diameter_m=0.1016, length_m=2000.0,
                               absolute_roughness_m=4.6e-5, heat_transfer_coefficient_w_per_m2k=5.0,
                               ambient_temperature_c=15.0)))
    res = SteadyCompositionalSolver().solve(case)
    assert res.converged
    m_well, m_line = (cf.mass_flow_kg_per_s for cf in res.component_flows)
    m_gl = res.well_states[0].gas_lift_kg_per_s
    assert m_line == pytest.approx(m_well + m_gl, rel=NODAL)
    mws = fluid.component_mws
    n_res = m_well / sum(z * m for z, m in zip(zs, mws))
    n_gl = m_gl / sum(z * m for z, m in zip(lift_zs, mws))
    expected_c1 = (n_res * zs[0] + n_gl * lift_zs[0]) / (n_res + n_gl)
    assert res.component_flows[1].zs[0] == pytest.approx(expected_c1, rel=NODAL)
    assert abs(res.global_energy_balance.energy_error_kw) < m_line * 2500.0 * 0.01 / 1e3
