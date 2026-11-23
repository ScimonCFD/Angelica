"""
Black-oil energy balance: Joule-Thomson cooling of the gas, friction heating
of the liquids, potential energy, isenthalpic chokes (v1.20.0), and the
gas-gravity hydrate check.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary  # noqa: E402
from angelica.core.components import Choke, Pipe  # noqa: E402
from angelica.properties.black_oil import BlackOilFluid, _dh_dp  # noqa: E402
from angelica.solvers import BlackOilSolverSettings, SteadyBlackOilSolver  # noqa: E402

# ── (∂h/∂P)_T against NIST reference equations of state ──────────────────────
# Data READ IN SOURCE: NIST Chemistry WebBook isobaric tables (methane:
# Setzmann & Wagner 1991, Cp uncertainty ~1 %; water: IAPWS-95, 0.1 %).
# Reference value −(∂h/∂P)_T = μ_JT·cp, independent of the gas cp parameter.
# Errors obtained (v1.20.0): methane −2.0 to −2.5 % at 5 MPa, −3.9 to −6.1 %
# at 10 MPa; water −1.9 to −2.8 %.  Tolerance 10 %: the gas term is the
# temperature derivative of a gravity-based Z correlation (Hall-Yarborough with
# Sutton pseudo-criticals), for which a few-percent error is expected; this
# tolerance was set with the results known and is documented as such.
NIST_METHANE = {  # (MPa, K): (mu_JT K/MPa, cp kJ/kg/K)
    (5, 280): (4.5723, 2.6226), (5, 300): (3.9428, 2.5703), (5, 320): (3.4198, 2.5551),
    (10, 280): (3.6727, 3.2424), (10, 300): (3.2606, 3.0023), (10, 320): (2.8769, 2.8780),
}
NIST_WATER_10MPA = {300: (-0.21999, 4.1536), 320: (-0.20827, 4.1580),
                    340: (-0.19760, 4.1672), 360: (-0.18704, 4.1810)}


@pytest.mark.parametrize("key", sorted(NIST_METHANE))
def test_gas_isothermal_throttling_coefficient_vs_nist(key):
    P_mpa, T_k = key
    mu, cp = NIST_METHANE[key]
    gas = BlackOilFluid(api_gravity=35.0, gas_gravity=16.043 / 28.964, gor_sc_m3_per_m3=1e6, wor_sc_m3_per_m3=0.0)
    ref = -mu * cp * 1e-3                                   # m³/kg
    assert _dh_dp(gas._pvt_params(), P_mpa * 1e6, T_k - 273.15) == pytest.approx(ref, rel=0.10)


@pytest.mark.parametrize("T_k", sorted(NIST_WATER_10MPA))
def test_water_isothermal_throttling_coefficient_vs_nist(T_k):
    mu, cp = NIST_WATER_10MPA[T_k]
    water = BlackOilFluid(api_gravity=35.0, gas_gravity=0.7, gor_sc_m3_per_m3=0.0, wor_sc_m3_per_m3=1e6)
    assert _dh_dp(water._pvt_params(), 10e6, T_k - 273.15) == pytest.approx(-mu * cp * 1e-3, rel=0.10)


# ── first law along a pipe (consistency with the model's own enthalpy) ────────

def _case(fluid, dz=0.0, U=0.0, seg=10, p_out=30e5):
    return NetworkCase(
        name="bo-energy", fluid_model=fluid,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=80e5),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=p_out),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=60.0, bc_type="fixed_temperature"),),
        components=(Pipe(start_node=1, end_node=2, diameter_m=0.15, length_m=10_000.0,
                         absolute_roughness_m=4.6e-5, height_change_m=dz, n_thermal_segments=seg,
                         heat_transfer_coefficient_w_per_m2k=U, ambient_temperature_c=15.0),),
    )


DEAD_OIL = BlackOilFluid(api_gravity=30.0, gas_gravity=0.65, gor_sc_m3_per_m3=0.0, wor_sc_m3_per_m3=0.0)
GASSY = BlackOilFluid(api_gravity=40.0, gas_gravity=0.70, gor_sc_m3_per_m3=300.0, wor_sc_m3_per_m3=0.0)


@pytest.mark.parametrize("fluid, dz, seg", [(DEAD_OIL, 0.0, 10), (DEAD_OIL, 200.0, 1), (GASSY, 0.0, 10), (GASSY, 300.0, 1)])
def test_adiabatic_pipe_conserves_h_plus_gz(fluid, dz, seg):
    """Standing: CONSISTENCY with the model's own enthalpy (exact first law).
    Errors obtained: ≤ 0.002 K.  Tolerance 0.01 K = outer-loop tolerance."""
    res = SteadyBlackOilSolver().solve(_case(fluid, dz=dz, seg=seg))
    assert res.converged
    T_out = res.node_temperatures_c[2]
    h_in = fluid.enthalpy_for_link(None, 80e5, 60.0)
    T_exact = fluid.temperature_from_enthalpy_for_link(None, h_in - 9.81 * dz, 30e5, T_out)
    assert T_out == pytest.approx(T_exact, abs=0.01)
    m = res.component_flows[0].mass_flow_kg_per_s
    assert abs(res.global_energy_balance.energy_error_kw) < m * 2000.0 * 0.01 / 1000.0


def test_dead_oil_warms_gassy_oil_cools():
    """Adiabatic liquid: ΔT = ΔP·(1 − T·β)/(ρ·cp), between 0 and the full
    dissipation bound ΔP/(ρ·cp) since 0 < T·β < 1 (physics bound, not a fit;
    obtained 2.03 K against a 2.93 K bound, i.e. T·β ≈ 0.31).  Free gas cools."""
    oil = SteadyBlackOilSolver().solve(_case(DEAD_OIL)).node_temperatures_c[2]
    gas = SteadyBlackOilSolver().solve(_case(GASSY)).node_temperatures_c[2]
    rho = DEAD_OIL.pvt(55e5, 61.0).mixture_density_kg_per_m3
    cp = DEAD_OIL.pvt(55e5, 61.0).mixture_specific_heat_j_per_kg_k
    bound = 50e5 / (rho * cp)
    assert 0.5 * bound < oil - 60.0 < bound
    assert gas < 60.0 - 1.0


def test_switch_off_keeps_temperature():
    s = SteadyBlackOilSolver(black_oil_settings=BlackOilSolverSettings(pipe_enthalpy_terms=False))
    res = s.solve(_case(DEAD_OIL))
    assert res.node_temperatures_c[2] == pytest.approx(60.0, abs=1e-6)


def test_choke_is_isenthalpic_in_black_oil():
    fluid = GASSY
    case = NetworkCase(
        name="bo-choke", fluid_model=fluid,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=80e5),),
        pressure_outlets=(PressureBoundary(node_id=3, pressure_pa=20e5),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=60.0, bc_type="fixed_temperature"),),
        components=(Pipe(start_node=1, end_node=2, diameter_m=0.1, length_m=100.0, absolute_roughness_m=4.6e-5),
                    Choke(start_node=2, end_node=3, diameter_m=0.1, orifice_diameter_m=0.02)),
    )
    res = SteadyBlackOilSolver().solve(case)
    assert res.converged
    P2, P3 = res.node_pressures_pa[2], res.node_pressures_pa[3]
    T2, T3 = res.node_temperatures_c[2], res.node_temperatures_c[3]
    assert P2 - P3 > 20e5 and T3 < T2 - 1.0
    assert fluid.enthalpy_for_link(None, P3, T3) == pytest.approx(fluid.enthalpy_for_link(None, P2, T2), abs=50.0)


def test_per_pipe_fluids_keep_the_case_correlation():
    """Fixed in v1.20.0: per-inlet fluids fell back to the Standing
    correlation and the default gas cp."""
    from angelica.solvers.steady_black_oil import SteadyBlackOilSolver as S
    base = BlackOilFluid(api_gravity=35.0, gas_gravity=0.7, gor_sc_m3_per_m3=100.0, wor_sc_m3_per_m3=0.1,
                         correlation="VASQUEZ_BEGGS", gas_specific_heat_j_per_kg_k=2500.0)
    case = NetworkCase(
        name="per-pipe", fluid_model=base,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=60e5),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=40e5),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=50.0, bc_type="fixed_temperature"),),
        inlet_fluid_bcs=(InletFluidBC(node_id=1, api_gravity=30.0, gas_gravity=0.68,
                                      gor_sc_m3_per_m3=80.0, wor_sc_m3_per_m3=0.2),),
        components=(Pipe(start_node=1, end_node=2, diameter_m=0.1, length_m=1000.0, absolute_roughness_m=4.6e-5),),
    )
    from angelica.core.network import build_network_state
    from angelica.properties.black_oil import BlackOilComposition
    ns = build_network_state(case)
    comp = BlackOilComposition(api_gravity=30.0, gas_gravity=0.68, gor_sc_m3_per_m3=80.0, wor_sc_m3_per_m3=0.2)
    eff = S()._build_mixed_fluid(ns, {0: comp}, comp, 50.0, base_fluid=base)
    f = eff.fluid_for_link(ns.components[0])
    assert f.correlation == "VASQUEZ_BEGGS" and f.gas_specific_heat_j_per_kg_k == 2500.0


# ── gas-gravity hydrate chart ─────────────────────────────────────────────────

def test_katz_chart_reproduces_its_digitised_values():
    """Standing: the chart values are an engineering chart (Katz 1945, GPSA),
    READ IN A SECONDARY SOURCE (Janna 2015, UTP dissertation, Tables 4.1a–4.5a).
    The interpolation must return them exactly at the nodes."""
    from angelica.properties.hydrate import katz_hydrate_pressure, katz_hydrate_temperature
    psi = 6894.757
    for sg, t_f, p_psia in ((0.6, 50.0, 510), (0.7, 70.0, 1937), (0.8, 40.0, 143), (1.0, 75.0, 2400)):
        assert katz_hydrate_pressure(sg, (t_f - 32) / 1.8) / psi == pytest.approx(p_psia, rel=1e-9)
        assert katz_hydrate_temperature(sg, p_psia * psi).temperature_c * 1.8 + 32 == pytest.approx(t_f, abs=1e-6)
    # monotonic in gravity: heavier gas forms hydrate at lower pressure
    assert katz_hydrate_pressure(0.65, 10.0) > katz_hydrate_pressure(0.75, 10.0) > katz_hydrate_pressure(0.85, 10.0)


def test_black_oil_results_carry_hydrate_check():
    fluid = BlackOilFluid(api_gravity=35.0, gas_gravity=0.70, gor_sc_m3_per_m3=150.0, wor_sc_m3_per_m3=0.3)
    case = NetworkCase(
        name="cold-line", fluid_model=fluid,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=80e5),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=60e5),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=30.0, bc_type="fixed_temperature"),),
        components=(Pipe(start_node=1, end_node=2, diameter_m=0.2, length_m=20_000.0, absolute_roughness_m=4.6e-5,
                         heat_transfer_coefficient_w_per_m2k=20.0, ambient_temperature_c=4.0),),
    )
    from angelica.properties.hydrate import katz_hydrate_temperature
    res = SteadyBlackOilSolver().solve(case)
    cf = res.component_flows[0]
    assert cf.hydrate_temperature_c is not None and cf.free_water_fraction > 0.0
    P2, T2 = res.node_pressures_pa[2], res.node_temperatures_c[2]
    assert cf.hydrate_subcooling_k >= katz_hydrate_temperature(0.70, P2).temperature_c - T2 - 1e-9
    assert cf.hydrate_subcooling_k > 0.0                   # cold seabed line: inside the hydrate region
