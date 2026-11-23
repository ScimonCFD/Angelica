"""
Tubing temperature profile of the ``Well`` component (roadmap item 4.1).

Standing of the checks
----------------------
- Ramey's solution: the differential equation and its closed form are derived
  in ``angelica.closures.well_model`` (linear ODE; own derivation).  The
  relaxation distance and the time function f(t) were READ IN A SECONDARY
  SOURCE (OpenGeoSys benchmark documentation "Wellbore Heat Transport -
  EUBHE", which reproduces Ramey 1962).  The comparison below is therefore a
  verification of the integrator against the analytical solution, with a
  numerical tolerance (0.001 K; error obtained < 1e-9 K, the exponential
  integrator being exact for constant coefficients).
- Black-oil wells: first-law consistency only.  NOT VALIDATED against
  measured flowing temperatures (Sagar, Doty & Schmidt 1991 field cases were
  not accessible).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures.well_model import (  # noqa: E402
    TubingThermal,
    formation_conductance_w_per_mk,
    integrate_tubing,
    ramey_time_function,
    well_curve,
    well_tubing_state,
)
from angelica.core.case import InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary  # noqa: E402
from angelica.core.components import Pipe, Well, WellInflow  # noqa: E402
from angelica.properties.black_oil import BlackOilFluid  # noqa: E402
from angelica.solvers import SteadyBlackOilSolver  # noqa: E402

BAR = 1e5
G = 9.80665
RICH = {"api_gravity": 40.0, "gas_gravity": 0.70, "gor_sc_m3_per_m3": 150.0, "wor_sc_m3_per_m3": 0.1}


def _well(J=6.0, D=0.0762, L=2100.0, dz=2000.0, **kw) -> Well:
    return Well(start_node=1, end_node=2, diameter_m=D, productivity_index_m3_per_day_per_bar=J,
                tubing_length_m=L, tubing_height_change_m=dz, **kw)


# ── Ramey's analytical solution ──────────────────────────────────────────────

@pytest.mark.parametrize("m_dot, n_seg", [(5.0, 25), (0.5, 25), (0.01, 10), (20.0, 50)])
def test_incompressible_static_column_reproduces_ramey(m_dot, n_seg):
    """Incompressible fluid, static column: T(s) = T_e(s) + γ·A·(1 − e^(−s/A)),
    A = ṁ·cp/C, starting at the formation temperature."""
    rho, cp, L, dz = 850.0, 2100.0, 2500.0, 2300.0
    C, gamma, t_res = 3.0, 0.03 * dz / L, 90.0
    thermal = TubingThermal(C, gamma, t_res)
    res = integrate_tubing(
        300 * BAR, t_res, m_dot, 0.0762, L, dz, n_seg,
        gradient=lambda P, T: rho * G * dz / L,
        thermal=thermal, dh_dp=lambda P, T: 1.0 / rho, cp=lambda P, T: cp, keep_profile=True)
    A = m_dot * cp / C
    for s, _P, T, T_e in res.profile:
        assert T == pytest.approx(T_e + gamma * A * (1.0 - math.exp(-s / A)), abs=1e-3)
    # heat given to the formation = sensible heat lost (dh = cp·dT + dP/ρ, and
    # dP/ρ = −g·dz exactly cancels the potential energy)
    assert res.formation_heat_w == pytest.approx(m_dot * cp * (t_res - res.wellhead_temperature_c), rel=1e-9)


def test_time_function_and_conductance():
    # the two published branches do not meet exactly at t_D = 1.5: 0.874 and
    # 0.853, a 2.5 % step that belongs to the published formula
    assert ramey_time_function(1.5) == pytest.approx(0.8740, abs=1e-4)
    assert ramey_time_function(1.5 + 1e-9) == pytest.approx(0.8526, abs=1e-4)
    f = [ramey_time_function(t) for t in (0.01, 0.1, 1.0, 10.0, 1e3, 1e5)]
    assert all(b > a for a, b in zip(f, f[1:]))
    # resistances in series: 1/C = 1/(2π·r·U) + f/(2π·k_e)
    D, Dw, U, k, alpha, t = 0.0762, 0.2159, 40.0, 2.4, 1e-6, 365 * 86400.0
    ft = ramey_time_function(alpha * t / (0.5 * Dw) ** 2)
    C = formation_conductance_w_per_mk(D, Dw, U, k, alpha, t)
    assert 1.0 / C == pytest.approx(1.0 / (math.pi * D * U) + ft / (2 * math.pi * k), rel=1e-12)
    assert formation_conductance_w_per_mk(D, Dw, float("inf"), k, alpha, t) == pytest.approx(2 * math.pi * k / ft)
    assert formation_conductance_w_per_mk(D, Dw, U, 0.0, alpha, t) == 0.0
    # the formation resists more as time goes on
    assert formation_conductance_w_per_mk(D, Dw, U, k, alpha, 10 * t) < C


# ── black-oil wells ──────────────────────────────────────────────────────────

def test_black_oil_tubing_satisfies_the_first_law():
    """ṁ·[h(P_res, T_res) − h(P_wh, T_wh)] = Q_formation + ṁ·g·Δz.
    Tolerance ṁ·cp·0.01 K (50–150 W); residuals obtained are printed by
    ``pytest -s``."""
    fluid = BlackOilFluid(**RICH)
    well = _well()
    for m_dot in (2.0, 4.0, 6.0):
        r = well_tubing_state(well, fluid, 250 * BAR, 85.0, m_dot)
        h_res = fluid.enthalpy_for_link(None, 250 * BAR, 85.0)
        h_wh = fluid.enthalpy_for_link(None, r.wellhead_pressure_pa, r.wellhead_temperature_c)
        residual = m_dot * (h_res - h_wh) - r.formation_heat_w - m_dot * G * 2000.0
        cp = fluid.pvt(r.wellhead_pressure_pa, r.wellhead_temperature_c).mixture_specific_heat_j_per_kg_k
        print(f"first law, m={m_dot}: residual {residual:.3f} W, tolerance {m_dot * cp * 0.01:.0f} W")
        assert abs(residual) < m_dot * cp * 0.01


def test_wellhead_temperature_trends():
    fluid = BlackOilFluid(**RICH)
    t_surface = 85.0 - 0.03 * 2000.0
    T = [well_tubing_state(_well(), fluid, 250 * BAR, 85.0, m).wellhead_temperature_c for m in (1.0, 2.5, 4.0, 6.0)]
    assert all(b > a for a, b in zip(T, T[1:]))            # faster flow arrives hotter
    assert all(t_surface < t < 95.0 for t in T)
    base = well_tubing_state(_well(), fluid, 250 * BAR, 85.0, 5.0)
    assert base.formation_heat_w > 0.0
    # insulating wellbore (finite U) and a longer production time both keep the heat in
    assert well_tubing_state(_well(tubing_u_w_per_m2k=10.0), fluid, 250 * BAR, 85.0, 5.0
                             ).wellhead_temperature_c > base.wellhead_temperature_c
    assert well_tubing_state(_well(production_time_days=3650.0), fluid, 250 * BAR, 85.0, 5.0
                             ).wellhead_temperature_c > base.wellhead_temperature_c
    # adiabatic tubing: no heat to the formation
    adiabatic = well_tubing_state(_well(formation_conductivity_w_per_mk=0.0), fluid, 250 * BAR, 85.0, 5.0)
    assert adiabatic.formation_heat_w == 0.0
    assert adiabatic.wellhead_temperature_c > base.wellhead_temperature_c


def test_very_low_rate_follows_the_formation():
    """Stiff limit (relaxation distance ≪ segment): the fluid arrives at the
    surface formation temperature and the march stays stable."""
    fluid = BlackOilFluid(api_gravity=30.0, gas_gravity=0.65, gor_sc_m3_per_m3=0.0, wor_sc_m3_per_m3=0.0)
    r = well_tubing_state(_well(), fluid, 250 * BAR, 85.0, 1e-4)
    assert r.wellhead_temperature_c == pytest.approx(85.0 - 0.03 * 2000.0, abs=0.05)


def test_cooler_tubing_weighs_more():
    """The deliverability curve uses the tubing temperature: with a steeper
    geothermal gradient the column is colder and denser, so the wellhead
    pressure at a given rate is lower."""
    fluid = BlackOilFluid(**RICH)
    warm = well_curve(_well(geothermal_gradient_k_per_m=0.0), fluid, 250 * BAR, 85.0)
    cold = well_curve(_well(geothermal_gradient_k_per_m=0.035), fluid, 250 * BAR, 85.0)
    assert cold.wellhead_pressures_pa[0] < warm.wellhead_pressures_pa[0]


def test_component_validation():
    assert _well().wellbore_diameter_m == pytest.approx(0.2159)
    assert _well(D=0.2).wellbore_diameter_m == pytest.approx(0.3)
    for bad in ({"wellbore_diameter_m": 0.05}, {"geothermal_gradient_k_per_m": -0.01},
                {"formation_conductivity_w_per_mk": -1.0}, {"production_time_days": 0.0},
                {"tubing_u_w_per_m2k": -5.0}):
        with pytest.raises(ValueError):
            _well(**bad)


# ── network ──────────────────────────────────────────────────────────────────

def _network(component) -> NetworkCase:
    return NetworkCase(
        name="well-thermal",
        fluid_model=BlackOilFluid(**RICH),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=250 * BAR),),
        pressure_outlets=(PressureBoundary(node_id=3, pressure_pa=15 * BAR),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=85.0, bc_type="fixed_temperature"),),
        inlet_fluid_bcs=(InletFluidBC(node_id=1, **RICH),),
        components=(component,
                    Pipe(start_node=2, end_node=3, diameter_m=0.1016, length_m=1500.0,
                         absolute_roughness_m=4.6e-5)),
    )


def test_well_delivers_its_wellhead_temperature_to_the_network():
    """Residual obtained: 2e-4 kW against ṁ·cp·0.01 K ≈ 0.16 kW."""
    well = _well()
    res = SteadyBlackOilSolver().solve(_network(well))
    assert res.converged
    m = res.component_flows[0].mass_flow_kg_per_s
    state = well_tubing_state(well, BlackOilFluid(**RICH), 250 * BAR, 85.0, m)
    # same enthalpy, at the network's wellhead pressure
    assert res.node_temperatures_c[2] == pytest.approx(state.wellhead_temperature_c, abs=0.2)
    assert res.node_temperatures_c[2] < 80.0
    geb = res.global_energy_balance
    assert not hasattr(geb, "isothermal_device_heat_kw")
    assert geb.heat_wall_loss_kw == pytest.approx(state.formation_heat_w / 1e3, rel=1e-3)
    assert geb.potential_energy_kw == pytest.approx(-m * G * 2000.0 / 1e3, rel=1e-6)
    assert abs(geb.energy_error_kw) < m * 2500.0 * 0.01 / 1e3


def test_ipr_link_is_isenthalpic():
    """A WellInflow link (reservoir → bottom hole) is an isenthalpic expansion:
    liquids warm slightly on the drawdown."""
    ipr = WellInflow(start_node=1, end_node=2, diameter_m=0.0762, productivity_index_m3_per_day_per_bar=6.0)
    res = SteadyBlackOilSolver().solve(_network(ipr))
    assert res.converged
    fluid = BlackOilFluid(**RICH)
    h_res = fluid.enthalpy_for_link(None, 250 * BAR, 85.0)
    h_bh = fluid.enthalpy_for_link(None, res.node_pressures_pa[2], res.node_temperatures_c[2])
    assert h_bh == pytest.approx(h_res, abs=2500.0 * 0.01)
    assert abs(res.global_energy_balance.energy_error_kw) < res.component_flows[0].mass_flow_kg_per_s * 25.0 / 1e3
