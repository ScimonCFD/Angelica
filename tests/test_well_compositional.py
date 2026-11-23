"""
``Well`` and ``WellInflow`` with compositional fluids (roadmap item 4.3).

Standing of the checks
----------------------
- PVT table: the tubing march on the interpolated table is compared with the
  same march flashing the EOS at every evaluation.  Errors obtained: P_wh
  within 0.05 %, T_wh within 0.01 K (gas and condensate wells, 10–90 % of the
  open flow).  Tolerances 0.2 % and 0.1 K: the table must stay well below the
  ±0.3 % P_wh scatter already caused by the tubing segmentation (ROADMAP,
  stage 4 limitations).  A first, geometric-only pressure grid (55 bar
  between nodes at 250 bar) gave 0.45 % and 2.9 kJ/kg enthalpy errors in
  dense gas; the grid now has at most 10 bar between nodes.
- Single-phase gas: tubing pressure against an independent integration
  (hydrostatic + Haaland friction with EOS density and viscosity).  Error
  obtained 0.045 %; tolerance 0.5 % for the different friction-factor
  formulas (Churchill in Beggs-Brill, Haaland here), as in the black-oil test.
- Standard volumes: mass per stock-tank m³ against the black-oil formula
  ρ_o + GOR·ρ_g with the same GOR and densities.  Error obtained 0.1 %;
  tolerance 0.5 %: black-oil takes the gas as ideal at standard conditions
  (Z ≈ 0.997 for these gases).
- First law: network balance closes (by construction for a table fluid: the
  well's heat is taken from the EOS enthalpies at both ends); the independent
  check is that this heat agrees with the table march's own heat
  (obtained 0.07 %).
- NOT VALIDATED against a published gas-well worked example (none
  accessible).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures import well_model as wm  # noqa: E402
from angelica.closures.well_fluid import CompositionalWellFluid, well_fluid_for  # noqa: E402
from angelica.closures.well_inflow import IPR, ipr_of, ipr_rate  # noqa: E402
from angelica.core.case import InletCompositionBC, NetworkCase, PressureBoundary, ThermalBoundary  # noqa: E402
from angelica.core.components import Pipe, Well, WellInflow  # noqa: E402
from angelica.properties import compositional_fluid as cf  # noqa: E402
from angelica.properties.compositional_fluid import CompositionalFluid  # noqa: E402
from angelica.solvers import SteadyCompositionalSolver  # noqa: E402

BAR = 1e5
NAMES = ("methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane", "n-decane", "nitrogen",
         "carbon dioxide")
GAS = (0.86, 0.06, 0.03, 0.015, 0.005, 0.005, 0.005, 0.01, 0.01)
CONDENSATE = (0.74, 0.08, 0.05, 0.03, 0.02, 0.02, 0.04, 0.01, 0.01)


def _gas_well(**kw) -> Well:
    args = {"component_id": "gw", "start_node": 1, "end_node": 2, "diameter_m": 0.0762,
            "ipr_model": "backpressure", "ipr_coefficient": 150.0, "ipr_exponent": 0.8,
            "tubing_length_m": 3000.0, "tubing_height_change_m": 3000.0, "tubing_correlation": "gray", **kw}
    return Well(**args)


class _Direct(CompositionalWellFluid):
    """Reference adapter: an EOS flash at every evaluation."""

    def phase_properties(self, P, T):
        return cf._beggs_brill_properties(self.component_names, round(P, -1), round(T, 3), self.zs, self.eos_name)

    def enthalpy(self, P, T):
        return cf._mixture_enthalpy(self.component_names, round(P, -1), round(T, 4), self.zs, self.eos_name)

    def cp(self, P, T):
        return (self.enthalpy(P, T + 0.25) - self.enthalpy(P, T - 0.25)) / 0.5

    def dh_dp(self, P, T):
        return (self.enthalpy(P + 2e4, T) - self.enthalpy(P - 2e4, T)) / 4e4

    def temperature_from_h(self, h, P, T0):
        T = T0
        for _ in range(8):
            f = self.enthalpy(P, T) - h
            if abs(f) < 0.05:
                break
            T -= f / self.cp(P, T)
        return T


@pytest.mark.parametrize("zs, fraction", [(GAS, 0.5), (GAS, 0.9), (CONDENSATE, 0.1)])
def test_pvt_table_matches_direct_flash(zs, fraction):
    fluid = CompositionalFluid(NAMES, zs)
    well = _gas_well()
    table = well_fluid_for(fluid, None, well, 250 * BAR, 100.0)
    ref = _Direct(table.component_names, table.zs, table.eos_name, table.p_max_pa, table.t_min_c, table.t_max_c)
    ipr = ipr_of(well)
    q = fraction * wm.well_curve(well, fluid, 250 * BAR, 100.0).oil_rates_sc_m3_per_day[-1]
    args = (250 * BAR, 100.0, q, ipr, 250 * BAR, table.mass_per_base_unit_for(ipr), 0.0762, 3000.0, 3000.0,
            4.6e-5, 25, "gray", wm._thermal_key(well), True)
    a, b = wm._march(table, *args), wm._march(ref, *args)
    assert a.wellhead_pressure_pa == pytest.approx(b.wellhead_pressure_pa, rel=2e-3)
    assert a.wellhead_temperature_c == pytest.approx(b.wellhead_temperature_c, abs=0.1)


def test_dry_gas_tubing_matches_single_phase_integration():
    names, zs = ("methane", "ethane", "nitrogen"), (0.95, 0.03, 0.02)
    fluid = CompositionalFluid(names, zs)
    D, L, eps, T, m_dot, p_wf = 0.0762, 3000.0, 4.6e-5, 80.0, 6.0, 200 * BAR
    wf = well_fluid_for(fluid, None, _gas_well(), 250 * BAR, 90.0)
    p_wh = wm.tubing_wellhead_pressure(p_wf, m_dot, wf, T, D, L, L, eps, 200)
    area, P, n = math.pi * D * D / 4.0, p_wf, 3000
    for _ in range(n):
        rho, mu = cf._flash_properties(names, round(P, -1), T, zs, "PR")[:2]
        v = m_dot / (rho * area)
        f = (-1.8 * math.log10((eps / D / 3.7) ** 1.11 + 6.9 / (rho * v * D / mu))) ** -2   # Haaland
        P -= L / n * (rho * 9.80665 + f * rho * v * v / (2 * D))
    assert wf.phase_properties(150 * BAR, T)[0] == 0.0          # no liquid anywhere
    assert p_wh == pytest.approx(P, rel=5e-3)


def test_standard_volumes_agree_with_black_oil_formula():
    """Volatile oil: mass per stock-tank m³ from the EOS stream equals
    ρ_o + GOR·ρ_g,sc with the stream's own GOR, oil density and gas gravity."""
    from angelica.properties.compositional_fluid import _hc_phase_split

    names = ("methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane", "n-decane", "n-hexadecane")
    zs = (0.45, 0.08, 0.06, 0.04, 0.03, 0.04, 0.20, 0.10)
    stream = CompositionalWellFluid(names, zs, "PR")
    _m, v_gas, v_oil = stream._standard_volumes
    rho_o = stream.stock_tank_oil_density_kg_per_m3
    mw_g = _hc_phase_split(names, 101_325.0, 15.56, zs, "PR")[3]
    black_oil = rho_o + v_gas / v_oil * (mw_g / 28.97 * 1.2250)
    assert stream.mass_per_rate(IPR("linear", productivity_index=1.0)) == pytest.approx(black_oil, rel=5e-3)
    # Peneloux-translated stock-tank oil: n-C16-rich liquid near 760 kg/m³
    # (untranslated PR gives 671)
    assert 740.0 < rho_o < 790.0


def test_ipr_rate_basis_checks():
    dry = CompositionalWellFluid(("methane", "ethane"), (0.97, 0.03), "PR")
    with pytest.raises(ValueError, match="no liquid at standard conditions"):
        dry.mass_per_rate(IPR("fetkovich", coefficient=1.0))
    assert dry.mass_per_rate(IPR("backpressure", coefficient=1.0)) == pytest.approx(0.72, rel=0.05)
    assert math.isinf(dry.bubble_point_pa(200 * BAR, 80.0))
    oil = CompositionalWellFluid(("methane", "n-decane"), (0.4, 0.6), "PR")
    assert 5 * BAR < oil.bubble_point_pa(300 * BAR, 90.0) < 300 * BAR


def _gas_network(well, p_res=250 * BAR, p_sep=50 * BAR, U=5.0, zs=GAS) -> NetworkCase:
    return NetworkCase(
        name="gas well", fluid_model=CompositionalFluid(NAMES, zs),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=p_res),),
        pressure_outlets=(PressureBoundary(node_id=3, pressure_pa=p_sep),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=100.0, bc_type="fixed_temperature"),),
        inlet_composition_bcs=(InletCompositionBC(node_id=1, zs=zs),),
        components=(well, Pipe(component_id="fl", start_node=2, end_node=3, diameter_m=0.1524, length_m=5000.0,
                               absolute_roughness_m=4.6e-5, heat_transfer_coefficient_w_per_m2k=U,
                               ambient_temperature_c=15.0)))


def test_gas_well_in_a_network():
    """Operating point on the curve; first law closed; the heat from the
    table march agrees with the heat from the EOS enthalpies at both ends
    (obtained 0.07 %; tolerance 0.5 %: enthalpy interpolation errors up to
    0.2 kJ/kg against ~20 kJ/kg lost per kg of gas)."""
    well = _gas_well()
    res = SteadyCompositionalSolver().solve(_gas_network(well, U=0.0))     # adiabatic flowline
    assert res.converged
    m = res.component_flows[0].mass_flow_kg_per_s
    assert res.well_states[0].state == "producing" and m > 3.0
    fluid = CompositionalFluid(NAMES, GAS)
    curve = wm.well_curve(well, fluid, 250 * BAR, 100.0)
    assert m == pytest.approx(curve.mass_flow(res.node_pressures_pa[2])[0], rel=1e-3)
    geb = res.global_energy_balance
    assert abs(geb.energy_error_kw) < m * 2500.0 * 0.01 / 1e3              # ṁ·cp·0.01 K
    state = wm.well_tubing_state(well, fluid, 250 * BAR, 100.0, m)
    assert geb.heat_wall_loss_kw == pytest.approx(state.formation_heat_w / 1e3, rel=5e-3)
    assert geb.potential_energy_kw == pytest.approx(-m * 9.80665 * 3000.0 / 1e3, rel=1e-6)
    # the gas cools on the drawdown (Joule-Thomson), unlike a black-oil liquid
    assert state.bottomhole_temperature_c < 100.0


def test_gas_well_dead_by_back_pressure():
    res = SteadyCompositionalSolver().solve(_gas_network(_gas_well(), p_sep=200 * BAR))
    assert res.converged
    assert res.well_states[0].state == "dead_back_pressure"
    assert res.component_flows[0].mass_flow_kg_per_s == 0.0


def test_ipr_only_link_with_compositional_fluid():
    """WellInflow (no tubing): the solver reproduces the back-pressure law,
    ṁ = C·(P_r² − P_wf²)ⁿ · (kg per standard m³ of gas)."""
    ipr_link = WellInflow(component_id="ipr", start_node=1, end_node=2, diameter_m=0.0762,
                          ipr_model="backpressure", ipr_coefficient=150.0, ipr_exponent=0.8)
    res = SteadyCompositionalSolver().solve(_gas_network(ipr_link, p_res=150 * BAR, p_sep=100 * BAR))
    assert res.converged
    m = res.component_flows[0].mass_flow_kg_per_s
    stream = CompositionalWellFluid(NAMES, GAS, "PR")
    ipr = IPR("backpressure", coefficient=150.0, exponent=0.8)
    expected = ipr_rate(ipr, 150 * BAR, res.node_pressures_pa[2]) * stream.mass_per_rate(ipr) / 86_400.0
    assert m == pytest.approx(expected, rel=1e-3)
