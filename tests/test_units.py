"""
Oil-field units and standard-volume rates (roadmap item 7.1, v1.31.0).

Validation
----------
Conversion factors — READ IN SOURCE: NIST SP 811 (2008 edition), Appendix
B.8, nist.gov (retrieved 2026-11-27), printed to 7 significant figures:
barrel (42 U.S. gal) 1.589 873 E-01 m³, gallon (U.S.) 3.785 412 E-03 m³,
cubic foot 2.831 685 E-02 m³, foot 3.048 E-01 m, inch 2.54 E-02 m,
pound (avdp) 4.535 924 E-01 kg, psi 6.894 757 E+03 Pa, centipoise
1.0 E-03 Pa·s.  Angelica's factors are the exact definitions; tolerance:
half a unit in the 7th printed digit (the table's rounding).

Standard-volume rates — BY DEFINITION in black-oil: q_o = ṁ/(ρ_o,sc +
GOR·ρ_g,sc + WOR·ρ_w,sc), q_g = GOR·q_o, q_w = WOR·q_o, recomputed here
independently for every link (the first version took gas from the split at
standard conditions and lost Standing's residual Rs there: −0.7 % of gas in
tutorial 14 — found in review).  CONSISTENCY: the stock-tank volumes do not
depend on the separation path, so the stream feeding a separator carries
the separator's own rates (black-oil tutorial 14):
obtained oil 705.49 = 705.49 Sm³/d (≤ 1e-6), water 102.1589 against 102.1587 m³/d
(1.6e-6);
gas through the separator conserved within 0.11 % (the gas outlet's
black-oil representation, GOR 10⁶, carries 1e-6 Sm³ oil per Sm³ gas).
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from angelica import units  # noqa: E402

NIST = [  # (Angelica factor, NIST SP 811 printed value, last printed digit's unit)
    (units.BBL_M3, 1.589873e-01, 1e-7),
    (units.GAL_M3, 3.785412e-03, 1e-9),
    (units.FT3_M3, 2.831685e-02, 1e-8),
    (units.FT_M, 3.048e-01, 1e-4),
    (units.IN_M, 2.54e-02, 1e-4),
    (units.LB_KG, 4.535924e-01, 1e-7),
    (units.PSI_PA, 6.894757e03, 1e-3),
    (units.CP_PA_S, 1.0e-03, 1e-4),
]


@pytest.mark.parametrize("value, nist, ulp", NIST)
def test_conversion_factors_against_nist_sp811(value, nist, ulp):
    assert abs(value - nist) <= 0.5 * ulp


def test_round_trips_and_fixed_points():
    assert units.c_to_f(100.0) == pytest.approx(212.0) and units.c_to_f(-40.0) == pytest.approx(-40.0)
    assert units.f_to_c(units.c_to_f(37.3)) == pytest.approx(37.3)
    assert units.pa_to_psig(units.ATM_PA) == pytest.approx(0.0, abs=1e-12)
    assert units.psig_to_pa(units.pa_to_psig(5.0e6)) == pytest.approx(5.0e6)
    assert units.ATM_PSI == pytest.approx(14.69595, abs=1e-5)
    assert units.sm3_to_bbl(units.BBL_M3) == pytest.approx(1.0)
    assert units.sm3_to_mmscf(1e6 * units.FT3_M3) == pytest.approx(1.0)


@pytest.fixture(scope="module")
def tutorial14():
    from angelica.closures import BeggsBrillCorrelation
    from angelica.solvers import SteadyBlackOilSolver
    spec = importlib.util.spec_from_file_location(
        "t14", ROOT / "src/angelica/tutorials/steady_black_oil/14_separator/run.py")
    t14 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(t14)
    return SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation()).solve(t14.build_case())


def test_black_oil_standard_rates_match_the_separator(tutorial14):
    res = tutorial14
    flows = {f.label.split(":")[-1]: f for f in res.component_flows}
    sep = res.separators[0]
    trunk = flows["trunk"]
    assert trunk.std_oil_sm3_per_day == pytest.approx(sep.oil_rate_sm3_per_day, rel=1e-6)
    assert trunk.std_water_sm3_per_day == pytest.approx(sep.water_rate_m3_per_day, rel=1e-5)   # obtained 1.6e-6
    out_gas = flows["gas_line"].std_gas_sm3_per_day + flows["oil_line"].std_gas_sm3_per_day
    assert out_gas == pytest.approx(trunk.std_gas_sm3_per_day, rel=2e-3)
    # the wells' streams add up at the manifold (standard volumes are conserved by the mixing)
    wells = flows["flowline_Well1"].std_oil_sm3_per_day + flows["flowline_Well2"].std_oil_sm3_per_day
    # within the hydraulic criterion on the mass flows (2e-3); obtained 3.6e-6
    assert wells == pytest.approx(trunk.std_oil_sm3_per_day, rel=2e-3)


def test_black_oil_rates_by_definition(tutorial14):
    """Every link's rates against the definition, written independently."""
    from angelica.properties.black_oil import _AIR_DENSITY_SC_KG_M3
    from angelica.properties.dead_oil import dead_oil_density_kg_per_m3
    spec = importlib.util.spec_from_file_location(
        "t14b", ROOT / "src/angelica/tutorials/steady_black_oil/14_separator/run.py")
    t14 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(t14)
    res = tutorial14
    flows = {f.label.split(":")[-1]: f for f in res.component_flows}
    for i, (name, _p, _t, api, gg, gor, wor, *_rest) in enumerate(t14.WELLS):
        f = flows[f"flowline_{name.replace(' ', '')}"]
        q_o = f.mass_flow_kg_per_s * 86_400.0 / (dead_oil_density_kg_per_m3(api)
                                                 + gor * gg * _AIR_DENSITY_SC_KG_M3 + wor * 1025.0)
        assert f.std_oil_sm3_per_day == pytest.approx(q_o, rel=1e-12)
        assert f.std_gas_sm3_per_day == pytest.approx(gor * q_o, rel=1e-12)     # the GOR, exactly
        assert f.std_water_sm3_per_day == pytest.approx(wor * q_o, rel=1e-12)
    gas_in = flows["flowline_Well1"].std_gas_sm3_per_day + flows["flowline_Well2"].std_gas_sm3_per_day
    assert gas_in == pytest.approx(flows["trunk"].std_gas_sm3_per_day, rel=2e-3)   # mixing conserves it


def test_compositional_water_counts_all_the_water():
    """Condensate with water: the water rate is all the stream's water as
    liquid at standard conditions, also what the flash leaves in the gas
    (review: the first version dropped 24 % of it)."""
    from types import SimpleNamespace as NS

    from angelica.closures.well_fluid import T_SC_C, _molar_masses
    from angelica.properties.compositional_fluid import CompositionalFluid
    from angelica.properties.free_water import water_liquid_density_kg_m3
    from angelica.solvers.base import BaseSolver
    names = ("methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane", "water")
    zs = (0.55, 0.10, 0.10, 0.08, 0.05, 0.05, 0.07)
    fluid = CompositionalFluid(list(names), list(zs))
    link = NS(mass_flow_kg_per_s=10.0, zs=zs, component=None)
    rates = BaseSolver._standard_rates(link, fluid)
    mws = _molar_masses(names)
    x_w = zs[6] * mws[6] / sum(z * m for z, m in zip(zs, mws))
    assert rates["std_water_sm3_per_day"] == pytest.approx(
        10.0 * 86_400.0 * x_w / water_liquid_density_kg_m3(T_SC_C + 273.15), rel=1e-12)
    assert rates["std_gas_sm3_per_day"] > 0.0   # (all vapour at 1 atm: 5 % C5 and C6 in a lean gas)


def test_compositional_dry_gas_has_no_liquid_at_standard_conditions():
    from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
    from angelica.core.components import Pipe
    from angelica.properties.compositional_fluid import CompositionalFluid
    from angelica.solvers import SteadyCompositionalSolver
    fluid = CompositionalFluid(["methane", "ethane", "nitrogen"], [0.93, 0.05, 0.02])
    case = NetworkCase(name="g", fluid_model=fluid, pressure_inlets=(PressureBoundary(1, 50e5),),
                       pressure_outlets=(PressureBoundary(2, 30e5),), thermal_inlets=(ThermalBoundary(1, 20.0),),
                       components=(Pipe(component_id="p", start_node=1, end_node=2, diameter_m=0.2,
                                        length_m=5000.0, absolute_roughness_m=46e-6),))
    res = SteadyCompositionalSolver().solve(case)
    f = res.component_flows[0]
    assert f.std_oil_sm3_per_day == pytest.approx(0.0, abs=1e-9) and f.std_water_sm3_per_day == 0.0
    # The gas volume is the EOS's real-gas volume at 60 °F, 1 atm: compare
    # with an independent thermo Peng-Robinson flash of the same gas.
    # Tolerance 2e-3: the Peneloux translation applied to the standard
    # volumes, not to this plain PR flash.
    from thermo import PRMIX, CEOSGas, CEOSLiquid, ChemicalConstantsPackage, FlashVL

    from angelica.closures.well_fluid import P_SC_PA, T_SC_C
    consts, props = ChemicalConstantsPackage.from_IDs(["methane", "ethane", "nitrogen"])
    kw = {"eos_kwargs": {"Tcs": consts.Tcs, "Pcs": consts.Pcs, "omegas": consts.omegas},
          "HeatCapacityGases": props.HeatCapacityGases}
    flash = FlashVL(consts, props, liquid=CEOSLiquid(PRMIX, **kw), gas=CEOSGas(PRMIX, **kw))
    st = flash.flash(T=T_SC_C + 273.15, P=P_SC_PA, zs=[0.93, 0.05, 0.02])
    rho_pr = st.rho_mass()
    rho_std = f.mass_flow_kg_per_s * 86_400.0 / f.std_gas_sm3_per_day
    assert rho_std == pytest.approx(rho_pr, rel=2e-3)      # obtained −2.2e-4


def test_exports_carry_standard_rates(tmp_path, tutorial14):
    from angelica.io.reporting import export_solve_result_csv
    export_solve_result_csv(tutorial14, tmp_path / "r.csv")
    text = (tmp_path / "r.csv").read_text(encoding="utf-8")
    assert "STB/d" in text and "MMscf/d" in text


def test_gui_oil_field_units_convert_temperature_and_pressure():
    pytest.importorskip("sv_ttk")                    # the GUI's optional dependency
    from angelica.gui.app import NetSimGui
    me = SimpleNamespace(UNIT_SYSTEMS=NetSimGui.UNIT_SYSTEMS, _unit_system_key="field")
    me._unit_quantities = lambda: NetSimGui._unit_quantities(me)
    assert NetSimGui._si_to_display(me, "100.0", "temperature") == "212"
    assert float(NetSimGui._display_to_si(me, "212", "temperature")) == pytest.approx(100.0)
    assert float(NetSimGui._si_to_display(me, str(units.PSI_PA * 1000.0), "pressure")) == pytest.approx(1000.0)
    me._unit_system_key = "si"
    assert NetSimGui._si_to_display(me, "100.0", "temperature") == "100"
