"""
Gas-hydrate formation conditions (van der Waals–Platteeuw, Munck constants).

Validation data and their standing
----------------------------------
Hydrate quadruple points as tabulated from Sloan & Koh, *Clathrate Hydrates
of Natural Gases* — SECONDARY SOURCE, seen only through a search-engine
summary of a ScienceDirect topic page (the table itself was not read):
    CH4  Q1 272.9 K, 2.563 MPa
    C2H6 Q1 273.1 K, 0.530 MPa;  Q2 287.8 K, 3.390 MPa
    C3H8 Q1 273.1 K, 0.172 MPa;  Q2 278.8 K, 0.556 MPa
    CO2  Q1 273 K,   1.256 MPa   (other sources give 1.04 MPa at 271.6 K)
Points at 273.1 K are compared with the liquid-water branch at 273.15 K.
Q2 lies where the Lw-H-V curve turns almost vertical (the former liquefies),
so it is compared in temperature at the measured pressure.

Tolerances
----------
Pressure: 11 %, the accuracy Antunes et al. (2018, Braz. J. Chem. Eng. 35,
265 — READ IN SOURCE) report for single-guest systems with a vdW-P model and
Munck's constants.  Temperature at Q2: 1 K (≈ 11 % in pressure on this part
of the curve; the roadmap's criterion).  These tolerances were not chosen
from the results below.

Methane solubility in water — READ IN SOURCE: Böttger, Kamps & Maurer
(2016), Fluid Phase Equilib. 407, 209, from the NIST ThermoML archive
(trc.nist.gov/ThermoML/10.1016/j.fluid.2015.03.041), with the expanded
(95 %) uncertainty of each point.

Hydrate temperatures at 10–100 MPa — SECONDARY SOURCE: Sloan & Koh values
as quoted by Algaba et al. (2024), arXiv:2408.01819, Table III.

Katz chart — CROSS-CHECK, not validation: Towler-Mokhatab temperatures as
tabulated by Janna (2015), Tables 4.1b–4.5b (read in the dissertation).

NOT VALIDATED (no measured data accessible): pure methane between 285 K and
295 K and natural-gas mixtures.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.properties.hydrate import (  # noqa: E402
    hydrate_curve,
    hydrate_formation_pressure,
    hydrate_formation_temperature,
)

P_TOL = 0.11       # Antunes et al. 2018, single-guest vdW-P + Munck accuracy
T_TOL_K = 1.0


# Results obtained (v1.19.1, with Henry's-law water activity):
#   methane  2.751 MPa (+7.3 %, vs Q1 at 272.9 K — part of it is the 0.25 K offset)
#   ethane   0.496 MPa (−6.5 %)    propane 0.166 MPa (−3.6 %)    CO2 1.240 MPa (−1.2 %)
@pytest.mark.parametrize("gas, p_q1_mpa, structure", [
    ("methane", 2.563, "sI"), ("ethane", 0.530, "sI"), ("propane", 0.172, "sII"),
    ("carbon dioxide", 1.256, "sI"),
])
def test_lower_quadruple_points(gas, p_q1_mpa, structure):
    pt = hydrate_formation_pressure([gas], [1.0], 0.0)
    assert pt.structure == structure
    assert pt.pressure_pa / 1e6 == pytest.approx(p_q1_mpa, rel=P_TOL)


# Results obtained: ethane 286.98 K (−0.82 K), propane 278.52 K (−0.28 K)
@pytest.mark.parametrize("gas, t_q2_k, p_q2_mpa", [("ethane", 287.8, 3.390), ("propane", 278.8, 0.556)])
def test_upper_quadruple_points(gas, t_q2_k, p_q2_mpa):
    pt = hydrate_formation_temperature([gas], [1.0], p_q2_mpa * 1e6)
    assert pt.temperature_c + 273.15 == pytest.approx(t_q2_k, abs=T_TOL_K)


def test_dissolved_gas_lowers_water_activity():
    """Consistency (not validation): CO2 dissolves ~25× more than methane, so
    its water activity drops more and its formation pressure rises."""
    from angelica.properties.hydrate import water_mole_fraction
    x_co2 = water_mole_fraction(278.15, {"124-38-9": 20.0})
    x_c1 = water_mole_fraction(278.15, {"74-82-8": 20.0})
    assert 0.97 < x_co2 < 0.985 and x_c1 > 0.999
    assert water_mole_fraction(278.15, {}) == 1.0


# (T K, P kPa, molality mol/kg, expanded uncertainty mol/kg)
_BOTTGER_2016 = [
    (283.14, 1151, 0.0214, 0.0043), (283.14, 1967, 0.0352, 0.0044), (283.14, 2456, 0.0447, 0.0044),
    (283.14, 3342, 0.0576, 0.0044), (283.14, 3876, 0.0669, 0.0044), (283.14, 4673, 0.0784, 0.0044),
    (283.14, 5701, 0.0920, 0.0044), (298.15, 3200, 0.0424, 0.0043), (298.15, 4111, 0.0527, 0.0043),
    (298.15, 4989, 0.0638, 0.0043), (298.15, 6539, 0.0812, 0.0043), (298.15, 8675, 0.0993, 0.0043),
    (298.15, 10360, 0.1141, 0.0043),
]


@pytest.mark.parametrize("T, p_kpa, molality, u", _BOTTGER_2016)
def test_methane_solubility_against_measurements(T, p_kpa, molality, u):
    """Validation (data read in source).  Errors obtained: −4.1 … +3.5 %
    (largest in absolute terms +0.0040 mol/kg at 10.36 MPa).  Tolerance: the
    expanded uncertainty reported with each point."""
    from angelica.properties.hydrate import _M_WATER, _former_fugacities, water_mole_fraction
    fug = _former_fugacities(("methane",), (1.0,), T, p_kpa * 1e3, "PR")
    x_w = water_mole_fraction(T, fug)
    m = (1.0 - x_w) / (x_w * _M_WATER)
    assert abs(m - molality) <= u


# Errors obtained: methane 40 MPa −1.56 K; CO2 10 MPa +0.36 K, 40 MPa −1.35 K.
# Tolerance 2 K: a regression guard chosen after seeing the results — the
# data come through a secondary source and no accuracy is published for
# vdW-P + Munck in this range, so this is NOT a validation.
@pytest.mark.parametrize("gas, p_mpa, t_k", [
    ("methane", 40.0, 297.0), ("carbon dioxide", 10.0, 283.6), ("carbon dioxide", 40.0, 286.2),
])
def test_high_pressure_hydrate_temperature_guard(gas, p_mpa, t_k):
    pt = hydrate_formation_temperature([gas], [1.0], p_mpa * 1e6)
    assert pt.temperature_c + 273.15 == pytest.approx(t_k, abs=2.0)


@pytest.mark.xfail(strict=True, reason="Known limitation: Henry's law without a Poynting factor "
                   "over-predicts dissolved CO2 at 100 MPa; model 280.8 K vs 289.7 K")
def test_co2_hydrate_at_100_mpa_known_limitation():
    pt = hydrate_formation_temperature(["carbon dioxide"], [1.0], 100e6)
    assert pt.temperature_c + 273.15 == pytest.approx(289.7, abs=2.0)


# Towler-Mokhatab temperature (°F) at 500 psia, Janna (2015) Tables 4.1b–4.5b.
# Differences obtained (chart − correlation): −1.45, −0.13, −0.77, −1.11,
# −1.81 °F.  Tolerance 2 °F: regression guard chosen after the results; a
# cross-check of two representations of gas-gravity data, not a validation
# against measurements.  At 2000 psia the two differ by 2.4 °F (gravity 0.6)
# to 8.5 °F (gravity 1.0): the log-linear correlation does not follow the
# steepening chart curve.
@pytest.mark.parametrize("gravity, t_f", [(0.6, 51.17), (0.7, 54.85), (0.8, 58.04), (0.9, 60.85), (1.0, 63.36)])
def test_katz_chart_against_towler_mokhatab(gravity, t_f):
    from angelica.properties.hydrate import katz_hydrate_temperature
    pt = katz_hydrate_temperature(gravity, 500.0 * 6894.757)
    assert pt.temperature_c * 1.8 + 32.0 == pytest.approx(t_f, abs=2.0)


def test_propane_promotes_structure_ii_in_natural_gas():
    lean = (["methane", "ethane", "propane", "n-butane"], [0.90, 0.06, 0.03, 0.01])
    p_ng = hydrate_formation_pressure(*lean, 10.0)
    p_c1 = hydrate_formation_pressure(["methane"], [1.0], 10.0)
    assert p_ng.structure == "sII" and p_c1.structure == "sI"
    assert p_ng.pressure_pa < 0.5 * p_c1.pressure_pa


def test_curve_monotonic_and_inverse_consistent():
    names, z = ["methane", "ethane", "propane"], [0.85, 0.10, 0.05]
    curve = hydrate_curve(names, z, 0.0, 20.0, 9)
    p = [pt.pressure_pa for pt in curve]
    assert all(b > a for a, b in zip(p, p[1:]))
    mid = curve[4]
    back = hydrate_formation_temperature(names, z, mid.pressure_pa)
    assert back.temperature_c == pytest.approx(mid.temperature_c, abs=0.02)


def test_water_and_non_formers():
    # water in the composition is ignored (dry basis, water in excess)
    dry = hydrate_formation_pressure(["methane", "propane"], [0.97, 0.03], 5.0)
    wet = hydrate_formation_pressure(["methane", "propane", "water"], [0.97 * 0.9, 0.03 * 0.9, 0.1], 5.0)
    assert wet.pressure_pa == pytest.approx(dry.pressure_pa, rel=1e-6)
    # n-pentane alone is not a hydrate former here
    assert hydrate_formation_pressure(["n-pentane"], [1.0], 5.0) is None
    # no hydrate above 0 °C at 1 bar for methane
    assert hydrate_formation_temperature(["methane"], [1.0], 1e5) is None
    with pytest.raises(ValueError):
        hydrate_formation_pressure(["methane"], [1.0], -5.0)


def test_compositional_results_report_hydrate_margin(tmp_path):
    """Bare subsea line of tutorial 08: hydrate temperature and subcooling per
    pipe, consistent with the stand-alone function, and exported to CSV."""
    import importlib.util

    from angelica.io.reporting import export_solve_result_csv
    from angelica.solvers import SteadyCompositionalSolver

    spec = importlib.util.spec_from_file_location(
        "tutorial_compositional_08",
        SRC_ROOT / "angelica/tutorials/steady_compositional/08_hydrate_risk_subsea/run.py")
    t08 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(t08)
    case = t08.build_case(15.0)
    res = SteadyCompositionalSolver().solve(case)
    assert res.converged
    last = res.component_flows[-1]
    P_out = res.node_pressures_pa[len(res.component_flows) + 1]
    T_out = res.node_temperatures_c[len(res.component_flows) + 1]
    ref = hydrate_formation_temperature(t08.COMPONENTS, t08.ZS, P_out)
    assert last.hydrate_temperature_c == pytest.approx(ref.temperature_c, abs=1e-6)
    assert last.hydrate_subcooling_k == pytest.approx(ref.temperature_c - T_out, abs=1e-6)
    assert last.hydrate_subcooling_k > 5.0 and last.hydrate_structure == "sII"
    assert res.component_flows[0].hydrate_subcooling_k < 0.0          # warm inlet section
    out = tmp_path / "r.csv"
    export_solve_result_csv(res, str(out))
    assert "Hydrate subcooling (K)" in out.read_text(encoding="utf-8")


# Semenov et al. (2024), Data in Brief 53, 110138, Table 1, system 1 (pure
# water) — READ IN SOURCE (Europe PMC full text, PMC10877679): (T K, P MPa).
SEMENOV_2024_WATER = [(274.25, 2.891), (276.96, 3.785), (279.26, 4.771), (281.11, 5.777), (282.62, 6.77),
                      (283.87, 7.776), (285.01, 8.839), (286.14, 10.027), (287.03, 11.143), (287.65, 12.009),
                      (287.85, 12.316), (289.02, 14.222)]


def test_methane_against_semenov_2024():
    """Pure methane, 274–289 K: the model is colder than the measurements by
    0.60 K at 2.9 MPa growing to 0.91 K at 14.2 MPa (the warm-end deviation
    of the module docstring).  Guard 1.0 K: set after the results, not an
    accuracy claim."""
    errors = [hydrate_formation_temperature(["methane"], [1.0], p * 1e6).temperature_c + 273.15 - t
              for t, p in SEMENOV_2024_WATER]
    assert all(-1.0 < e < -0.55 for e in errors), errors
