"""
Gas-hydrate formation conditions (van der Waals–Platteeuw, Munck constants), v1.19.0.

Validation data: hydrate quadruple points as tabulated from Sloan & Koh,
*Clathrate Hydrates of Natural Gases* (Q1: ice–water–hydrate–gas; Q2:
water–liquid hydrocarbon–gas–hydrate):
    CH4  Q1 272.9 K, 2.563 MPa
    C2H6 Q1 273.1 K, 0.530 MPa;  Q2 287.8 K, 3.390 MPa
    C3H8 Q1 273.1 K, 0.172 MPa;  Q2 278.8 K, 0.556 MPa
Points at 273.1 K are compared with the liquid-water branch at 273.15 K.
Q2 lies where the Lw-H-V curve turns almost vertical (the former liquefies),
so it is compared in temperature at the measured pressure.
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


@pytest.mark.parametrize("gas, p_q1_mpa, structure", [
    ("methane", 2.563, "sI"), ("ethane", 0.530, "sI"), ("propane", 0.172, "sII"),
])
def test_lower_quadruple_points(gas, p_q1_mpa, structure):
    pt = hydrate_formation_pressure([gas], [1.0], 0.0)
    assert pt.structure == structure
    # within 8 % in pressure (≈ 0.5 K on the Lw-H-V curve near 273 K)
    assert pt.pressure_pa / 1e6 == pytest.approx(p_q1_mpa, rel=0.08)


@pytest.mark.parametrize("gas, t_q2_k, p_q2_mpa", [("ethane", 287.8, 3.390), ("propane", 278.8, 0.556)])
def test_upper_quadruple_points_within_1_kelvin(gas, t_q2_k, p_q2_mpa):
    pt = hydrate_formation_temperature([gas], [1.0], p_q2_mpa * 1e6)
    assert pt.temperature_c + 273.15 == pytest.approx(t_q2_k, abs=1.0)


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
