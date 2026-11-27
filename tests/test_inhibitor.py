"""
Thermodynamic hydrate inhibitors, methanol and MEG (roadmap item 2b.3, v1.37.0).

Data — READ IN SOURCE (full text through Europe PMC, 2026-11-27):
  [S23] Semenov, Gong, Medvedev, Stoporev, Istomin, Vinokurov, Li (2023),
        Data in Brief 46, 108892, Table 1 (PMC9876827): methane hydrate
        V–Lw–H points with 10–50 mass % methanol, MEG, DEG, 4.9–8.4 MPa,
        and the depression ΔT_h against pure water.
  [S24] Semenov, Mendgaziev, Istomin, Sergeeva, Vinokurov, Gong, Li,
        Stoporev (2024), Data in Brief 53, 110138, Tables 1 and 10
        (PMC10877679): methane hydrate with 5–60 mass % methanol, 3–13 MPa,
        new samples; pure water; ice points of the methanol solutions.

Model: the water activity of the inhibited water, R·T·ln γ_w = x²(A + B·x),
fitted through the vdW-P model to the mean ΔT_h of [S23] at 6.5 MPa
(methanol: 10–50 %; MEG: 10, 30, 50 %).  The depression is compared as the
model's own pure-water temperature minus its inhibited one, so the pure-
water offset of the vdW-P model (below) does not enter.

  - Methanol against [S24]: new samples FROM THE SAME LABORATORY (same
    apparatus and method; near 5 MPa its depressions repeat [S23] within
    0.02 K, and 5 of its 12 pure-water points are [S23]'s).  What is new:
    5 % and 44.7 %, 3 MPa and 10–13 MPa.  Obtained |ΔT_h error| ≤ 0.07 K at
    5 %, 0.19 at 10 %, 0.39 at 20 %, 0.78 at 30 %, 0.98 at 40 %, 0.94 at
    44.7 %, 1.39 K at 50 %.  Sign: at 3 MPa the model under-predicts the
    depression (safe side, the −1.39 K), above ~7 MPa it over-predicts it
    (unsafe side, +0.3 to +1.02 K at 9.7–12.8 MPa, 20–50 %).  Measured
    depressions are flat in pressure (31.53–31.64 K at 50 %); the model's
    grow (30.2 K at 3 MPa, 32.7 at 13, 34.4 at 30): above 13 MPa
    unvalidated and likely more on the unsafe side.  Tolerance per concentration:
    a regression guard set after the results (≈ 1.2× the obtained), not a
    claimed accuracy; [S24] gives no uncertainty for ΔT_h.
  - Methanol and MEG against [S23], the data they were fitted to: the
    fitted concentrations are consistency only; the MEG 20 % and 40 % points
    were left out of the fit (obtained ≤ 0.13 and ≤ 0.41 K): a check within
    the same dataset, not independent validation.
  - Ice points of methanol solutions [S24, Table 10] from the same activity:
    within 0.07 K up to 10 %, too high by 1.0–8.9 K at 20–50 % (the activity
    is fitted at 250–285 K and extrapolated to 218–258 K); reported, and
    the ice point is not used as a limit for that reason.
  - Pure water, methane [S24, Table 1, 2.9–14.2 MPa]: the vdW-P model is
    0.60–0.91 K colder than measured (in test_hydrate.py).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT / "src", ROOT / "tests"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from angelica.properties.hydrate import (  # noqa: E402
    hydrate_formation_pressure,
    hydrate_formation_temperature,
    inhibitor_depression_k,
)
from angelica.properties.inhibitor import check_inhibitor, freezing_point_k, water_activity  # noqa: E402

CH4, ONE = ["methane"], [1.0]

# [S24] Table 1, methanol only: mass fraction → (T K, P MPa, ΔT_h K)
S24 = {
    0.0500: [(272.83, 3.10, 2.12), (278.10, 5.27, 2.12), (281.04, 7.19, 2.14), (283.76, 9.79, 2.15),
             (286.01, 12.77, 2.14)],
    0.1001: [(270.26, 3.01, 4.41), (275.40, 5.05, 4.42), (281.50, 9.82, 4.46), (283.68, 12.79, 4.48)],
    0.2000: [(264.73, 2.99, 9.86), (269.83, 4.98, 9.84), (276.55, 10.37, 9.88), (278.19, 12.71, 9.91)],
    0.3001: [(257.88, 2.95, 16.58), (263.12, 4.90, 16.41), (269.53, 9.78, 16.38), (271.69, 12.75, 16.45)],
    0.3999: [(251.13, 3.09, 23.77), (255.69, 4.76, 23.54), (262.32, 9.73, 23.55), (264.48, 12.76, 23.66)],
    0.4468: [(252.75, 5.09, 27.14), (256.15, 7.29, 27.15), (258.89, 9.98, 27.18), (260.95, 12.94, 27.31)],
    0.4999: [(243.02, 2.98, 31.53), (248.42, 5.03, 31.38), (254.50, 9.76, 31.43), (256.49, 12.76, 31.64)],
}
GUARD_S24 = {0.0500: 0.1, 0.1001: 0.25, 0.2000: 0.5, 0.3001: 0.95, 0.3999: 1.2, 0.4468: 1.15, 0.4999: 1.7}

# [S23] Table 1: mass fraction → (P MPa, ΔT_h K)
S23_MEG = {
    0.1: [(4.99, 2.20), (5.99, 2.24), (6.86, 2.18), (7.43, 2.21), (8.18, 2.20)],
    0.2: [(5.02, 5.43), (6.04, 5.44), (6.87, 5.45), (7.45, 5.45), (8.21, 5.46)],
    0.3: [(4.97, 9.32), (5.97, 9.32), (6.88, 9.35), (7.46, 9.35), (8.15, 9.34)],
    0.4: [(4.94, 14.92), (5.96, 14.91), (6.97, 14.93), (7.43, 14.91), (8.18, 14.93)],
    0.5: [(4.99, 22.20), (5.98, 22.23), (6.89, 22.27), (7.36, 22.30), (8.27, 22.32)],
}
S23_MEOH = {
    0.1: [(4.97, 4.41), (5.99, 4.44), (6.90, 4.43), (7.45, 4.44), (8.19, 4.44)],
    0.2: [(4.99, 9.84), (5.99, 9.83), (6.91, 9.85), (7.43, 9.83), (8.15, 9.86)],
    0.3: [(4.97, 16.41), (5.94, 16.38), (7.01, 16.35), (7.50, 16.37), (8.14, 16.35)],
    0.4: [(4.94, 23.54), (5.84, 23.51), (5.87, 23.56), (6.91, 23.47), (6.93, 23.45), (7.49, 23.53), (8.35, 23.54)],
    0.5: [(4.99, 31.36), (5.93, 31.37), (5.94, 31.36), (6.90, 31.37), (7.52, 31.31), (7.53, 31.34), (7.99, 31.42)],
}


def depression(name: str, w: float, p_mpa: float) -> float:
    pure = hydrate_formation_temperature(CH4, ONE, p_mpa * 1e6).temperature_c
    return pure - hydrate_formation_temperature(CH4, ONE, p_mpa * 1e6, inhibitor=(name, w)).temperature_c


@pytest.mark.parametrize("w", sorted(S24))
def test_methanol_against_independent_samples(w):
    errors = [depression("methanol", w, p) - dt for _t, p, dt in S24[w]]
    assert max(abs(e) for e in errors) < GUARD_S24[w], errors


@pytest.mark.parametrize("w", sorted(S23_MEOH))
def test_methanol_reproduces_its_fitting_data(w):
    """Consistency: obtained |error| ≤ 0.41 K (the fit is at 6.5 MPa; the
    error grows towards 5 and 8 MPa).  Guard 0.5 K."""
    assert all(abs(depression("methanol", w, p) - dt) < 0.5 for p, dt in S23_MEOH[w])


@pytest.mark.parametrize("w, guard", [(0.1, 0.25), (0.2, 0.2), (0.3, 0.15), (0.4, 0.5), (0.5, 0.25)])
def test_meg_against_its_dataset(w, guard):
    """Fitted at 10, 30, 50 %; 20 and 40 % left out (obtained max 0.13 and
    0.41 K).  Guards ≈ 1.2× the obtained (regression guards)."""
    assert all(abs(depression("MEG", w, p) - dt) < guard for p, dt in S23_MEG[w])


def test_ice_points_of_methanol_solutions():
    """[S24] Table 10.  Close at low concentration; too high at high
    concentration (recorded, see the module docstring)."""
    measured = {0.05: 270.045, 0.1001: 266.575, 0.2: 257.687, 0.3001: 246.20, 0.3999: 233.12, 0.4999: 217.605}
    got = {w: freezing_point_k("methanol", w) - t for w, t in measured.items()}
    assert abs(got[0.05]) < 0.1 and abs(got[0.1001]) < 0.1
    assert got[0.4999] == pytest.approx(8.87, abs=0.05)                     # the recorded discrepancy


def test_inhibitor_inputs():
    assert check_inhibitor(None) is None and check_inhibitor(("methanol", 0.0)) is None
    assert check_inhibitor(("MeOH", 0.2)) == ("methanol", 0.2)
    assert check_inhibitor(("ethylene glycol", 0.3)) == ("MEG", 0.3)
    for bad in (("methanol", 0.6), ("MEG", -0.1), ("ethanol", 0.2), ("methanol", float("nan"))):
        with pytest.raises(ValueError):
            check_inhibitor(bad)
    assert water_activity("methanol", 0.0, 270.0) == 1.0
    assert water_activity("MEG", 0.5, 260.0) < water_activity("MEG", 0.3, 260.0) < 1.0


def test_below_zero_with_an_inhibitor():
    """With an inhibitor the liquid-water branch extends below 0 °C (to
    −60 °C); without it, below 0 °C stays an error."""
    pt = hydrate_formation_pressure(CH4, ONE, -10.0, inhibitor=("methanol", 0.4))
    assert pt is not None and pt.pressure_pa > 0.0
    with pytest.raises(ValueError):
        hydrate_formation_pressure(CH4, ONE, -10.0)
    with pytest.raises(ValueError):
        hydrate_formation_pressure(CH4, ONE, -70.0, inhibitor=("methanol", 0.4))
    t = hydrate_formation_temperature(CH4, ONE, pt.pressure_pa, inhibitor=("methanol", 0.4)).temperature_c
    assert t == pytest.approx(-10.0, abs=1e-3)                    # the two solves are inverse


def test_black_oil_depression_is_the_smaller_of_two_gases():
    """Black-oil (gas-gravity chart) shift: min of the methane (sI) and
    C1/C2/C3 (sII) depressions.  Obtained at 70 bar, 30 % methanol: 16.36
    and 14.91 K."""
    names, zs = ["methane", "ethane", "propane"], [0.90, 0.07, 0.03]
    inh = ("methanol", 0.3)
    sii = (hydrate_formation_temperature(names, zs, 70e5).temperature_c
           - hydrate_formation_temperature(names, zs, 70e5, inhibitor=inh).temperature_c)
    si = depression("methanol", 0.3, 7.0)
    assert inhibitor_depression_k(70e5, inh) == pytest.approx(min(si, sii), abs=1e-9)
    assert si == pytest.approx(16.36, abs=0.01) and sii == pytest.approx(14.91, abs=0.01)
    # At 1 MPa the pipe pressure is used (v1.37.0 review: holding the 3 MPa
    # value over-shifted by 0.14 K here, the unsafe side).  Obtained 14.24 K.
    low = (hydrate_formation_temperature(names, zs, 1e6).temperature_c
           - hydrate_formation_temperature(names, zs, 1e6, inhibitor=inh).temperature_c)
    assert inhibitor_depression_k(1e6, inh) == pytest.approx(low, abs=1e-9)
    assert low == pytest.approx(14.24, abs=0.01)
    assert inhibitor_depression_k(70e5, None) == 0.0


def _pipe_case(fluid, boundaries):
    from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
    from angelica.core.components import Pipe
    return NetworkCase(name="p", fluid_model=fluid, pressure_inlets=(PressureBoundary(1, boundaries[0]),),
                       pressure_outlets=(PressureBoundary(2, boundaries[1]),),
                       thermal_inlets=(ThermalBoundary(1, 20.0),),
                       components=(Pipe(component_id="p", start_node=1, end_node=2, diameter_m=0.2,
                                        length_m=3000.0, absolute_roughness_m=46e-6,
                                        heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=4.0),))


def test_compositional_solver_uses_the_inhibitor():
    """The per-pipe hydrate check equals the direct call with the inhibitor
    (tolerance 1e-9 K: the same function at the same state)."""
    from angelica.properties.compositional_fluid import CompositionalFluid
    from angelica.solvers import SteadyCompositionalSolver
    names, zs = ["methane", "ethane", "propane", "water"], [0.88, 0.07, 0.03, 0.02]
    base = SteadyCompositionalSolver().solve(_pipe_case(CompositionalFluid(names, zs), (80e5, 70e5)))
    res = SteadyCompositionalSolver().solve(
        _pipe_case(CompositionalFluid(names, zs, hydrate_inhibitor=("MEG", 0.3)), (80e5, 70e5)))
    f0, f1 = base.component_flows[0], res.component_flows[0]
    p = res.node_pressures_pa
    direct = max(hydrate_formation_temperature(names, zs, p[n], inhibitor=("MEG", 0.3)).temperature_c
                 - res.node_temperatures_c[n] for n in (1, 2))
    assert f1.hydrate_subcooling_k == pytest.approx(direct, abs=1e-9)
    assert f1.hydrate_temperature_c < f0.hydrate_temperature_c - 8.0          # ~9 K lower with 30 % MEG


def test_black_oil_solver_uses_the_inhibitor():
    """Each end's Katz hydrate temperature is shifted by the depression at
    its pressure; the report keeps the worse end (tolerance 1e-9 K: the
    same functions at the same state)."""
    from angelica.properties.black_oil import BlackOilFluid
    from angelica.properties.hydrate import katz_hydrate_temperature
    from angelica.solvers import SteadyBlackOilSolver
    kw = {"api_gravity": 35.0, "gas_gravity": 0.7, "gor_sc_m3_per_m3": 150.0, "wor_sc_m3_per_m3": 0.2}
    inh = ("methanol", 0.25)
    base = SteadyBlackOilSolver().solve(_pipe_case(BlackOilFluid(**kw), (60e5, 50e5))).component_flows[0]
    res = SteadyBlackOilSolver().solve(_pipe_case(BlackOilFluid(**kw, hydrate_inhibitor=inh), (60e5, 50e5)))
    flow = res.component_flows[0]
    expected = max(katz_hydrate_temperature(0.7, res.node_pressures_pa[n]).temperature_c
                   - inhibitor_depression_k(res.node_pressures_pa[n], inh) - res.node_temperatures_c[n]
                   for n in (1, 2))
    assert flow.hydrate_subcooling_k == pytest.approx(expected, abs=1e-9)
    assert flow.hydrate_temperature_c < base.hydrate_temperature_c - 10.0


def test_gui_material_carries_the_inhibitor():
    from angelica.gui.io import build_network_case_from_scene, load_scene_from_file, material_inhibitor
    assert material_inhibitor({}) is None and material_inhibitor({"hydrate_inhibitor": "none"}) is None
    assert material_inhibitor({"hydrate_inhibitor": "MEG", "hydrate_inhibitor_mass_percent": "30"}) == ("MEG", 0.3)
    for bad in ({"hydrate_inhibitor": "MEG", "hydrate_inhibitor_mass_percent": "x"},
                {"hydrate_inhibitor": "MEG", "hydrate_inhibitor_mass_percent": "70"}):
        with pytest.raises(ValueError, match="Material"):
            material_inhibitor(bad)
    scene = load_scene_from_file(ROOT / "src/angelica/tutorials/steady_black_oil/07_hot_oil_pipeline/"
                                 "07_hot_oil_pipeline.gui.json")
    scene.update_material({**scene.material, "hydrate_inhibitor": "methanol",
                           "hydrate_inhibitor_mass_percent": "20"})
    assert build_network_case_from_scene(scene).fluid_model.hydrate_inhibitor == ("methanol", 0.2)


def test_hydrate_partition_with_an_inhibitor():
    """The inhibitor concentrates as water goes into hydrate; the leftover
    gas then sits on its hydrate curve with the final concentration.
    Obtained (C1/C2/C3 + 50 % water, 5 °C, 70 bar, 20 % methanol): 0.1069
    mol of water in hydrate, methanol 20 % → 24.1 % (without inhibitor all
    0.49995 mol of water goes into hydrate)."""
    from angelica.properties.hydrate_flash import hydrate_partition
    names = ["methane", "ethane", "propane", "water"]
    zs = [x * 0.5 for x in (0.90, 0.07, 0.03)] + [0.5]
    r = hydrate_partition(names, zs, 5.0, 70e5, inhibitor=("methanol", 0.2))
    assert r.limited_by == "equilibrium"
    assert r.hydrate_water_mol == pytest.approx(0.1069, abs=1e-4)
    assert r.inhibitor_mass_fraction == pytest.approx(0.2413, abs=1e-4)
    gas = [r.fluid_mol[n] for n in names[:3]]
    t = hydrate_formation_temperature(names[:3], gas, 70e5, inhibitor=("methanol", r.inhibitor_mass_fraction))
    assert t.temperature_c == pytest.approx(5.0, abs=1e-5)
    assert r.fluid_mol["water"] + r.hydrate_water_mol == pytest.approx(0.5, abs=1e-12)


@pytest.mark.parametrize("inh", [("methanol", 1e-6), ("methanol", 0.3), ("MEG", 0.1)])
def test_hydrate_partition_without_free_water_ignores_the_inhibitor(inh):
    """An undersaturated gas (50 ppm water, 0 °C, 100 bar) has no free water
    for the inhibitor to act in: the same hydrate as without it (v1.37.0
    review: it returned none, the unsafe side).  Obtained 4.3755e-6 mol."""
    from angelica.properties.hydrate_flash import hydrate_partition
    names = ["methane", "ethane", "propane", "water"]
    z_w = 5e-5
    zs = [x * (1.0 - z_w) for x in (0.90, 0.07, 0.03)] + [z_w]
    plain = hydrate_partition(names, zs, 0.0, 100e5)
    r = hydrate_partition(names, zs, 0.0, 100e5, inhibitor=inh)
    assert r.limited_by == plain.limited_by == "water"
    assert r.hydrate_water_mol == pytest.approx(plain.hydrate_water_mol, rel=1e-12)
    assert r.hydrate_water_mol == pytest.approx(4.3755e-6, rel=1e-4)
    assert r.inhibitor_mass_fraction is None


def test_absolute_error_against_the_2024_data_at_high_pressure():
    """What the hydrate margin sees is the absolute temperature: the
    inhibited model against [S24] at ≥ 9.7 MPa is too cold by up to 1.89 K
    (50 %, 12.76 MPa) — the pure-water offset (−0.9 K) and the depression
    over-prediction (+1.0 K) add up, on the unsafe side.  Recorded, not a
    tolerance claim: the bounds bracket the obtained value."""
    worst = 0.0
    for w, pts in S24.items():
        for t, p, _dt in pts:
            if p < 9.7:
                continue
            got = hydrate_formation_temperature(CH4, ONE, p * 1e6, inhibitor=("methanol", w)).temperature_c + 273.15
            worst = min(worst, got - t)
    assert -2.2 < worst < -1.5
