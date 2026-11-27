"""
Water content of gas, water dew point and the hydrate check of an
undersaturated gas (roadmap item 2b.1, v1.35.0).

Validation
----------
Water saturation pressure — READ IN SOURCE: NIST Chemistry WebBook, water,
saturation properties (IAPWS-95), 1–91 °C (retrieved 2026-11-27).
Angelica's Wagner form: deviations −7.2e-5 … +4.9e-5.  Tolerance 1e-4: a
regression guard set after the results (no published accuracy of this
Wagner form was read).  With Raoult's law the water dew point is the
inverse of this curve: from NIST's pressures it returns NIST's
temperatures within 0.002 K (the 7e-5 above over d ln P_sat/dT ≈ 0.06 /K).

NOT VALIDATED: the water content of gas at pressure (the ideal-vapour model
neglects the non-ideal enhancement; no measured water-content data were
read), and hydrate–vapour equilibria (no measured data read).  The
hydrate check is checked by consistency: at saturation it is the
water-saturated curve exactly; below saturation it lies between the water
dew point and that curve, and a dry enough gas forms no hydrate above 0 °C.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from angelica.properties.free_water import water_psat_pa  # noqa: E402
from angelica.properties.hydrate import hydrate_formation_temperature  # noqa: E402
from angelica.properties.water_content import (  # noqa: E402
    mol_per_sm3,
    saturated_water_mole_fraction,
    water_dew_point_c,
    water_mole_fraction,
)

NIST_PSAT = [(1.0, 657.086), (11.0, 1312.97), (21.0, 2488.22), (31.0, 4496.93), (41.0, 7787.84),
             (51.0, 12978.1), (61.0, 20887.9), (71.0, 32575.2), (81.0, 49367.3), (91.0, 72889.8)]
NAMES, Z = ["methane", "ethane", "propane"], [0.90, 0.07, 0.03]
P = 70e5


@pytest.mark.parametrize("t_c, p_pa", NIST_PSAT)
def test_saturation_pressure_against_nist(t_c, p_pa):
    assert abs(water_psat_pa(t_c + 273.15) / p_pa - 1.0) < 1e-4


@pytest.mark.parametrize("t_c, p_pa", NIST_PSAT)
def test_dew_point_inverts_the_saturation_curve(t_c, p_pa):
    y = p_pa / 50e5                                 # a gas at 50 bar holding NIST's P_sat as partial pressure
    assert water_dew_point_c(y, 50e5) == pytest.approx(t_c, abs=0.002)


def test_water_content_inputs():
    assert water_mole_fraction(ppm_mol=250.0) == pytest.approx(250e-6)
    y = water_mole_fraction(mg_per_sm3=100.0)       # 100 mg/Sm³ = 5.551e-3 mol water per 42.21 mol gas
    assert y == pytest.approx(100e-6 / 0.018015 / mol_per_sm3())
    assert y * 1e6 == pytest.approx(131.5, abs=0.1)
    assert water_mole_fraction(saturation=0.5, temperature_c=20.0, pressure_pa=P) == pytest.approx(
        0.5 * saturated_water_mole_fraction(20.0, P))
    assert water_dew_point_c(0.0, P) is None
    for bad in ({}, {"ppm_mol": 1.0, "mg_per_sm3": 1.0}, {"ppm_mol": -1.0}, {"saturation": 1.2,
                "temperature_c": 20.0, "pressure_pa": P}, {"saturation": 0.5}):
        with pytest.raises(ValueError):
            water_mole_fraction(**bad)


def test_saturated_gas_reproduces_the_water_saturated_curve():
    base = hydrate_formation_temperature(NAMES, Z, P).temperature_c          # 17.03 °C
    for dew in (base, base + 1.0, 40.0):
        y = saturated_water_mole_fraction(dew, P)
        t = hydrate_formation_temperature(NAMES, Z, P, water_mole_fraction_gas=y).temperature_c
        assert t == pytest.approx(base, abs=1e-6)


@pytest.mark.parametrize("dew", [15.0, 10.0, 5.0, 0.0])
def test_undersaturated_gas_forms_hydrate_between_its_dew_point_and_the_saturated_curve(dew):
    """Obtained (70 bar): dew 15 → 15.43 °C, 5 → 7.22, 0 → 3.06 (saturated
    curve 17.03)."""
    base = hydrate_formation_temperature(NAMES, Z, P).temperature_c
    y = saturated_water_mole_fraction(dew, P)
    t = hydrate_formation_temperature(NAMES, Z, P, water_mole_fraction_gas=y).temperature_c
    assert dew < t < base


def test_dry_gas_forms_no_hydrate_above_zero():
    y = saturated_water_mole_fraction(-10.0, P)          # 41 ppm
    assert hydrate_formation_temperature(NAMES, Z, P, water_mole_fraction_gas=y) is None


def test_compositional_pipe_uses_the_stream_water_content():
    """An all-vapour gas line with 60 ppm of water: its hydrate margin is the
    hydrate–vapour one, not the saturated.  Obtained: 0.420 °C at the 80 bar
    end (subcooling −11.58 K), none above 0 °C at the 70 bar end.  The
    report swallows exceptions, so the value is compared with a direct call
    (tolerance 1e-9 K: the same function at the same state)."""
    from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
    from angelica.core.components import Pipe
    from angelica.properties.compositional_fluid import CompositionalFluid
    from angelica.solvers import SteadyCompositionalSolver

    def run(names, zs, full=False):
        case = NetworkCase(name="g", fluid_model=CompositionalFluid(names, zs),
                           pressure_inlets=(PressureBoundary(1, 80e5),), pressure_outlets=(PressureBoundary(2, 70e5),),
                           thermal_inlets=(ThermalBoundary(1, 12.0),),
                           components=(Pipe(component_id="p", start_node=1, end_node=2, diameter_m=0.2,
                                            length_m=3000.0, absolute_roughness_m=46e-6,
                                            heat_transfer_coefficient_w_per_m2k=2.0, ambient_temperature_c=8.0),))
        res = SteadyCompositionalSolver().solve(case)
        return res if full else res.component_flows[0]
    wet = run(NAMES, Z)                                                    # no water component: in excess
    y = 60e-6
    names, zs = NAMES + ["water"], [z * (1 - y) for z in Z] + [y]
    res = run(names, zs, full=True)
    dry = res.component_flows[0]
    assert wet.hydrate_temperature_c is not None and wet.hydrate_subcooling_k > 0.0     # inside the curve
    direct = hydrate_formation_temperature(names, zs, 80e5, water_mole_fraction_gas=y)
    assert direct is not None and direct.temperature_c == pytest.approx(0.420, abs=1e-3)
    assert dry.hydrate_temperature_c == pytest.approx(direct.temperature_c, abs=1e-9)
    assert dry.hydrate_subcooling_k == pytest.approx(direct.temperature_c - 12.0, abs=1e-9)
    assert hydrate_formation_temperature(names, zs, 70e5, water_mole_fraction_gas=y) is None


def test_invalid_water_inputs():
    """Inputs outside the model's range raise instead of returning a clamped
    or extrapolated number (v1.35.0 review)."""
    for bad in (lambda: water_dew_point_c(1e-6, 1e5),                 # dew point below −50 °C
                lambda: water_dew_point_c(float("nan"), P),
                lambda: water_dew_point_c(1.0, 300e5),                # above the critical point
                lambda: water_mole_fraction(saturation=0.5, temperature_c=20.0, pressure_pa=0.0),
                lambda: water_mole_fraction(ppm_mol=2e6),
                lambda: water_mole_fraction(mg_per_sm3=1e9),
                lambda: hydrate_formation_temperature(NAMES, Z, P, water_mole_fraction_gas=-1e-6),
                lambda: hydrate_formation_temperature(NAMES, Z, P, water_mole_fraction_gas=float("nan"))):
        with pytest.raises(ValueError):
            bad()
    assert hydrate_formation_temperature(NAMES, Z, P, water_mole_fraction_gas=0.0) is None   # no water


def test_curve_at_water_content_follows_the_solver_rule():
    """GUI curve of a feed with water (v1.35.0 review: it used z_w regardless
    of the phases).  Obtained for 60 ppm: 77.26 bar at 0 °C and 116.45 at
    5 °C (hydrate–vapour, gas all vapour); at 10 °C the hydrate–vapour point
    (173 bar) is dense phase that the flash labels liquid, so the saturated
    point 25.51 bar.  A condensate feed keeps the saturated curve."""
    from angelica.properties.hydrate import hydrate_curve, hydrate_curve_at_water_content
    y = 60e-6
    names, zs = NAMES + ["water"], [z * (1 - y) for z in Z] + [y]
    got = hydrate_curve_at_water_content(names, zs, 0.0, 20.0, 5)
    sat = hydrate_curve(names, zs, 0.0, 20.0, 5)
    assert [p.pressure_pa / 1e5 for p in got[:3]] == pytest.approx([77.26, 116.45, 25.51], abs=0.01)
    assert got[2:] == sat[2:]
    cond_names = ["methane", "propane", "n-heptane", "water"]
    cond = [0.8 * (1 - 1e-4), 0.1 * (1 - 1e-4), 0.1 * (1 - 1e-4), 1e-4]
    assert hydrate_curve_at_water_content(cond_names, cond, 0.0, 20.0, 5) == hydrate_curve(cond_names, cond,
                                                                                         0.0, 20.0, 5)
    assert hydrate_curve_at_water_content(NAMES, Z, 0.0, 20.0, 5) == hydrate_curve(NAMES, Z, 0.0, 20.0, 5)
