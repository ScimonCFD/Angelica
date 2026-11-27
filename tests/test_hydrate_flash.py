"""
Hydrate as a phase: amount of hydrate formed (roadmap item 2b.2, v1.36.0).

Validation
----------
Methane hydrate composition — COMPARED, NOT VALIDATED (a recorded
discrepancy).  Read in a secondary source (a search-engine
excerpt of the USGS abstract of Sum, Burruss & Sloan 1997, J. Phys. Chem. B
101, 7371; the abstract page itself refused the fetch, HTTP 403): Raman
hydration number 6.04 ± 0.03.  The conditions of the measurement were not
read.  The model gives 5.895 at 2 °C and 70 bar (−2.4 %, outside the stated
±0.03); the test pins the obtained value and the stoichiometric bound of a
full sI lattice, 46 waters / 8 cages = 5.75 (structure data of the module).

Consistency (the Done-when of the roadmap item):
  - with excess water the incipient boundary is the vdW-P curve: just above
    the formation temperature no hydrate forms, just below a little does,
    and the amount grows continuously from zero;
  - mass balance of the water and of each former closes (tolerance 1e-12:
    sums of the same numbers);
  - at equilibrium the remaining gas is exactly at its own hydrate
    condition (the formation temperature of the leftover gas at P is T,
    to 1e-5 K: the formation solve stops at 1e-4 in T, obtained ≤ 1.5e-6);
  - when the free water runs out, the gas's water content puts it on the
    hydrate–vapour curve of roadmap 2b.1 at T (same tolerance).  Both use
    the same water-activity expression, so this checks that the solve
    converged to that condition, not the physics of it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from angelica.properties.hydrate import hydrate_formation_temperature  # noqa: E402
from angelica.properties.hydrate_flash import hydrate_partition  # noqa: E402

HC = ["methane", "ethane", "propane"]
NAMES = HC + ["water"]
P = 70e5


def feed(z_w: float) -> list[float]:
    return [x * (1.0 - z_w) for x in (0.90, 0.07, 0.03)] + [z_w]


def balance(r, names, zs):
    s = sum(zs)
    for n, z in zip(names, zs):
        held = r.hydrate_water_mol if n == "water" else r.hydrate_guests_mol.get(n, 0.0)
        assert r.fluid_mol[n] + held == pytest.approx(z / s, abs=1e-12)


def leftover_gas(r):
    return [r.fluid_mol[n] for n in HC]


def test_incipient_boundary_is_the_vdwp_curve():
    t_eq = hydrate_formation_temperature(HC, [0.9, 0.07, 0.03], P).temperature_c       # 17.03 °C
    assert hydrate_partition(NAMES, feed(0.5), t_eq + 0.05, P).limited_by == "none"
    amounts = [hydrate_partition(NAMES, feed(0.5), t_eq - d, P).hydrate_water_mol for d in (0.01, 0.05, 0.5)]
    assert 0.0 < amounts[0] < amounts[1] < amounts[2]
    # Linear in the subcooling near the curve, so it extrapolates to zero at
    # t_eq.  Obtained 0.1984 and 0.1976 mol/K (0.4 %); tolerance 2 %: the
    # curvature of a smooth function over 0.05 K.
    assert amounts[0] / 0.01 == pytest.approx(amounts[1] / 0.05, rel=0.02)


def test_equilibrium_with_free_water_left():
    """95 % water at 12 °C: the gas gives up propane (sII) until it sits on
    its own hydrate curve at 12 °C.  Obtained: 0.0888 mol water in hydrate,
    hydration number 6.156, propane in the gas 3.0 % → 0.27 %."""
    zs = feed(0.95)
    r = hydrate_partition(NAMES, zs, 12.0, P)
    assert r.limited_by == "equilibrium" and r.structure == "sII" and r.free_water_mol > 0.0
    balance(r, NAMES, zs)
    gas = leftover_gas(r)
    assert gas[2] / sum(gas) < 0.1 * 0.03                      # propane goes into the hydrate first
    assert hydrate_formation_temperature(HC, gas, P).temperature_c == pytest.approx(12.0, abs=1e-5)
    assert r.hydrate_water_mol == pytest.approx(0.08881, abs=1e-5)


def test_water_limited_gas_ends_on_the_hydrate_vapour_curve():
    """1 % water at 5 °C: all free water goes into hydrate; the gas keeps
    1.02e-4 mol of water vapour, which puts it on its hydrate–vapour curve."""
    zs = feed(0.01)
    r = hydrate_partition(NAMES, zs, 5.0, P)
    assert r.limited_by == "water" and r.free_water_mol == 0.0
    assert r.vapour_water_mol == pytest.approx(1.0210e-4, rel=1e-4)        # obtained 1.02096e-4
    balance(r, NAMES, zs)
    gas = leftover_gas(r)
    y_w = r.vapour_water_mol / (sum(gas) + r.vapour_water_mol)
    t = hydrate_formation_temperature(HC, gas, P, water_mole_fraction_gas=y_w).temperature_c
    assert t == pytest.approx(5.0, abs=1e-5)


def test_pure_methane_is_limited_by_the_gas():
    """Methane's fugacity does not fall as it is consumed: hydrate grows until
    the methane is used up.  Hydration number 5.895 (sI): see the module
    docstring for the comparison with Sum et al. (1997)."""
    zs = [0.1, 0.9]
    r = hydrate_partition(["methane", "water"], zs, 2.0, P)
    assert r.limited_by == "formers" and r.structure == "sI"
    assert r.other_structure_stable                  # sII methane is also below its curve at 70 bar
    balance(r, ["methane", "water"], zs)
    assert r.fluid_mol["methane"] < 1e-9
    assert r.hydration_number == pytest.approx(5.8949, abs=1e-4)
    assert r.hydration_number > 46.0 / 8.0
    assert r.hydrate_water_mol == pytest.approx(0.1 * r.hydration_number, rel=1e-8)


def test_a_former_running_out_flags_the_other_structure():
    """95 % water at 4 °C: sII takes the gas almost entirely; the methane-rich
    remainder would also form sI (two-hydrate equilibrium not computed:
    reported by the flag)."""
    zs = feed(0.95)
    r = hydrate_partition(NAMES, zs, 4.0, P)
    assert r.limited_by == "formers" and r.other_structure_stable
    balance(r, NAMES, zs)
    # Obtained 5.7e-6 of the gas left: the bisection stops when ξ is known to
    # 1e-10 and the last formers are nearly gone; the bound is 2× that.
    assert sum(leftover_gas(r)) < 1.2e-5 * sum(zs[:3])


def test_no_hydrate_and_invalid_inputs():
    assert hydrate_partition(NAMES, feed(0.5), 25.0, P).limited_by == "none"
    r = hydrate_partition(NAMES, feed(0.0), 5.0, P)                    # no water
    assert r.limited_by == "none" and r.hydrate_mol == 0.0
    assert hydrate_partition(["n-heptane", "water"], [0.5, 0.5], 2.0, P).limited_by == "none"
    for bad in (lambda: hydrate_partition(HC, [0.9, 0.07, 0.03], 5.0, P),          # no water component
                lambda: hydrate_partition(NAMES, feed(0.5), -1.0, P),             # ice branch
                lambda: hydrate_partition(NAMES, feed(0.5), 5.0, 0.0),
                lambda: hydrate_partition(NAMES, [0.5, 0.1, 0.1, float("nan")], 5.0, P)):
        with pytest.raises(ValueError):
            bad()


def test_oil_gas_water_feed():
    """Oil, gas and water (methane/propane/n-heptane + 50 % water) at 8 °C,
    70 bar: the hydrocarbons are two-phase (VF 0.43).  The water runs out;
    the hydrate takes methane and propane, the heptane stays in the fluid.
    Obtained: 0.499977 mol water in hydrate, 2.273e-5 left as vapour,
    G/RT = −7.7e-13 (a secant step on a bracket 1e-10 of ξ wide; tolerance
    1e-9, the bracket's width times dG/dξ ≈ 1/n_wv ≈ 4e4)."""
    from angelica.properties.compositional_fluid import _hc_phase_compositions
    names, zs = ["methane", "propane", "n-heptane", "water"], [0.30, 0.05, 0.15, 0.50]
    assert _hc_phase_compositions(tuple(names[:3]), P, 8.0, (0.6, 0.1, 0.3), "PR")[2] == pytest.approx(0.4313,
                                                                                                    abs=1e-4)
    r = hydrate_partition(names, zs, 8.0, P)
    assert r.limited_by == "water" and r.structure == "sII"
    balance(r, names, zs)
    assert r.fluid_mol["n-heptane"] == 0.15 and "n-heptane" not in r.hydrate_guests_mol
    assert r.vapour_water_mol == pytest.approx(2.273e-5, rel=1e-3)
    assert abs(r.driving_force) < 1e-9


@pytest.mark.parametrize("names, zs, p", [(["propane", "water"], [0.5, 0.5], 70e5),
                                          (["carbon dioxide", "water"], [0.5, 0.5], 70e5),
                                          (NAMES, feed(0.5), 300e5)])
def test_hydrocarbon_without_vapour_takes_all_the_water(names, zs, p):
    """Liquid propane, liquid CO2, and the C1/C2/C3 gas at 300 bar (dense,
    labelled liquid by PR): the free-water model gives the liquid no water,
    so the hydrate takes all of it and stays stable (G < 0).  Labelled
    "water" (v1.36.0 review: it said "equilibrium")."""
    r = hydrate_partition(names, zs, 2.0 if p < 1e7 else 5.0, p)
    assert r.limited_by == "water" and r.driving_force < 0.0
    assert r.hydrate_water_mol == pytest.approx(zs[-1] / sum(zs), rel=1e-11)
    balance(r, names, zs)


def test_little_gas_is_consumed_to_the_stated_accuracy():
    """95 %→99.999 % water: the bracket is narrowed by decades, so ξ is found
    to 1e-10 of itself however little gas there is (v1.36.0 review: with
    1e-8 of gas the old bisection left 0.4 % of it).  Obtained ξ/z_gas =
    6.0677 at 1e-5 and 1e-3, the leftover gas ≤ 6e-6 of the feed gas."""
    for z_g in (1e-3, 1e-5):
        r = hydrate_partition(NAMES, feed(1.0 - z_g), 5.0, P)
        assert r.limited_by == "formers"
        assert r.hydrate_water_mol / z_g == pytest.approx(6.0677, abs=1e-3)
        assert sum(leftover_gas(r)) < 1.2e-5 * z_g
