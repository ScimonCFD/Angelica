"""
Thermodynamic hydrate inhibitors in the free water (roadmap item 2b.3).

The inhibitor lowers the activity of water in the aqueous phase, and with it
the hydrate temperature: the van der Waals–Platteeuw water term uses
a_w = γ_w·x_w instead of x_w, with

    R·T·ln γ_w = x_i²·(A + B·x_i)        (x_i: inhibitor mole fraction, salt-free)

held independent of temperature (regular-solution assumption).  The
inhibitor is assumed to stay in the water (its loss to the gas and oil is
not modelled): the mass fraction is that of the free water.

A and B are fitted, through Angelica's vdW-P model, to measured methane
hydrate temperature depressions ΔT_h (the measured temperature with the
inhibitor minus that without, at the same pressure; the model's own
pure-water temperature is the reference, so its pure-water offset does not
enter): Semenov et al. (2023), Data in Brief 46, 108892, Table 1 (read in
source through Europe PMC, PMC9876827), 4.9–8.4 MPa, mean ΔT_h of each
sample at 6.5 MPa.  They are tied to this model: they are not general
activity coefficients.  Fitted with methane (sI); applied to any gas in
compositional runs without validation (at 7 MPa and 30 % methanol the
model's depression is 16.4 K for methane, 14.9 K for a C1/C2/C3 sII gas,
20.9 K for CO2).

  - methanol: fitted to 10–50 mass %; compared with new samples from the
    same laboratory, Semenov et al. (2024), Data in Brief 53, 110138,
    Table 1 (5–50 mass %, 3–13 MPa; PMC10877679): see the tests for the
    errors and their sign (over-predicted depression, the unsafe side,
    above ~7 MPa; the model's depression grows with pressure where the
    measured one is flat, so above 13 MPa it is unvalidated);
  - MEG: fitted to 10, 30 and 50 mass %; 20 and 40 % left out as a check
    (same dataset: consistency, not independent validation).

The ice point of the inhibited water (reported, not used as a limit: no ice
phase is modelled, and the hydrate module treats the water as liquid down
to −60 °C) follows from

    ln a_w(T_f) = −(ΔH_fus/R)(1/T_f − 1/T0) + (Δc_p/R)(T0/T_f − 1 + ln(T_f/T0))

with ΔH_fus = 6009.5 J/mol (the hydrate module's value) and
Δc_p = c_p,liquid − c_p,ice = 76.014 − 37.77 = 38.24 J/(mol·K), held
constant: liquid water 76.014 J/(mol·K) at 273.16 K, 0.1 MPa (NIST
Chemistry WebBook, IAPWS-95, read 2026-11-27); ice Ih 2096.71 J/(kg·K) at
the normal melting point (IAPWS R10-06(2009), Table 6, read 2026-11-27).
For methanol it is compared with the measured ice points of Semenov et al.
(2024), Table 10: within 0.1 K up to 10 %, but 1.0, 3.0, 5.3 and 8.9 K too
high at 20, 30, 40 and 50 % (the activity is fitted at the hydrate
temperatures, 250–285 K; the regular-solution extrapolation to the ice
points, 218–258 K, is poor).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

_R = 8.314462618
_T0 = 273.15
_H_FUS = 6009.5                         # J/mol, as in hydrate.py
_DCP = 76.014 - 2096.71 * 0.018015      # J/(mol·K), liquid water − ice Ih at 273.15 K
_M_WATER = 0.018015


@dataclass(frozen=True)
class Inhibitor:
    name: str
    molar_mass: float                   # kg/mol
    a: float                            # J/mol
    b: float                            # J/mol
    max_mass_fraction: float            # upper end of the data the fit used
    source: str


INHIBITORS = {
    "methanol": Inhibitor("methanol", 0.032042, -1942.01, 4859.27, 0.50,
                          "Semenov et al. 2023, Data in Brief 46, 108892, Table 1 (hydrate ΔT_h)"),
    "MEG": Inhibitor("MEG", 0.062068, -2511.60, -1338.21, 0.50,
                     "Semenov et al. 2023, Data in Brief 46, 108892, Table 1 (hydrate ΔT_h)"),
}
_ALIASES = {"meoh": "methanol", "methanol": "methanol", "meg": "MEG", "ethylene glycol": "MEG",
            "monoethylene glycol": "MEG", "ethane-1,2-diol": "MEG"}


def get_inhibitor(name: str) -> Inhibitor:
    key = _ALIASES.get(str(name).strip().lower())
    if key is None or key not in INHIBITORS:
        raise ValueError(f"unknown hydrate inhibitor {name!r} (available: {', '.join(INHIBITORS)})")
    return INHIBITORS[key]


def check_inhibitor(inhibitor) -> tuple[str, float] | None:
    """Validate ``(name, mass fraction in the free water)``; None or a zero
    fraction means no inhibitor.  Returns a hashable, rounded key."""
    if inhibitor is None:
        return None
    if isinstance(inhibitor, str) or len(inhibitor) != 2:
        raise ValueError(f"an inhibitor is (name, mass fraction), e.g. ('methanol', 0.3); got {inhibitor!r}")
    name, w = inhibitor
    inh = get_inhibitor(name)
    if not isinstance(w, (int, float)) or isinstance(w, bool):
        raise ValueError(f"inhibitor mass fraction must be a number (got {w!r})")
    w = float(w)
    if not (math.isfinite(w) and 0.0 <= w <= inh.max_mass_fraction + 1e-12):
        raise ValueError(f"{inh.name} mass fraction must be in [0, {inh.max_mass_fraction}] "
                         f"(the range of the data it is fitted to; got {w})")
    if w == 0.0:
        return None
    return (inh.name, round(w, 10))


def mole_fraction(name: str, mass_fraction: float) -> float:
    """Inhibitor mole fraction in the salt-free aqueous solution."""
    inh = get_inhibitor(name)
    n_i = mass_fraction / inh.molar_mass
    n_w = (1.0 - mass_fraction) / _M_WATER
    return n_i / (n_i + n_w)


def water_activity(name: str, mass_fraction: float, temperature_k: float) -> float:
    """a_w = γ_w·x_w of water in the inhibitor solution (no dissolved gas)."""
    if mass_fraction <= 0.0:
        return 1.0
    inh = get_inhibitor(name)
    x = mole_fraction(name, mass_fraction)
    ln_gamma = x * x * (inh.a + inh.b * x) / (_R * temperature_k)
    return (1.0 - x) * math.exp(ln_gamma)


def ice_water_activity(temperature_k: float) -> float:
    """Activity of liquid water in equilibrium with ice Ih at T (0.1 MPa)."""
    T = temperature_k
    return math.exp(-_H_FUS / _R * (1.0 / T - 1.0 / _T0) + _DCP / _R * (_T0 / T - 1.0 + math.log(T / _T0)))


def freezing_point_k(name: str, mass_fraction: float) -> float:
    """Ice point (K) of the inhibitor solution at 0.1 MPa."""
    if mass_fraction <= 0.0:
        return _T0
    lo, hi = 180.0, _T0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if water_activity(name, mass_fraction, mid) > ice_water_activity(mid):
            lo = mid                    # ice is stable: the ice point is higher
        else:
            hi = mid
        if hi - lo < 1e-9:
            break
    return 0.5 * (lo + hi)
