"""
Inner film coefficient of a pipe from the flow (roadmap item 6.1).

    h_i = Nu·k/D,   Re = ρ·|v|·D/μ = 4·|ṁ|/(π·D·μ),   Pr = μ·c_p/k

  - Re ≤ 2300: fully developed laminar flow at uniform wall temperature,
    Nu = 3.66 (same secondary source; 4.36 at uniform heat flux — the wall
    temperature case is taken; entrance effects are not modelled);
  - Re ≥ 3000: Gnielinski (1975, Forsch. Ing.-Wes. 41, 8–16), with the
    friction factor f = (0.79·ln Re − 1.64)⁻² (Filonenko),

        Nu = (f/8)(Re − 1000)Pr / (1 + 12.7·(f/8)^½·(Pr^⅔ − 1)),

    stated for 0.5 ≤ Pr ≤ 2000 and 3000 ≤ Re ≤ 5·10⁶, smooth tubes (read in
    a secondary source: Wikipedia, "Nusselt number", citing Incropera &
    DeWitt p. 490; read 2026-11-27).  Outside, it is used and
    ``InnerFilm.in_range`` is False;
  - in between: linear in Re between the two ends (no correlation covers
    the transition; this only avoids a jump).

Mixture properties are those the fluid model gives the pipe (homogeneous
for two-phase flow: the film of a two-phase flow is not modelled).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

NU_LAMINAR = 3.66
RE_LAMINAR, RE_TURBULENT = 2300.0, 3000.0


@dataclass(frozen=True)
class InnerFilm:
    h_w_per_m2k: float
    nusselt: float
    reynolds: float
    prandtl: float
    in_range: bool


def gnielinski_nusselt(reynolds: float, prandtl: float) -> float:
    f = (0.79 * math.log(reynolds) - 1.64) ** -2
    return (f / 8.0) * (reynolds - 1000.0) * prandtl / (1.0 + 12.7 * math.sqrt(f / 8.0) * (prandtl ** (2.0 / 3.0) - 1.0))


def nusselt(reynolds: float, prandtl: float) -> float:
    """Nusselt number of the inner film (see the module docstring)."""
    if reynolds <= RE_LAMINAR:
        return NU_LAMINAR
    if reynolds >= RE_TURBULENT:
        return gnielinski_nusselt(reynolds, prandtl)
    w = (reynolds - RE_LAMINAR) / (RE_TURBULENT - RE_LAMINAR)
    return (1.0 - w) * NU_LAMINAR + w * gnielinski_nusselt(RE_TURBULENT, prandtl)


def inner_film(mass_flow_kg_per_s: float, diameter_m: float, viscosity_pa_s: float,
               specific_heat_j_per_kg_k: float, conductivity_w_per_mk: float,
               length_m: float | None = None) -> InnerFilm:
    """Inner film coefficient of the flow ``mass_flow_kg_per_s`` in a bore of
    ``diameter_m`` with the given fluid properties.  ``in_range`` is False
    outside Gnielinski's stated range, in the 2300–3000 interpolation, and
    for laminar flow whose thermal entry length 0.05·Re·Pr·D exceeds the
    pipe ``length_m`` (developing flow: the real film is higher than 3.66)."""
    mu, cp, k = viscosity_pa_s, specific_heat_j_per_kg_k, conductivity_w_per_mk
    values = (mass_flow_kg_per_s, diameter_m, mu, cp, k)
    if not all(math.isfinite(v) for v in values) or not (mu > 0.0 and cp > 0.0 and k > 0.0 and diameter_m > 0.0):
        raise ValueError("inner film needs a finite flow and positive viscosity, specific heat, "
                         "conductivity and diameter")
    re = 4.0 * abs(mass_flow_kg_per_s) / (math.pi * diameter_m * mu)
    pr = mu * cp / k
    nu = nusselt(re, pr)
    if re <= RE_LAMINAR:
        ok = length_m is None or 0.05 * re * pr * diameter_m <= length_m
    elif re < RE_TURBULENT:
        ok = False
    else:
        ok = 0.5 <= pr <= 2000.0 and re <= 5e6
    return InnerFilm(nu * k / diameter_m, nu, re, pr, ok)
