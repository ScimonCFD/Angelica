"""
Vertical two-phase flow correlations for well tubing: Hagedorn-Brown and Gray.

Both return a ``BeggsBrillState``-compatible result so that the ``Well``
component can integrate its tubing with any of them (``beggs_brill_gradient``
is the third option).  SI units throughout; the dimensionless groups are the
SI forms (e.g. N_LV = v_sL·(ρ_L/(g·σ))^¼, numerically equal to the field-unit
form 1.938·v_sL·(ρ_L/σ)^¼).  The acceleration term is neglected, as in the
Beggs-Brill implementation.

Hagedorn & Brown (1965), JPT 17(4), 475–484, with the modifications of
Brill & Mukherjee (1999, §4.2.1): holdup not below the no-slip value, and
Griffith bubble flow (Orkiszewski boundary) at low gas fraction.
  - Groups: N_LV, N_GV, N_D = D·(ρ_L·g/σ)^½, N_L = μ_L·(g/(ρ_L·σ³))^¼.
  - Charts as polynomial fits.  The coefficients were read in two
    independent open implementations that agree to the digits each prints:
    pyResToolbox (pyrestoolbox/nodal/nodal.py, M. Burgoyne) and rNodal
    (R/Duns-Ros.R and R/Hagedorn-Brown.R, A. Reyes), the latter crediting
    them to Guo, Lyons & Ghalambor (2007), *Petroleum Production
    Engineering: A Computer-Assisted Approach*, and to Brown's *The
    Technology of Artificial Lift Methods*.  The book itself was not read,
    so that attribution is second-hand.  Only the numbers are used; the
    implementation is Angelica's own:
        CN_L     = 10^Σ c_i·x^i,  x = log10(N_L) + 3, N_L held in [0.002, 0.4]
        H_L/ψ    = Σ y_i·u^i,     u = log10(φ1) + 6 (held at the fit maximum),
                   φ1 = N_LV·(P/P_atm)^0.1·CN_L / (N_GV^0.575·N_D)
        ψ        = Σ s_i·φ2^i,    φ2 = N_GV·N_L^0.38 / N_D^2.14 in [0.012, 0.09];
                   ψ = 1 for φ2 ≤ 0.001
  - Friction: f_D(N_Re, ε/D)·G²/(2·ρ_s·D), G = ṁ/A,
    N_Re = G·D / (μ_L^H_L·μ_g^(1−H_L)) (Hagedorn & Brown 1965, Eq. 3).
  - Bubble flow (Griffith): λ_g < L_B = 1.071 − 0.2218·v_m²/D (v in ft/s, D in
    ft; L_B ≥ 0.13), slip velocity 0.8 ft/s,
        H_L = 1 − ½·[1 + v_m/v_s − √((1 + v_m/v_s)² − 4·v_sg/v_s)],
    friction on the liquid alone: ρ_L·v_L²·f_D/(2·D), v_L = v_sL/H_L
    (whitson wiki, Hagedorn-Brown page).

Gray (1974), API 14B user's manual (lean gas-condensate wells, mist flow),
equations as given by the whitson wiki Gray page:
    H_L = 1 − (1 − e^A)/(R + 1),  R = v_sL/v_sg
    A = −2.314·[N_v·(1 + 205/N_d)]^B,  B = 0.0814·[1 − 0.0554·ln(1 + 730·R/(1+R))]
    N_v = ρ_m²·v_m²/(g·σ·(ρ_L − ρ_g)),  N_d = g·(ρ_L − ρ_g)·D²/σ
    effective roughness k_e = k0 for R ≥ 0.007, else k + R·(k0 − k)/0.007,
    k0 = 28.5·σ/(ρ_m·v_m²), k_e ≥ 2.77e-5 ft (8.44e-6 m; the floor is in feet in
    the original program — confirmed by pyResToolbox)
    gravity on ρ_s, friction f_D(ρ_m·v_m·D/μ_m, k_e/D)·ρ_m·v_m²/(2·D).
"""
from __future__ import annotations

import math

from .beggs_brill import BeggsBrillState, _darcy_friction_factor
from .gravity import GRAVITY_M_PER_S2

_G = GRAVITY_M_PER_S2
_P_ATM = 101_325.0
_FT = 0.3048

# Hagedorn-Brown chart fits (see the module docstring for their origin)
_CNL = (-2.69851, 0.1584095, -0.5509976, 0.5478492, -0.1219458)
_HL_PSI = (-0.10306578, 0.617774, -0.632946, 0.29598, -0.0401)
_PSI = (0.9116257, -4.821756, 1232.25, -22253.58, 116174.3)
_NL_MIN, _NL_MAX = 0.002, 0.4
_U_MAX = 3.6722855            # log10(φ1) + 6 at the H_L/ψ fit maximum
_PHI2_MIN, _PHI2_MAX, _PHI2_FLOOR = 0.012, 0.09, 0.001
_V_BUBBLE = 0.8 * _FT         # Griffith slip velocity (m/s)

_GRAY_ROUGH_FLOOR = 2.77e-5 * _FT


def _poly(c, x):
    return sum(ci * x**i for i, ci in enumerate(c))


def _no_slip(v_m, C_L, rho_l, rho_g, mu_l, mu_g):
    rho_ns = C_L * rho_l + (1.0 - C_L) * rho_g
    mu_ns = C_L * mu_l + (1.0 - C_L) * mu_g
    return rho_ns, mu_ns, C_L * v_m, (1.0 - C_L) * v_m


def hagedorn_brown_gradient(
    v_m: float, C_L: float, rho_l: float, rho_g: float, mu_l: float, mu_g: float,
    D: float, roughness: float, theta_rad: float, sigma_l: float = 0.030,
    pressure_pa: float = _P_ATM,
) -> BeggsBrillState:
    """Pressure gradient (Pa/m) by Hagedorn-Brown (modified); ``pressure_pa``
    enters the holdup group φ1."""
    v_m = max(abs(v_m), 1e-6)
    rho_ns, mu_ns, v_sl, v_sg = _no_slip(v_m, C_L, rho_l, rho_g, mu_l, mu_g)
    G = rho_ns * v_m
    sigma = max(sigma_l, 1e-5)

    if C_L >= 1.0 - 1e-9 or v_sg <= 1e-9:            # single-phase liquid
        HL, regime = 1.0, "liquid"
        re = G * D / max(mu_l, 1e-12)
        f = _darcy_friction_factor(re, roughness / D)
        rho_s = rho_l
        fric = f * G * G / (2.0 * rho_s * D)
    elif C_L <= 1e-9:                                  # single-phase gas
        HL, regime = 0.0, "gas"
        re = G * D / max(mu_g, 1e-12)
        f = _darcy_friction_factor(re, roughness / D)
        rho_s = rho_g
        fric = f * G * G / (2.0 * rho_s * D)
    else:
        lam_g = 1.0 - C_L
        vm_ft, d_ft = v_m / _FT, D / _FT
        L_B = max(1.071 - 0.2218 * vm_ft * vm_ft / d_ft, 0.13)
        if lam_g < L_B:                                # Griffith bubble flow
            vs = _V_BUBBLE
            a = 1.0 + v_m / vs
            disc = max(a * a - 4.0 * v_sg / vs, 0.0)
            HL = min(max(1.0 - 0.5 * (a - math.sqrt(disc)), C_L), 1.0)
            regime = "bubble"
            rho_s = HL * rho_l + (1.0 - HL) * rho_g
            v_L = v_sl / HL
            re = rho_l * v_L * D / max(mu_l, 1e-12)
            f = _darcy_friction_factor(re, roughness / D)
            fric = f * rho_l * v_L * v_L / (2.0 * D)
        else:
            grp = (rho_l / (_G * sigma)) ** 0.25
            N_LV, N_GV = v_sl * grp, v_sg * grp
            N_D = D * (rho_l * _G / sigma) ** 0.5
            N_L = mu_l * (_G / (rho_l * sigma**3)) ** 0.25
            x = math.log10(min(max(N_L, _NL_MIN), _NL_MAX)) + 3.0
            CNL = 10.0 ** _poly(_CNL, x)
            phi1 = N_LV * (max(pressure_pa, 1.0) / _P_ATM) ** 0.1 * CNL / (max(N_GV, 1e-10) ** 0.575 * N_D)
            u = min(math.log10(max(phi1, 1e-30)) + 6.0, _U_MAX)
            hl_psi = max(_poly(_HL_PSI, u), 0.0)
            phi2 = N_GV * max(N_L, 1e-12) ** 0.38 / N_D**2.14
            psi = 1.0 if phi2 <= _PHI2_FLOOR else _poly(_PSI, min(max(phi2, _PHI2_MIN), _PHI2_MAX))
            HL = min(max(psi * hl_psi, C_L), 1.0)      # not below no-slip (Brill & Mukherjee)
            regime = "hagedorn-brown"
            rho_s = HL * rho_l + (1.0 - HL) * rho_g
            mu_tp = max(mu_l, 1e-12) ** HL * max(mu_g, 1e-12) ** (1.0 - HL)
            re = G * D / mu_tp
            f = _darcy_friction_factor(re, roughness / D)
            fric = f * G * G / (2.0 * rho_s * D)

    grav = rho_s * _G * math.sin(theta_rad)
    N_Fr = v_m * v_m / (_G * D)
    return BeggsBrillState(regime=regime, holdup=HL, C_L=C_L, N_Fr=N_Fr, N_LV=0.0,
                           rho_ns=rho_ns, rho_s=rho_s, f_n=f, f_tp=f,
                           dP_dz_friction=fric, dP_dz_gravity=grav, dP_dz_total=fric + grav)


def gray_gradient(
    v_m: float, C_L: float, rho_l: float, rho_g: float, mu_l: float, mu_g: float,
    D: float, roughness: float, theta_rad: float, sigma_l: float = 0.030,
    pressure_pa: float = _P_ATM,
) -> BeggsBrillState:
    """Pressure gradient (Pa/m) by Gray (gas and gas-condensate wells)."""
    v_m = max(abs(v_m), 1e-6)
    rho_ns, mu_ns, v_sl, v_sg = _no_slip(v_m, C_L, rho_l, rho_g, mu_l, mu_g)
    sigma = max(sigma_l, 1e-5)
    if C_L <= 1e-12 or v_sl <= 0.0:
        HL, k_e = 0.0, roughness
    elif v_sg <= 1e-9:
        HL, k_e = 1.0, roughness
    else:
        R = v_sl / v_sg
        drho = max(rho_l - rho_g, 1e-6)
        N_v = rho_ns**2 * v_m**2 / (_G * sigma * drho)
        N_d = _G * drho * D * D / sigma
        B = 0.0814 * (1.0 - 0.0554 * math.log(1.0 + 730.0 * R / (1.0 + R)))
        A = -2.314 * (N_v * (1.0 + 205.0 / N_d)) ** B
        HL = 1.0 - (1.0 - math.exp(A)) / (R + 1.0)
        HL = min(max(HL, C_L), 1.0)
        k0 = 28.5 * sigma / (rho_ns * v_m * v_m)
        k_e = k0 if R >= 0.007 else roughness + R * (k0 - roughness) / 0.007
        k_e = max(k_e, _GRAY_ROUGH_FLOOR)
    rho_s = HL * rho_l + (1.0 - HL) * rho_g
    re = rho_ns * v_m * D / max(mu_ns, 1e-12)
    f = _darcy_friction_factor(re, min(k_e / D, 0.5))
    fric = f * rho_ns * v_m * v_m / (2.0 * D)
    grav = rho_s * _G * math.sin(theta_rad)
    N_Fr = v_m * v_m / (_G * D)
    return BeggsBrillState(regime="gray", holdup=HL, C_L=C_L, N_Fr=N_Fr, N_LV=0.0,
                           rho_ns=rho_ns, rho_s=rho_s, f_n=f, f_tp=f,
                           dP_dz_friction=fric, dP_dz_gravity=grav, dP_dz_total=fric + grav)


TUBING_CORRELATIONS = ("beggs_brill", "hagedorn_brown", "gray")
