"""
Hagedorn-Brown (modified) and Gray tubing correlations (v1.21.0).

Standing of the checks
----------------------
No published worked example (Brill & Mukherjee, Economides) was accessible,
and installing an independent implementation (pyResToolbox) for a
cross-check was not permitted in this environment.  These tests are
therefore CONSISTENCY checks: a step-by-step hand calculation of the
equations as documented in ``closures/well_correlations.py`` (equations from
the whitson wiki and chart-fit coefficients from pyResToolbox — both read in
secondary sources), single-phase limits, and physical orderings.
NOT VALIDATED against measured well data.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures.beggs_brill import _darcy_friction_factor  # noqa: E402
from angelica.closures.well_correlations import gray_gradient, hagedorn_brown_gradient  # noqa: E402

G = 9.81
FT = 0.3048
# 6-in. tubing, oil-gas flow (the data set used for the Beggs-Brill check)
D, VSL, VSG = 6 * 0.0254, 3.97 * FT, 3.86 * FT
RL, RG = 47.61 * 16.0185, 5.88 * 16.0185
MUL, MUG, SIG, EPS = 0.97e-3, 0.016e-3, 8.41e-3, 0.00006 * FT
P = 2000 * 6894.757        # 2000 psia


def test_hagedorn_brown_hand_calculation():
    vm = VSL + VSG
    lam = VSL / vm
    # λ_g = 0.493; L_B = 1.071 − 0.2218·vm²/D (ft units) is negative → 0.13: not bubble flow
    grp = (RL / (G * SIG)) ** 0.25
    nlv, ngv = VSL * grp, VSG * grp
    nd = D * (RL * G / SIG) ** 0.5
    nl = MUL * (G / (RL * SIG**3)) ** 0.25
    x = math.log10(min(max(nl, 0.002), 0.4)) + 3
    cnl = 10 ** sum(c * x**i for i, c in enumerate((-2.69851, 0.1584095, -0.5509976, 0.5478492, -0.1219458)))
    phi1 = nlv * (P / 101325.0) ** 0.1 * cnl / (ngv**0.575 * nd)
    u = min(math.log10(phi1) + 6, 3.6722855)
    hl_psi = sum(c * u**i for i, c in enumerate((-0.10306578, 0.617774, -0.632946, 0.29598, -0.0401)))
    phi2 = ngv * nl**0.38 / nd**2.14
    psi = 1.0 if phi2 <= 0.001 else sum(
        c * min(max(phi2, 0.012), 0.09) ** i for i, c in enumerate((0.9116257, -4.821756, 1232.25, -22253.58, 116174.3)))
    hl = min(max(psi * hl_psi, lam), 1.0)
    rho_s = hl * RL + (1 - hl) * RG
    gmass = (lam * RL + (1 - lam) * RG) * vm
    re = gmass * D / (MUL**hl * MUG ** (1 - hl))
    f = _darcy_friction_factor(re, EPS / D)
    expected = rho_s * G + f * gmass**2 / (2 * rho_s * D)
    st = hagedorn_brown_gradient(vm, lam, RL, RG, MUL, MUG, D, EPS, math.pi / 2, SIG, P)
    assert st.regime == "hagedorn-brown"
    assert st.holdup == pytest.approx(hl, rel=1e-12)
    assert st.dP_dz_total == pytest.approx(expected, rel=1e-12)
    assert lam <= st.holdup <= 1.0


def test_hagedorn_brown_bubble_flow_griffith():
    """Little gas at low velocity: Griffith bubble flow, holdup from the slip
    velocity 0.8 ft/s, friction on the liquid alone."""
    vsl, vsg = 1.0 * FT, 0.1 * FT
    vm, lam = vsl + vsg, vsl / (vsl + vsg)
    st = hagedorn_brown_gradient(vm, lam, RL, RG, MUL, MUG, D, EPS, math.pi / 2, SIG, P)
    vs = 0.8 * FT
    a = 1 + vm / vs
    hl = 1 - 0.5 * (a - math.sqrt(a * a - 4 * vsg / vs))
    assert st.regime == "bubble" and st.holdup == pytest.approx(hl, rel=1e-12)
    assert st.holdup > lam                                  # gas slips upward


def test_gray_hand_calculation():
    vm = VSL + VSG
    lam = VSL / vm
    rho_ns, mu_ns = lam * RL + (1 - lam) * RG, lam * MUL + (1 - lam) * MUG
    R = VSL / VSG
    nv = rho_ns**2 * vm**2 / (G * SIG * (RL - RG))
    nd = G * (RL - RG) * D * D / SIG
    B = 0.0814 * (1 - 0.0554 * math.log(1 + 730 * R / (1 + R)))
    A = -2.314 * (nv * (1 + 205 / nd)) ** B
    hl = max(1 - (1 - math.exp(A)) / (R + 1), lam)
    k0 = 28.5 * SIG / (rho_ns * vm * vm)
    ke = max(k0 if R >= 0.007 else EPS + R * (k0 - EPS) / 0.007, 2.77e-5 * FT)
    f = _darcy_friction_factor(rho_ns * vm * D / mu_ns, ke / D)
    expected = (hl * RL + (1 - hl) * RG) * G + f * rho_ns * vm * vm / (2 * D)
    st = gray_gradient(vm, lam, RL, RG, MUL, MUG, D, EPS, math.pi / 2, SIG, P)
    assert st.holdup == pytest.approx(hl, rel=1e-12)
    assert st.dP_dz_total == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("fn", [hagedorn_brown_gradient, gray_gradient])
def test_single_phase_gas_limit(fn):
    """Dry gas: no holdup, gravity of the gas, Darcy friction on the gas."""
    v = 8.0
    st = fn(v, 0.0, RL, RG, MUL, MUG, D, EPS, math.pi / 2, SIG, P)
    f = _darcy_friction_factor(RG * v * D / MUG, EPS / D)
    assert st.holdup == 0.0
    assert st.dP_dz_total == pytest.approx(RG * G + f * RG * v * v / (2 * D), rel=1e-9)


def test_well_with_each_correlation():
    from angelica.closures.well_model import well_curve
    from angelica.core.components import Well
    from angelica.properties.black_oil import BlackOilFluid
    f = BlackOilFluid(api_gravity=40, gas_gravity=0.70, gor_sc_m3_per_m3=150, wor_sc_m3_per_m3=0.1)
    curves = {}
    for c in ("beggs_brill", "hagedorn_brown", "gray"):
        w = Well(start_node=1, end_node=2, diameter_m=0.0762, productivity_index_m3_per_day_per_bar=6,
                 tubing_length_m=2100, tubing_height_change_m=2000, tubing_correlation=c)
        curves[c] = well_curve(w, f, 250e5, 85.0)
    for c in ("beggs_brill", "hagedorn_brown"):
        # same static liquid column at zero rate (both have H_L → 1)
        assert curves[c].wellhead_pressures_pa[0] == pytest.approx(curves["beggs_brill"].wellhead_pressures_pa[0], rel=1e-6)
    # Results obtained: peak P_wh 11.43 (BB), 11.80 (HB), 12.69 MPa (Gray) — the
    # correlations differ by a few percent, as expected for a 2 km oil well.
    peaks = [curves[c].max_wellhead_pressure_pa for c in curves]
    assert max(peaks) / min(peaks) < 1.2
    with pytest.raises(ValueError):
        Well(start_node=1, end_node=2, diameter_m=0.07, productivity_index_m3_per_day_per_bar=5,
             tubing_length_m=1000, tubing_correlation="orkiszewski")
