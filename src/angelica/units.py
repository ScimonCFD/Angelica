"""
Oil-field units and standard-volume rates (roadmap item 7.1).

Conversion factors are the exact definitions; NIST SP 811 (2008), Appendix
B.8, prints them to 7 significant figures and the tests compare them
(``tests/test_units.py``):

  foot          = 0.3048 m                    (international foot, exact)
  inch          = 0.0254 m                    (exact)
  pound (avdp)  = 0.45359237 kg               (exact)
  pound-force   = 0.45359237 kg × 9.80665 m/s² (standard gravity, exact)
  psi           = lbf / in²  = 6 894.757 293 168 … Pa
  U.S. gallon   = 231 in³    = 3.785 411 784 × 10⁻³ m³ (exact)
  barrel (oil)  = 42 gallons = 0.158 987 294 928 m³ (exact)
  cubic foot    = 0.3048³ m³ = 0.028 316 846 592 m³ (exact)
  centipoise    = 10⁻³ Pa·s
  °F            : t/°C = (t/°F − 32)/1.8;  psig = psia − 14.695 95 (1 atm)

Standard conditions.  Angelica's standard cubic metre (Sm³) is at 60 °F
(15.56 °C) and 1 atm (101 325 Pa) in the black-oil and compositional
models — the conditions of the standard cubic foot (scf).  So
scf = Sm³ / 0.3048³ with no correction for temperature or pressure, and a
stock-tank barrel (STB) = Sm³ of oil / 0.158 987 294 928.
"""
from __future__ import annotations

FT_M = 0.3048
IN_M = 0.0254
LB_KG = 0.45359237
G0 = 9.80665
LBF_N = LB_KG * G0
PSI_PA = LBF_N / IN_M ** 2
GAL_M3 = 231.0 * IN_M ** 3
BBL_M3 = 42.0 * GAL_M3
FT3_M3 = FT_M ** 3
CP_PA_S = 1.0e-3
ATM_PA = 101_325.0
ATM_PSI = ATM_PA / PSI_PA
DAY_S = 86_400.0


def c_to_f(t_c: float) -> float:
    return t_c * 1.8 + 32.0


def f_to_c(t_f: float) -> float:
    return (t_f - 32.0) / 1.8


def pa_to_psia(p_pa: float) -> float:
    return p_pa / PSI_PA


def pa_to_psig(p_pa: float) -> float:
    return p_pa / PSI_PA - ATM_PSI


def psig_to_pa(p_psig: float) -> float:
    return (p_psig + ATM_PSI) * PSI_PA


def sm3_to_bbl(v_sm3: float) -> float:
    """Standard m³ of liquid → barrels (stock-tank barrels for oil)."""
    return v_sm3 / BBL_M3


def sm3_to_mmscf(v_sm3: float) -> float:
    """Standard m³ of gas → million standard cubic feet (same 60 °F, 1 atm)."""
    return v_sm3 / FT3_M3 / 1.0e6
