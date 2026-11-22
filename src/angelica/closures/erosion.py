"""
Erosional velocity, API RP 14E (5th ed., 1991, §2.5):

    V_e = C / √ρ_m        V_e in ft/s, ρ_m in lbm/ft³ (no-slip mixture density
                          at flowing conditions)

C = 100 for continuous service and 125 for intermittent service; for
solids-free, non-corrosive fluids C = 150–200 continuous may be justified.
In SI units V_e [m/s] = 0.3048·C / √(ρ_m [kg/m³] / 16.0185).

A line is usually sized so that the mixture velocity stays below V_e
(ratio v_m / V_e < 1).
"""
from __future__ import annotations

import math

API_14E_C_CONTINUOUS = 100.0
API_14E_C_INTERMITTENT = 125.0
_FT = 0.3048
_KG_M3_PER_LB_FT3 = 16.0185


def erosional_velocity(mixture_density_kg_per_m3: float, c_factor: float = API_14E_C_CONTINUOUS) -> float:
    """API RP 14E erosional velocity (m/s) for a no-slip mixture density (kg/m³)."""
    rho_lb_ft3 = max(mixture_density_kg_per_m3, 1e-9) / _KG_M3_PER_LB_FT3
    return _FT * c_factor / math.sqrt(rho_lb_ft3)
