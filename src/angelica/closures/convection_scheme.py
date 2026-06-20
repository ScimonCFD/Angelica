from __future__ import annotations

import math
from abc import ABC, abstractmethod


class ConvectionScheme(ABC):
    """Interface for FV convection-diffusion face interpolation schemes.

    Parameters passed to face_coefficients:
        F  — convective strength at the face: ṁ · c_p  (W/K, signed)
             positive means flow goes from W-node toward E-node.
        D  — diffusive conductance at the face: k · A / δx  (W/K, positive)

    Returns (a_W_contribution, a_E_contribution) to add to the FV coefficients
    for the west and east neighbours of node P.
    """

    @abstractmethod
    def face_coefficients(self, F: float, D: float) -> tuple[float, float]:
        """Return (a_W, a_E) contributions from one face pair (west, east)."""


class UpwindScheme(ConvectionScheme):
    """First-order upwind differencing scheme.

    Stable for all Péclet numbers. Introduces numerical diffusion proportional
    to |F|/2 per face, which is negligible when Pe >> 1 (typical pipe flow).
    """

    def face_coefficients(self, F: float, D: float) -> tuple[float, float]:
        a_W = D + max(F, 0.0)
        a_E = D + max(-F, 0.0)
        return a_W, a_E


class HybridScheme(ConvectionScheme):
    """Hybrid scheme (Spalding, 1972).

    Switches between central differencing (|Pe| ≤ 2) and upwind (|Pe| > 2).
    Second-order accurate at low Péclet numbers, first-order at high Pe.
    """

    def face_coefficients(self, F: float, D: float) -> tuple[float, float]:
        a_W = max(F, D - F / 2.0, 0.0)
        a_E = max(-F, D - F / 2.0, 0.0)
        return a_W, a_E


class PowerLawScheme(ConvectionScheme):
    """Power-law scheme (Patankar, 1980).

    Approximates the exact exponential profile. Approaches upwind for |Pe| > 10
    and central differencing for |Pe| → 0.  Recommended default in Malalasekera.
    """

    def face_coefficients(self, F: float, D: float) -> tuple[float, float]:
        pe = F / D if D > 0.0 else math.copysign(1e30, F)
        a_W = D * max(0.0, (1.0 - 0.1 * abs(pe)) ** 5) + max(F, 0.0)
        a_E = D * max(0.0, (1.0 - 0.1 * abs(pe)) ** 5) + max(-F, 0.0)
        return a_W, a_E
