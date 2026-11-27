"""
Gas compressor: thermodynamics of the compression and its hydraulic law.

Compression path
----------------
Along a reversible adiabatic (isentropic) path the Gibbs relation
dh = T·ds + v·dP reduces to dh = v·dP.  The isentropic work from P1 to P2 is
therefore W_s = ∫ v dP along the isentrope, which is marched in ln P
with the fluid's own h(P, T) and v(P, T) (Heun's predictor-corrector, the
temperature at each step from inverting h).  No entropy function is needed,
so the same march serves every fluid with an enthalpy model — gas (ideal or
Peng-Robinson), black-oil, compositional — and the result is an isentrope
of that fluid's model exactly when its h and v are thermodynamically
consistent (checked in the tests against the EOS entropy).

Efficiency (one of):
  - isentropic:  h2 = h1 + W_s / η_s
  - polytropic:  the march itself with dh = v·dP / η_p (the actual path,
    each small step at efficiency η_p).

Hydraulics
----------
The compressor imposes a pressure rise on the flow through it: a ratio
P_out/P_in, a rise ΔP, or a shaft power W (the rise is then the one at
which ṁ·(h2 − h1) = W, up to ``max_pressure_ratio``).  In the
pressure-correction solver it is a nearly transparent link around that
lagged rise; it does not run backwards (a compressor with the outlet above
what it can deliver passes no flow, like a check valve).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from types import SimpleNamespace

from .pressure_drop import PressureDropCorrelation

_K = 1.0             # loss coefficient of the link around the imposed rise (made up by the rise)
_V_FLOOR = 1e-4      # m/s, velocity floor of the coupling (see heat_source.py)
_CLOSED_COUPLING = 1e-12


@dataclass(frozen=True)
class CompressionResult:
    p_in_pa: float
    t_in_c: float
    p_out_pa: float
    t_out_c: float
    t_out_isentropic_c: float
    isentropic_work_j_per_kg: float      # h2s − h1
    actual_work_j_per_kg: float          # h2 − h1


class _StateAt:
    """A link state seen at a given (P, T): what the fluid models read."""

    __slots__ = ("start_node", "end_node", "temperature_c", "zs", "component")

    def __init__(self, link, P: float, T: float) -> None:
        node = SimpleNamespace(pressure_pa=P)
        self.start_node = node
        self.end_node = node
        self.temperature_c = T
        self.zs = getattr(link, "zs", ())
        self.component = getattr(link, "component", None)


class CompressorFluid:
    """h(P, T), v(P, T) and T(h, P) of the fluid in a link, at any (P, T)."""

    def __init__(self, fluid_model, link) -> None:
        from .choke import resolve_link_fluid
        self.fluid = resolve_link_fluid(fluid_model, link)
        self.link = link

    def h(self, P: float, T: float) -> float:
        f = self.fluid
        if hasattr(f, "enthalpy_for_link"):
            return float(f.enthalpy_for_link(_StateAt(self.link, P, T), P, T))
        if hasattr(f, "enthalpy_j_per_kg"):
            return float(f.enthalpy_j_per_kg(P, T))
        raise ValueError("A compressor needs a fluid with an enthalpy model "
                         "(compressible gas, black-oil or compositional).")

    def v(self, P: float, T: float) -> float:
        return 1.0 / float(self.fluid.density_for_link(_StateAt(self.link, P, T)))

    def liquid_at(self, P: float, T: float) -> str | None:
        """A note when the feed is not (almost) all gas at (P, T): the
        compression model is that of a single-phase fluid's isentrope."""
        f = self.fluid
        if hasattr(f, "pvt"):
            gas = float(f.pvt(P, T).holdup_gas)
            return None if gas >= 0.99 else f"{(1.0 - gas) * 100:.1f} % liquid by volume"
        if hasattr(f, "vapor_fraction_for_link"):
            vf = float(f.vapor_fraction_for_link(_StateAt(self.link, P, T)))
            return None if vf >= 0.99 else f"vapour mole fraction {vf:.3f}"
        return None

    def T_from_h(self, h: float, P: float, T_guess: float) -> float:
        f = self.fluid
        if hasattr(f, "temperature_from_enthalpy_for_link"):
            return float(f.temperature_from_enthalpy_for_link(_StateAt(self.link, P, T_guess), h, P,
                                                               T_guess_c=T_guess))
        return float(f.temperature_from_enthalpy(h, P, T_guess_c=T_guess))


def _march(fl: CompressorFluid, P1: float, T1: float, P2: float, eta: float, n: int) -> tuple[float, float]:
    """(h, T) at P2 after dh = v·dP/η from (P1, T1): Heun's method in ln P
    (dh = P·v·d(ln P)/η; P·v varies far less than v along a compression)."""
    h, T, x = fl.h(P1, T1), T1, math.log(P1)
    dx = (math.log(P2) - x) / n
    for _ in range(n):
        P0, P1_ = math.exp(x), math.exp(x + dx)
        g0 = P0 * fl.v(P0, T)
        h_pred = h + g0 * dx / eta
        T_pred = fl.T_from_h(h_pred, P1_, T)
        g1 = P1_ * fl.v(P1_, T_pred)
        h = h + 0.5 * (g0 + g1) * dx / eta
        x = x + dx
        T = fl.T_from_h(h, P1_, T_pred)
    return h, T


# Results by (fluid object, composition, rounded state): every energy pass and
# every step of the power search asks again.  The fluid object is kept in the
# value so its id cannot be reused while the entry lives.
_CACHE: dict = {}
_CACHE_MAX = 20000


def compress(fl: CompressorFluid, P1: float, T1: float, P2: float, efficiency: float,
             efficiency_type: str = "isentropic", n_steps: int = 20) -> CompressionResult:
    """State after compressing from (P1, T1) to P2 (cached on the state
    rounded to 1 Pa and 1e-4 K)."""
    key = (id(fl.fluid), tuple(round(z, 9) for z in (getattr(fl.link, "zs", ()) or ())),
           round(P1), round(T1, 4), round(P2), efficiency, efficiency_type, n_steps)
    hit = _CACHE.get(key)
    if hit is not None:
        return hit[1]
    res = _compress(fl, float(round(P1)), round(T1, 4), float(round(P2)), efficiency, efficiency_type, n_steps)
    if len(_CACHE) >= _CACHE_MAX:
        _CACHE.clear()
    _CACHE[key] = (fl.fluid, res)
    return res


def _compress(fl: CompressorFluid, P1: float, T1: float, P2: float, efficiency: float,
              efficiency_type: str, n_steps: int) -> CompressionResult:
    h1 = fl.h(P1, T1)
    if P2 <= P1:
        return CompressionResult(P1, T1, P1, T1, T1, 0.0, 0.0)
    h2s, T2s = _march(fl, P1, T1, P2, 1.0, n_steps)
    if efficiency_type == "polytropic":
        h2, T2 = _march(fl, P1, T1, P2, efficiency, n_steps)
    else:
        h2 = h1 + (h2s - h1) / efficiency
        T2 = fl.T_from_h(h2, P2, T2s)
    return CompressionResult(P1, T1, P2, T2, T2s, h2s - h1, h2 - h1)


def outlet_pressure_for_power(fl: CompressorFluid, comp, P1: float, T1: float, m_dot: float) -> float:
    """Outlet pressure at which ṁ·(h2 − h1) equals the shaft power, capped at
    ``max_pressure_ratio`` (zero or tiny flow: the cap)."""
    r_max = float(comp.max_pressure_ratio)
    W = float(comp.power_w)
    m = abs(m_dot)

    def power_at(r: float) -> float:
        res = compress(fl, P1, T1, r * P1, comp.efficiency, comp.efficiency_type, comp.n_steps)
        return m * res.actual_work_j_per_kg

    if W <= 0.0:
        return P1                             # no power, no compression
    g_hi = power_at(r_max) - W if m >= 1e-12 else 1.0
    if m < 1e-12 or g_hi <= 0.0:
        return r_max * P1
    # W(r) increases with r: false position with the Illinois modification.
    lo, g_lo, hi, side = 1.0, -W, r_max, 0
    r = hi
    for _ in range(60):
        r = hi - g_hi * (hi - lo) / (g_hi - g_lo)
        g = power_at(r) - W
        if abs(g) <= 1e-7 * W or hi - lo < 1e-9:
            break
        if g < 0.0:
            lo, g_lo = r, g
            if side == -1:
                g_hi *= 0.5
            side = -1
        else:
            hi, g_hi = r, g
            if side == 1:
                g_lo *= 0.5
            side = 1
    return r * P1


class CompressorModel(PressureDropCorrelation):
    """Hydraulic law of a ``Compressor`` link: a link with K = 1 around the
    imposed pressure rise ``link_state.compressor_rise_pa``, which the solver
    adjusts each outer iteration until the node pressures meet the
    specification (so the K loss is made up); no back-flow.  (K = 0.01, a
    stiffer link, left 0.15 % of mass imbalance at a converged boundary.)"""

    @staticmethod
    def _rise(link_state) -> float:
        # Set by the solver once per outer iteration (lagged) so that the node
        # pressures meet the specification: a ratio taken from the inlet
        # pressure inside the pressure correction fed back on itself (the
        # coupling cannot carry dRise/dP_in) and the flow ran away to
        # 5 000 kg/s; the link's own loss is made up by the same adjustment.
        return float(getattr(link_state, "compressor_rise_pa", 0.0) or 0.0)

    def calculate_velocity(self, link_state, delta_p: float, density: float, viscosity: float,
                           tolerance: float | None = None) -> float:
        dp_eff = delta_p + self._rise(link_state)
        if dp_eff <= 0.0:
            return 0.0                        # cannot deliver against that outlet pressure
        return math.sqrt(2.0 * dp_eff / (_K * density))

    def calculate_coupling(self, link_state, density: float, viscosity: float) -> float:
        dp_eff = (float(link_state.start_node.pressure_pa) - float(link_state.end_node.pressure_pa)
                  + self._rise(link_state))
        if dp_eff <= 0.0 and abs(link_state.velocity_m_per_s) < 1e-12:
            return -_CLOSED_COUPLING
        return -2.0 * link_state.area_m2 / (_K * max(abs(link_state.velocity_m_per_s), _V_FLOOR))
