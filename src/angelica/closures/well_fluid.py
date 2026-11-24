"""
Fluid adapters for the well models (``WellInflow`` and ``Well``).

The IPR, the tubing march and the deliverability curve need the same few
things from the produced fluid, whatever its model:

    phase_properties(P, T)  (C_L, ρ_l, ρ_g, μ_l, μ_g, σ) for the tubing
                            correlations (liquid = oil/condensate + water)
    enthalpy(P, T), cp(P, T), dh_dp(P, T), temperature_from_h(h, P, T0)
                            for the first law along the tubing
    mass_per_rate(ipr)      kg produced per standard m³ of the IPR's
                            reference phase (oil, or gas for back-pressure)
    bubble_point_pa(P, T)   for the composite Vogel IPR

Two adapters:

``BlackOilWellFluid`` wraps the black-oil correlations exactly as the well
model used them before (same results).

``CompositionalWellFluid`` wraps the equation of state.  A PT flash costs
10–15 ms, and a deliverability curve needs tens of thousands of property
evaluations, so the well's fluid is flashed once on a (P, T) grid covering
the tubing and interpolated bilinearly — the PVT-table approach of
commercial network simulators.  Pressure nodes from 1 bar to just above the
reservoir pressure, geometric (ratio 1.35) with at most 10 bar between
nodes; temperature nodes every ~6 K from 40 K below the formation
temperature at the wellhead to 10 K above the reservoir's (about 900
flashes, 5–10 s per well and composition).  cp and (∂h/∂P)_T are the
derivatives of the interpolated enthalpy, so the march is consistent with
it.  Against a march flashing at every step the table gives P_wh within
0.05 % and T_wh within 0.01 K (tests/test_well_compositional.py).

Standard conditions: 101 325 Pa and 15.56 °C (60 °F).  For a compositional
stream the standard volumes come from a single-stage flash of the
hydrocarbons to standard conditions: the gas rate is the vapour volume, the
oil (condensate) rate the liquid volume; water adds mass but no volume, as
the WOR does in black-oil.
"""
from __future__ import annotations

import bisect
import math
from functools import cached_property, lru_cache

from .well_inflow import IPR, mass_per_oil_sc_m3, mass_per_rate_unit

P_SC_PA = 101_325.0
T_SC_C = 15.56

# Table resolution of the compositional adapter (see the module docstring and
# tests/test_well_compositional.py for the interpolation error it gives).
_P_RATIO = 1.35        # pressure nodes: geometric from 1 bar ...
_DP_MAX_PA = 10.0e5    # ... with at most 10 bar between nodes (dense gas: h curves with P)
_DT_K = 6.0            # temperature spacing
_T_BELOW_TOP_K = 40.0  # margin below the formation temperature at the wellhead
_T_ABOVE_RES_K = 10.0  # margin above the reservoir temperature


class BlackOilWellFluid:
    """Black-oil produced fluid (see the module docstring)."""

    exact_enthalpy = True       # enthalpy is the fluid's own h(P, T), not a table

    def __init__(self, fluid) -> None:
        from angelica.properties.black_oil import BlackOilFluid

        self.key = ("black_oil", float(fluid.api_gravity), float(fluid.gas_gravity),
                    float(fluid.gor_sc_m3_per_m3), float(fluid.wor_sc_m3_per_m3),
                    str(getattr(fluid, "correlation", "STANDING")),
                    float(getattr(fluid, "gas_specific_heat_j_per_kg_k", 2_200.0)),
                    float(getattr(fluid, "gas_thermal_conductivity_w_per_mk", 0.035)))
        # Rebuilt from the PVT parameters only (as the well model always did),
        # so that equal fluids share one cached curve.
        _, api, gg, gor, wor, corr, gas_cp, gas_k = self.key
        self.fluid = BlackOilFluid(api_gravity=api, gas_gravity=gg, gor_sc_m3_per_m3=gor,
                                   wor_sc_m3_per_m3=wor, correlation=corr,
                                   gas_specific_heat_j_per_kg_k=gas_cp,
                                   gas_thermal_conductivity_w_per_mk=gas_k)
        self.params = self.fluid._pvt_params()

    def __hash__(self) -> int:
        return hash(self.key)

    def __eq__(self, other) -> bool:
        return isinstance(other, BlackOilWellFluid) and other.key == self.key

    # ── tubing ────────────────────────────────────────────────────────────────
    def phase_properties(self, P: float, T: float) -> tuple[float, float, float, float, float, float]:
        from .beggs_brill import phase_properties_from_pvt, surface_tension_from_pvt

        pvt = self.fluid.pvt(P, T)
        return (*phase_properties_from_pvt(pvt), surface_tension_from_pvt(pvt, float(self.fluid.api_gravity)))

    def enthalpy(self, P: float, T: float) -> float:
        from angelica.properties.black_oil import black_oil_enthalpy
        return black_oil_enthalpy(self.params, P, T)

    def cp(self, P: float, T: float) -> float:
        return self.fluid.pvt(P, T).mixture_specific_heat_j_per_kg_k

    def dh_dp(self, P: float, T: float) -> float:
        from angelica.properties.black_oil import _dh_dp
        return _dh_dp(self.params, P, T)

    def temperature_from_h(self, h: float, P: float, T0: float) -> float:
        # two secant steps from the cp-based estimate (already within ~0.1 K)
        h0 = self.enthalpy(P, T0)
        T1 = T0 + (h - h0) / self.cp(P, T0)
        if abs(T1 - T0) < 1e-6:
            return T1
        h1 = self.enthalpy(P, T1)
        if h1 == h0:
            return T1
        T2 = T1 + (h - h1) * (T1 - T0) / (h1 - h0)
        if abs(T2 - T1) < 1e-4:
            return T2
        h2 = self.enthalpy(P, T2)
        return T2 if h2 == h1 else T2 + (h - h2) * (T2 - T1) / (h2 - h1)

    # ── inflow ────────────────────────────────────────────────────────────────
    def mass_per_base_unit_for(self, ipr: IPR) -> float:
        """kg per stock-tank m³ of oil: the curve of a black-oil well is
        tabulated in oil rate whatever the IPR."""
        return mass_per_oil_sc_m3(self.fluid)

    def mass_per_rate(self, ipr: IPR) -> float:
        return mass_per_rate_unit(self.fluid, ipr)

    def bubble_point_pa(self, P: float, T: float) -> float:
        return self.fluid.pvt(P, T).bubble_point_pa


class CompositionalWellFluid:
    """Compositional produced fluid on a PVT table (see the module docstring).

    ``p_max_pa``, ``t_min_c``, ``t_max_c``: table range; None for an IPR-only
    use (``WellInflow``), which needs no tubing properties.
    """

    exact_enthalpy = False

    def __init__(self, component_names: tuple[str, ...], zs: tuple[float, ...], eos_name: str,
                 p_max_pa: float | None = None, t_min_c: float | None = None,
                 t_max_c: float | None = None) -> None:
        self.component_names = tuple(component_names)
        self.zs = tuple(zs)
        self.eos_name = eos_name
        self.p_max_pa, self.t_min_c, self.t_max_c = p_max_pa, t_min_c, t_max_c
        self.key = ("compositional", self.component_names, self.zs, eos_name, p_max_pa, t_min_c, t_max_c)

    def __hash__(self) -> int:
        return hash(self.key)

    def __eq__(self, other) -> bool:
        return isinstance(other, CompositionalWellFluid) and other.key == self.key

    # ── composition helpers ───────────────────────────────────────────────────
    @cached_property
    def _split(self):
        """(hydrocarbon names, normalised hydrocarbon zs, HC mole fraction, water mole fraction)."""
        from angelica.properties.free_water import find_water_index

        w = find_water_index(self.component_names)
        if w is None:
            return self.component_names, self.zs, 1.0, 0.0
        names = tuple(n for i, n in enumerate(self.component_names) if i != w)
        raw = tuple(z for i, z in enumerate(self.zs) if i != w)
        s = sum(raw)
        return names, tuple(z / s for z in raw) if s > 0 else raw, s, self.zs[w]

    @cached_property
    def _standard_volumes(self) -> tuple[float, float, float]:
        """(kg, gas m³, oil m³) at standard conditions per mole of feed.

        Single-stage flash of the hydrocarbons to standard conditions.  The
        phase volumes use the Peneloux volume translation (the correction the
        surface-tension model already applies): untranslated PR stock-tank
        liquid volumes are 10–15 % too large for heavy oils, which would
        overstate oil rates by as much."""
        from angelica.properties.compositional_fluid import _get_flash_obj, _peneloux_molar_density, _tp_flash

        names, hc_zs, s_hc, _z_w = self._split
        mws = _molar_masses(self.component_names)
        mass = sum(z * m for z, m in zip(self.zs, mws)) / 1000.0
        if s_hc <= 0.0:
            return mass, 0.0, 0.0
        const = _get_flash_obj(names, self.eos_name).constants
        res = _tp_flash(names, P_SC_PA, T_SC_C, hc_zs, self.eos_name)
        phases = list(res.phases)
        if len(phases) == 1:
            ph = phases[0]
            v = s_hc / (_peneloux_molar_density(const, ph) * 1e6)          # m³ per mole of feed
            gas_like = float(ph.rho_mass()) < 100.0
            return (mass, v, 0.0) if gas_like else (mass, 0.0, v)
        i_g = min(range(len(phases)), key=lambda i: phases[i].rho_mass())
        v_gas = v_oil = 0.0
        for i, ph in enumerate(phases):
            v = s_hc * float(res.betas[i]) / (_peneloux_molar_density(const, ph) * 1e6)
            if i == i_g:
                v_gas += v
            else:
                v_oil += v
        return mass, v_gas, v_oil

    @cached_property
    def standard_gas_composition(self) -> tuple[float, ...]:
        """Mole fractions (all components, water 0) of the gas of the stream
        flashed to standard conditions — the default lift gas."""
        from angelica.properties.compositional_fluid import _tp_flash
        from angelica.properties.free_water import find_water_index

        names, hc_zs, _s, _w = self._split
        res = _tp_flash(names, P_SC_PA, T_SC_C, hc_zs, self.eos_name)
        gas = min(res.phases, key=lambda ph: ph.rho_mass())
        ys = list(gas.zs)
        w = find_water_index(self.component_names)
        if w is not None:
            ys.insert(w, 0.0)
        s = sum(ys)
        return tuple(y / s for y in ys)

    @property
    def stock_tank_oil_density_kg_per_m3(self) -> float:
        """Density of the standard-conditions liquid (Peneloux-translated)."""
        from angelica.properties.compositional_fluid import _get_flash_obj, _peneloux_molar_density, _tp_flash

        names, hc_zs, _s, _w = self._split
        const = _get_flash_obj(names, self.eos_name).constants
        res = _tp_flash(names, P_SC_PA, T_SC_C, hc_zs, self.eos_name)
        liq = max(res.phases, key=lambda ph: ph.rho_mass())
        return float(liq.MW()) / 1000.0 * _peneloux_molar_density(const, liq) * 1e6

    @property
    def gas_oil_ratio_sc(self) -> float:
        """Standard gas volume per standard oil volume (inf for a dry gas)."""
        _m, v_g, v_o = self._standard_volumes
        return v_g / v_o if v_o > 0.0 else math.inf

    # ── inflow ────────────────────────────────────────────────────────────────
    def mass_per_rate(self, ipr: IPR) -> float:
        mass, v_gas, v_oil = self._standard_volumes
        if ipr.rate_basis == "gas":
            if v_gas <= 0.0:
                raise ValueError("The back-pressure IPR gives a gas rate, but this stream has no gas "
                                 "at standard conditions.")
            return mass / v_gas
        if v_oil <= 1e-12 * max(v_gas, 1e-30):
            raise ValueError(f"The {ipr.model!r} IPR gives a stock-tank liquid rate, but this stream "
                             "has no liquid at standard conditions: use ipr_model='backpressure' "
                             "for a gas well.")
        return mass / v_oil

    def mass_per_base_unit_for(self, ipr: IPR) -> float:
        """The curve of a compositional well is tabulated in the rate of its
        IPR's reference phase."""
        return self.mass_per_rate(ipr)

    def bubble_point_pa(self, P: float, T: float) -> float:
        """Bubble-point pressure at T (VF = 0 flash of the hydrocarbons); inf
        when there is none (a gas), so the composite Vogel IPR then reduces
        to its saturated form below P_r."""
        return _bubble_point(self._split[0], self._split[1], self.eos_name, round(T, 2))

    # ── PVT table ─────────────────────────────────────────────────────────────
    @cached_property
    def _table(self):
        if self.p_max_pa is None or self.t_min_c is None or self.t_max_c is None:
            raise ValueError("CompositionalWellFluid has no table range (IPR-only use).")
        return _build_table(self.component_names, self.zs, self.eos_name,
                            self.p_max_pa, self.t_min_c, self.t_max_c)

    def _cell(self, P: float, T: float):
        p_nodes, t_nodes, rows = self._table
        i = min(max(bisect.bisect_right(p_nodes, P) - 1, 0), len(p_nodes) - 2)
        j = min(max(bisect.bisect_right(t_nodes, T) - 1, 0), len(t_nodes) - 2)
        u = (P - p_nodes[i]) / (p_nodes[i + 1] - p_nodes[i])
        v = (T - t_nodes[j]) / (t_nodes[j + 1] - t_nodes[j])
        return i, j, u, v, rows

    def _interp(self, P: float, T: float, k: int) -> float:
        i, j, u, v, rows = self._cell(P, T)
        a, b = rows[i][j][k], rows[i][j + 1][k]
        c, d = rows[i + 1][j][k], rows[i + 1][j + 1][k]
        return (1 - u) * ((1 - v) * a + v * b) + u * ((1 - v) * c + v * d)

    def phase_properties(self, P: float, T: float) -> tuple[float, float, float, float, float, float]:
        i, j, u, v, rows = self._cell(P, T)
        out = []
        for k in range(6):
            a, b = rows[i][j][k], rows[i][j + 1][k]
            c, d = rows[i + 1][j][k], rows[i + 1][j + 1][k]
            out.append((1 - u) * ((1 - v) * a + v * b) + u * ((1 - v) * c + v * d))
        C_L = min(max(out[0], 0.0), 1.0)
        return (C_L, max(out[1], 1e-3), max(out[2], 1e-3), max(out[3], 1e-7), max(out[4], 1e-7),
                max(out[5], 1e-5))

    def enthalpy(self, P: float, T: float) -> float:
        return self._interp(P, T, 6)

    def cp(self, P: float, T: float) -> float:
        i, j, u, _v, rows = self._cell(P, T)
        t_nodes = self._table[1]
        dh = ((1 - u) * (rows[i][j + 1][6] - rows[i][j][6])
              + u * (rows[i + 1][j + 1][6] - rows[i + 1][j][6]))
        return max(dh / (t_nodes[j + 1] - t_nodes[j]), 100.0)

    def dh_dp(self, P: float, T: float) -> float:
        i, j, _u, v, rows = self._cell(P, T)
        p_nodes = self._table[0]
        dh = ((1 - v) * (rows[i + 1][j][6] - rows[i][j][6])
              + v * (rows[i + 1][j + 1][6] - rows[i][j + 1][6]))
        return dh / (p_nodes[i + 1] - p_nodes[i])

    def temperature_from_h(self, h: float, P: float, T0: float) -> float:
        """Exact inverse of the interpolated h(P, ·), piecewise linear in T."""
        p_nodes, t_nodes, rows = self._table
        i = min(max(bisect.bisect_right(p_nodes, P) - 1, 0), len(p_nodes) - 2)
        u = (P - p_nodes[i]) / (p_nodes[i + 1] - p_nodes[i])
        H = [(1 - u) * rows[i][j][6] + u * rows[i + 1][j][6] for j in range(len(t_nodes))]
        n = len(t_nodes)
        for j in range(n - 1):
            lo, hi = H[j], H[j + 1]
            if (lo - h) * (hi - h) <= 0.0 and hi != lo:
                return t_nodes[j] + (h - lo) / (hi - lo) * (t_nodes[j + 1] - t_nodes[j])
        # outside the table: extrapolate the nearest end segment
        j = 0 if abs(h - H[0]) < abs(h - H[-1]) else n - 2
        lo, hi = H[j], H[j + 1]
        if hi == lo:
            return T0
        return t_nodes[j] + (h - lo) / (hi - lo) * (t_nodes[j + 1] - t_nodes[j])


@lru_cache(maxsize=64)
def _molar_masses(component_names: tuple[str, ...]) -> tuple[float, ...]:
    from thermo import ChemicalConstantsPackage
    constants, _ = ChemicalConstantsPackage.from_IDs(list(component_names))
    return tuple(constants.MWs)


@lru_cache(maxsize=256)
def _bubble_point(hc_names: tuple[str, ...], hc_zs: tuple[float, ...], eos_name: str, T_c: float) -> float:
    from angelica.properties.compositional_fluid import _get_flash_obj

    try:
        res = _get_flash_obj(hc_names, eos_name).flash(T=T_c + 273.15, VF=0.0, zs=list(hc_zs))
        P = float(res.P)
        return P if math.isfinite(P) and P > 0.0 else math.inf
    except Exception:      # no bubble point at this temperature (a gas)
        return math.inf


@lru_cache(maxsize=64)
def _build_table(component_names, zs, eos_name, p_max, t_min, t_max):
    from angelica.properties.compositional_fluid import _beggs_brill_properties, _mixture_enthalpy

    p_nodes_l = [1.0e5]
    while p_nodes_l[-1] < p_max:       # geometric at low pressure, then steps of _DP_MAX_PA
        p_nodes_l.append(round(min(p_nodes_l[-1] * _P_RATIO, p_nodes_l[-1] + _DP_MAX_PA), -1))
    p_nodes = tuple(p_nodes_l)
    n_t = max(4, int(math.ceil((t_max - t_min) / _DT_K)) + 1)
    t_nodes = tuple(round(t_min + (t_max - t_min) * k / (n_t - 1), 3) for k in range(n_t))
    rows = []
    for P in p_nodes:
        row = []
        for T in t_nodes:
            props = _beggs_brill_properties(component_names, P, T, zs, eos_name)
            row.append((*props, _mixture_enthalpy(component_names, P, T, zs, eos_name)))
        rows.append(tuple(row))
    return p_nodes, t_nodes, tuple(rows)


def well_fluid_for(fluid, link_state=None, component=None, p_res: float | None = None,
                   t_res: float | None = None):
    """Adapter for the fluid of a well link; None if the fluid model is not
    supported (neither black-oil nor compositional).

    For a compositional ``Well`` the table range follows the reservoir
    pressure and temperature and the tubing's geothermal drop."""
    if fluid is None:
        return None
    if hasattr(fluid, "pvt"):
        return BlackOilWellFluid(fluid)
    from angelica.properties.compositional_fluid import _ZS_DIGITS, CompositionalFluid

    if not isinstance(fluid, CompositionalFluid):
        return None
    zs = tuple(round(z, _ZS_DIGITS) for z in fluid._get_zs(link_state))
    if component is None or not getattr(component, "is_well_model", False) or p_res is None or t_res is None:
        return CompositionalWellFluid(fluid.component_names, zs, fluid.eos_name)
    t_top = t_res - float(component.geothermal_gradient_k_per_m) * max(float(component.tubing_height_change_m), 0.0)
    return CompositionalWellFluid(
        fluid.component_names, zs, fluid.eos_name,
        p_max_pa=round(1.05 * p_res, -5) + 1e5,
        t_min_c=round(t_top - _T_BELOW_TOP_K, 0),
        t_max_c=round(t_res + _T_ABOVE_RES_K, 0),
    )


# ── gas lift ──────────────────────────────────────────────────────────────────

class GasLiftMixture:
    """Produced fluid plus lift gas above the injection point, mixed without
    mass transfer: the lift gas stays gas.  ``w`` is the lift-gas mass
    fraction.  Phase volumes add (liquid with liquid, gas with gas);
    enthalpy, cp and (∂h/∂P)_T are mass-weighted."""

    exact_enthalpy = True

    def __init__(self, base, gas, w: float) -> None:
        self.base, self.gas, self.w = base, gas, min(max(w, 0.0), 1.0)
        self.key = ("lift_mix", base.key, gas.key, round(self.w, 12))

    def __hash__(self) -> int:
        return hash(self.key)

    def __eq__(self, other) -> bool:
        return isinstance(other, GasLiftMixture) and other.key == self.key

    def phase_properties(self, P: float, T: float) -> tuple[float, float, float, float, float, float]:
        w = self.w
        parts = []
        if w < 1.0:
            parts.append((1.0 - w, self.base.phase_properties(P, T)))
        if w > 0.0:
            parts.append((w, self.gas.phase_properties(P, T)))
        v_l = v_g = m_l = m_g = mu_l = mu_g = 0.0
        sigma = parts[0][1][5]
        for mass, (C_L, rho_l, rho_g, m_ul, m_ug, _s) in parts:
            v = mass / (C_L * rho_l + (1.0 - C_L) * rho_g)
            vl, vg = C_L * v, (1.0 - C_L) * v
            v_l += vl
            v_g += vg
            m_l += vl * rho_l
            m_g += vg * rho_g
            mu_l += vl * m_ul
            mu_g += vg * m_ug
        base_props = parts[0][1]
        rho_l = m_l / v_l if v_l > 0.0 else base_props[1]
        rho_g = m_g / v_g if v_g > 0.0 else base_props[2]
        return (v_l / (v_l + v_g), rho_l, rho_g,
                mu_l / v_l if v_l > 0.0 else base_props[3],
                mu_g / v_g if v_g > 0.0 else base_props[4], sigma)

    def _mix(self, f, P: float, T: float) -> float:
        a = f(self.base, P, T) if self.w < 1.0 else 0.0
        b = f(self.gas, P, T) if self.w > 0.0 else 0.0
        return (1.0 - self.w) * a + self.w * b

    def enthalpy(self, P: float, T: float) -> float:
        return self._mix(lambda fl, p, t: fl.enthalpy(p, t), P, T)

    def cp(self, P: float, T: float) -> float:
        return self._mix(lambda fl, p, t: fl.cp(p, t), P, T)

    def dh_dp(self, P: float, T: float) -> float:
        return self._mix(lambda fl, p, t: fl.dh_dp(p, t), P, T)

    def temperature_from_h(self, h: float, P: float, T0: float) -> float:
        T = T0
        for _ in range(6):
            f = self.enthalpy(P, T) - h
            if abs(f) < 0.01:
                break
            T -= f / self.cp(P, T)
        return T


def lift_gas_fluid(wf, component, fluid=None, link_state=None):
    """Adapter of a ``Well``'s lift gas, built like the produced fluid's."""
    if isinstance(wf, BlackOilWellFluid):
        from angelica.properties.black_oil import _GAS_ONLY_GOR, BlackOilFluid

        _, api, gg, _gor, _wor, corr, gas_cp, gas_k = wf.key
        gg_lift = float(component.gas_lift_gas_gravity)
        if math.isnan(gg_lift):
            gg_lift = gg
        return BlackOilWellFluid(BlackOilFluid(
            api_gravity=api, gas_gravity=gg_lift, gor_sc_m3_per_m3=_GAS_ONLY_GOR, wor_sc_m3_per_m3=0.0,
            correlation=corr, gas_specific_heat_j_per_kg_k=gas_cp, gas_thermal_conductivity_w_per_mk=gas_k))
    zs = tuple(component.gas_lift_zs) or wf.standard_gas_composition
    if len(zs) != len(wf.component_names):
        raise ValueError(f"Well gas_lift_zs has {len(zs)} mole fractions; the fluid has "
                         f"{len(wf.component_names)} components.")
    return CompositionalWellFluid(wf.component_names, tuple(round(z, 6) for z in zs), wf.eos_name,
                                  wf.p_max_pa, wf.t_min_c, wf.t_max_c)


def lift_gas_mass_per_sc_m3(gas) -> float:
    """kg per standard m³ of lift gas."""
    if isinstance(gas, BlackOilWellFluid):
        from angelica.properties.black_oil import _AIR_DENSITY_SC_KG_M3
        return gas.key[2] * _AIR_DENSITY_SC_KG_M3
    return gas.mass_per_rate(IPR("backpressure", coefficient=1.0))
