"""
Separator phase split (``SeparatorOutlet`` links).

The streams arriving at a separator node are flashed at the node's pressure
and temperature into gas, oil (hydrocarbon liquid) and free water.  Each
``SeparatorOutlet`` takes the phases it names; the remaining phases — the
gas, or gas + oil for a free-water knockout — leave through the node's other
links.

Black-oil (per stock-tank m³ of oil of the feed):

    oil phase    ρ_o,sc + Rs·ρ_g,sc          (oil with its dissolved gas)
    gas phase    (GOR − Rs)·ρ_g,sc           (free gas)
    water phase  WOR·ρ_w,sc

with Rs = Rs(P, T) of the fluid's correlation, capped at the GOR.  This is
the phase split the black-oil model itself uses in the pipes, so the
separator's gas mass fraction equals the in-situ free-gas mass fraction at
the separator's conditions.  The outlet streams are black-oil compositions:
oil with GOR = Rs; gas with the gas-only GOR of
:func:`~angelica.properties.black_oil.mix_black_oil_streams`; water with a
water-only WOR (``_WATER_ONLY_WOR``).

Compositional: PT flash of the hydrocarbons (the lighter phase is the gas);
free water by the immiscible-water model of the fluid, the water vapour
staying in the gas.

Standard volumes per kg of feed: black-oil from the GOR, Rs and WOR;
compositional by flashing each phase to standard conditions (Peneloux
volumes, ``angelica.closures.well_fluid``) — the gas rate is the vapour of
the gas phase, the oil rate the stock-tank liquid of the oil phase.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

PHASE_SETS = {
    "liquid": frozenset({"oil", "water"}),
    "oil": frozenset({"oil"}),
    "water": frozenset({"water"}),
}
ALL_PHASES = frozenset({"gas", "oil", "water"})

_WATER_ONLY_WOR = 1.0e6     # WOR given to a stream of water with no oil
_WATER_DENSITY_SC = 1_025.0  # as in the black-oil PVT


def gas_outlet_phases(outlet_phases) -> frozenset:
    """Phases leaving through a separator's other links, given its
    ``SeparatorOutlet`` phases."""
    taken: frozenset = frozenset()
    for p in outlet_phases:
        taken |= PHASE_SETS[p]
    return ALL_PHASES - taken


def validate_outlet_phases(node_id: int, phases: list[str]) -> None:
    """A separator has one ``"liquid"`` outlet, or ``"oil"`` + ``"water"``,
    or a ``"water"`` outlet alone (free-water knockout)."""
    if len(set(phases)) != len(phases):
        raise ValueError(f"Separator at node {node_id}: two SeparatorOutlet links take the same phase.")
    if "liquid" in phases and len(phases) > 1:
        raise ValueError(f"Separator at node {node_id}: a 'liquid' outlet already takes oil and water.")
    if "oil" in phases and "water" not in phases:
        raise ValueError(f"Separator at node {node_id}: an 'oil' outlet needs a 'water' outlet "
                         "(three-phase separator); use phase='liquid' for a two-phase one.")


@dataclass(frozen=True)
class SeparatorSplit:
    """Phase split of a separator feed.

    ``fractions``: mass fraction of the feed in each phase ("gas", "oil",
    "water").  ``std_per_kg``: standard volume per kg of feed (m³) of the
    gas phase's gas, the oil phase's stock-tank liquid and the water.
    ``streams``: composition of each combination of phases that leaves
    through one outlet (a ``BlackOilComposition`` or a mole-fraction tuple).
    """

    fractions: dict
    std_per_kg: dict
    streams: dict = field(default_factory=dict)

    def fraction(self, phases) -> float:
        return sum(self.fractions[p] for p in phases)

    def stream(self, phases):
        return self.streams[frozenset(phases)]


# ── black-oil ─────────────────────────────────────────────────────────────────

def black_oil_split(feed, pressure_pa: float, temperature_c: float, correlation: str = "STANDING",
                    gas_cp: float = 2_200.0, gas_k: float = 0.035) -> SeparatorSplit:
    """Split of a black-oil feed (``BlackOilComposition``) at (P, T)."""
    return _black_oil_split(feed.api_gravity, feed.gas_gravity, feed.gor_sc_m3_per_m3,
                            feed.wor_sc_m3_per_m3, round(float(pressure_pa), 0),
                            round(float(temperature_c), 3), correlation, gas_cp, gas_k)


@lru_cache(maxsize=4096)
def _black_oil_split(api, gg, gor, wor, P, T, correlation, gas_cp, gas_k) -> SeparatorSplit:
    from angelica.properties.black_oil import (
        _AIR_DENSITY_SC_KG_M3,
        _GAS_ONLY_GOR,
        BlackOilComposition,
        compute_pvt,
    )
    from angelica.properties.dead_oil import dead_oil_density_kg_per_m3

    rs = gor
    if gor > 0.0:
        rs = min(compute_pvt(P, T, api, gg, gor, wor, correlation, gas_cp, gas_k).rs_m3_per_m3, gor)
    rho_o = dead_oil_density_kg_per_m3(api)
    rho_g = gg * _AIR_DENSITY_SC_KG_M3
    m = {"oil": rho_o + rs * rho_g, "gas": (gor - rs) * rho_g, "water": wor * _WATER_DENSITY_SC}
    total = sum(m.values())

    def comp(phases: frozenset):
        if "oil" in phases:
            return BlackOilComposition(api_gravity=api, gas_gravity=gg,
                                       gor_sc_m3_per_m3=gor if "gas" in phases else rs,
                                       wor_sc_m3_per_m3=wor if "water" in phases else 0.0)
        if "gas" in phases:
            return BlackOilComposition(api_gravity=api, gas_gravity=gg,
                                       gor_sc_m3_per_m3=_GAS_ONLY_GOR, wor_sc_m3_per_m3=0.0)
        return BlackOilComposition(api_gravity=api, gas_gravity=gg,
                                   gor_sc_m3_per_m3=0.0, wor_sc_m3_per_m3=_WATER_ONLY_WOR)

    return SeparatorSplit(
        fractions={p: v / total for p, v in m.items()},
        std_per_kg={"gas": (gor - rs) / total, "oil": 1.0 / total, "water": wor / total},
        streams={s: comp(s) for s in _STREAM_SETS},
    )


# Every phase combination an outlet can take (see validate_outlet_phases).
_STREAM_SETS = (
    frozenset({"gas"}), frozenset({"oil"}), frozenset({"water"}),
    frozenset({"oil", "water"}), frozenset({"gas", "oil"}),
)


# ── compositional ─────────────────────────────────────────────────────────────

def compositional_split(fluid, zs, pressure_pa: float, temperature_c: float) -> SeparatorSplit:
    """Split of a compositional feed (mole fractions ``zs`` of ``fluid``'s
    components) at (P, T)."""
    return _compositional_split(tuple(fluid.component_names), tuple(float(z) for z in zs), fluid.eos_name,
                                round(float(pressure_pa), 0), round(float(temperature_c), 3))


def _hc_split(names, hc_zs, P, T, eos):
    """(β_gas, y, x) of the hydrocarbon PT flash; the lighter phase is the gas."""
    from angelica.properties.compositional_fluid import _get_flash_obj, _tp_flash, _wilson_psat

    try:
        res = _tp_flash(names, P, T, hc_zs, eos)
    except ZeroDivisionError:           # single component (as in _hc_phase_split)
        const = _get_flash_obj(names, eos).constants
        T_K = T + 273.15
        is_gas = T_K >= const.Tcs[0] or P <= _wilson_psat(const.Tcs[0], const.Pcs[0], const.omegas[0], T_K)
        return (1.0 if is_gas else 0.0), hc_zs, hc_zs
    phases = list(res.phases)
    if len(phases) == 1:
        from angelica.properties.compositional_fluid import _gas_like
        const = _get_flash_obj(names, eos).constants
        beta = 1.0 if _gas_like(phases[0], T + 273.15, const) else 0.0
        return beta, hc_zs, hc_zs
    i_g = min(range(len(phases)), key=lambda i: phases[i].rho_mass())
    i_l = max(range(len(phases)), key=lambda i: phases[i].rho_mass())
    beta = float(res.betas[i_g])
    return beta, tuple(float(y) for y in phases[i_g].zs), tuple(float(x) for x in phases[i_l].zs)


@lru_cache(maxsize=1024)
def _compositional_split(names, zs, eos, P, T) -> SeparatorSplit:
    from angelica.properties.free_water import find_water_index, free_water_split, water_liquid_density_kg_m3

    from .well_fluid import T_SC_C, CompositionalWellFluid, _molar_masses

    mws = _molar_masses(names)
    n = len(names)
    w = find_water_index(names)
    hc_idx = [i for i in range(n) if i != w]
    s_hc = sum(zs[i] for i in hc_idx)
    z_w = zs[w] if w is not None else 0.0

    moles = {"gas": [0.0] * n, "oil": [0.0] * n, "water": [0.0] * n}
    if s_hc > 1e-12:
        hc_names = tuple(names[i] for i in hc_idx)
        hc_zs = tuple(zs[i] / s_hc for i in hc_idx)
        beta, y, x = _hc_split(hc_names, hc_zs, P, T, eos)
        for k, i in enumerate(hc_idx):
            moles["gas"][i] = s_hc * beta * y[k]
            moles["oil"][i] = s_hc * (1.0 - beta) * x[k]
        if w is not None:
            n_wl, n_wv = free_water_split(z_w, beta, s_hc, T + 273.15, P)
            moles["gas"][w] = n_wv
            moles["water"][w] = n_wl
    elif w is not None:
        moles["water"][w] = z_w

    mass = {p: sum(m * mw for m, mw in zip(v, mws)) for p, v in moles.items()}   # g per mole of feed
    total = sum(mass.values())

    def zs_of(phases: frozenset) -> tuple[float, ...]:
        v = [sum(moles[p][i] for p in phases) for i in range(n)]
        s = sum(v)
        if s <= 0.0:                    # an empty phase: any valid composition
            return tuple(zs)
        return tuple(c / s for c in v)

    def std_volume(phase: str, which: int) -> float:
        """Standard m³ of gas (which=1) or oil (2) of a phase, per kg of feed."""
        n_p = sum(moles[phase])
        if n_p <= 0.0:
            return 0.0
        kg_per_mol, v_gas, v_oil = CompositionalWellFluid(names, zs_of(frozenset({phase})), eos)._standard_volumes
        return n_p * (v_gas, v_oil)[which - 1] / (total / 1000.0)

    rho_w_sc = water_liquid_density_kg_m3(T_SC_C + 273.15)
    return SeparatorSplit(
        fractions={p: m / total for p, m in mass.items()},
        std_per_kg={"gas": std_volume("gas", 1), "oil": std_volume("oil", 2),
                    "water": mass["water"] / total / rho_w_sc},
        streams={s: zs_of(s) for s in _STREAM_SETS},
    )


# ── flow of an outlet ─────────────────────────────────────────────────────────

def separator_feed_mass_flow(outlet_state) -> float:
    """Mass flow arriving at the separator node of a ``SeparatorOutlet``
    link: inflows of its other links plus internal sources."""
    node = outlet_state.start_node
    m = max(float(node.source_mass_flow_kg_per_s), 0.0)
    for link in outlet_state.separator_feed_links:
        q = float(link.mass_flow_kg_per_s)
        if link.end_node is node:
            m += max(q, 0.0)
        elif link.start_node is node:
            m += max(-q, 0.0)
    return m


def separator_outlet_mass_flow(outlet_state) -> float:
    return outlet_state.separator_fraction * separator_feed_mass_flow(outlet_state)
