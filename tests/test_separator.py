"""
Separator node (``SeparatorOutlet``, roadmap item 5.1).

Standing of the checks
----------------------
- Black-oil split: the separator's gas mass fraction equals the in-situ
  free-gas mass fraction of the black-oil PVT at the same (P, T) — the
  phase split the pipes use (identity, round-off; obtained < 1e-12).
  Standard volumes: the oil rate equals the feed's stock-tank oil rate and
  free gas + dissolved gas = GOR × oil (identities).
- Compositional split: species are conserved (round-off) and the gas mass
  fraction equals the one built from the fluid's own phase-split routine
  (``_hc_phase_split`` + immiscible water), an independent code path
  (obtained < 1e-9).
- Energy: the outlet streams' enthalpies add up to the feed's at the
  separator's (P, T): within 0.0015 K·cp black-oil (obtained 2e-4 to 1.5e-3
  K at 0.5–8 MPa; the gas-only black-oil stream carries 0.1 % phantom oil)
  and 1e-5 K·cp compositional (obtained 6e-6 K).
- Network: outlet flow = fraction × feed and gas = feed − liquid within the
  solver's nodal tolerance; mass and energy balances close.
- Not compared with published separator tests: the split is the existing
  black-oil PVT and PR flash, validated where those are.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures.separator import (  # noqa: E402
    black_oil_split,
    compositional_split,
    gas_outlet_phases,
)
from angelica.core.case import (  # noqa: E402
    InletCompositionBC,
    InletFluidBC,
    NetworkCase,
    PressureBoundary,
    ThermalBoundary,
)
from angelica.core.components import Pipe, SeparatorOutlet  # noqa: E402
from angelica.properties.black_oil import (  # noqa: E402
    BlackOilComposition,
    BlackOilFluid,
    black_oil_enthalpy,
    compute_pvt,
)
from angelica.properties.compositional_fluid import (  # noqa: E402
    CompositionalFluid,
    _hc_phase_split,
    _mixture_enthalpy,
    _reference_enthalpy,
)
from angelica.properties.dead_oil import dead_oil_density_kg_per_m3  # noqa: E402
from angelica.solvers import (  # noqa: E402
    SteadyBlackOilSolver,
    SteadyCompositionalSolver,
    SteadyIsothermalIncompressibleSolver,
)

# The solver's nodal tolerance is 1e-3 of (inflow + outflow), i.e. 2e-3 of the
# flow through a node.
NODAL = 2e-3
CRUDE = {"api_gravity": 35.0, "gas_gravity": 0.70, "gor_sc_m3_per_m3": 120.0, "wor_sc_m3_per_m3": 0.3}
NAMES = ("methane", "ethane", "propane", "n-butane", "n-pentane", "n-decane", "water")
ZS = (0.70, 0.08, 0.05, 0.03, 0.02, 0.07, 0.05)


def _params(c: BlackOilComposition) -> tuple:
    return (c.api_gravity, c.gas_gravity, c.gor_sc_m3_per_m3, c.wor_sc_m3_per_m3, "STANDING", 2_200.0, 0.035)


# ── black-oil split ───────────────────────────────────────────────────────────

@pytest.mark.parametrize("P,T", [(0.5e6, 30.0), (2.26e6, 59.27), (8.0e6, 80.0)])
def test_black_oil_gas_fraction_is_the_in_situ_free_gas(P, T):
    feed = BlackOilComposition(**CRUDE)
    split = black_oil_split(feed, P, T)
    pvt = compute_pvt(P, T, *_params(feed))
    w_gas = pvt.holdup_gas * pvt.density_gas_kg_per_m3 / pvt.mixture_density_kg_per_m3
    w_wat = pvt.holdup_water * pvt.density_water_kg_per_m3 / pvt.mixture_density_kg_per_m3
    assert split.fractions["gas"] == pytest.approx(w_gas, rel=1e-12)
    assert split.fractions["water"] == pytest.approx(w_wat, rel=1e-12)
    assert sum(split.fractions.values()) == pytest.approx(1.0, abs=1e-15)


def test_black_oil_standard_volumes_and_streams():
    feed = BlackOilComposition(**CRUDE)
    P, T = 2.0e6, 50.0
    split = black_oil_split(feed, P, T)
    rs = compute_pvt(P, T, *_params(feed)).rs_m3_per_m3
    m_per_oil = (dead_oil_density_kg_per_m3(35.0) + 120.0 * 0.70 * 1.2250 + 0.3 * 1_025.0)
    assert split.std_per_kg["oil"] == pytest.approx(1.0 / m_per_oil, rel=1e-12)
    assert split.std_per_kg["gas"] + rs * split.std_per_kg["oil"] == pytest.approx(120.0 / m_per_oil, rel=1e-12)
    assert split.std_per_kg["water"] == pytest.approx(0.3 / m_per_oil, rel=1e-12)
    liquid = split.stream({"oil", "water"})
    assert liquid.gor_sc_m3_per_m3 == pytest.approx(rs) and liquid.wor_sc_m3_per_m3 == 0.3
    assert split.stream({"oil"}).wor_sc_m3_per_m3 == 0.0
    assert split.stream({"gas", "oil"}).gor_sc_m3_per_m3 == 120.0
    assert gas_outlet_phases(["water"]) == frozenset({"gas", "oil"})
    assert gas_outlet_phases(["oil", "water"]) == frozenset({"gas"})


@pytest.mark.parametrize("P,T", [(0.5e6, 30.0), (2.26e6, 59.27), (8.0e6, 80.0)])
def test_black_oil_split_conserves_enthalpy(P, T):
    feed = BlackOilComposition(**CRUDE)
    split = black_oil_split(feed, P, T)
    h_feed = black_oil_enthalpy(_params(feed), P, T)
    h_out = sum(split.fractions[p] * black_oil_enthalpy(_params(split.stream({p})), P, T)
                for p in ("gas", "oil", "water"))
    assert abs(h_out - h_feed) / 2_000.0 < 0.0015      # K·cp, cp ≈ 2 kJ/(kg·K)


# ── compositional split ───────────────────────────────────────────────────────

def _fluid():
    return CompositionalFluid(NAMES, ZS)


@pytest.mark.parametrize("P,T", [(10e5, 20.0), (45.88e5, 43.83), (80e5, 60.0)])
def test_compositional_split_conserves_species_and_matches_phase_split(P, T):
    from angelica.properties.free_water import free_water_split

    fl = _fluid()
    mws = fl.component_mws
    split = compositional_split(fl, ZS, P, T)
    mw_feed = sum(z * m for z, m in zip(ZS, mws))
    # species: Σ_phases (mass fraction / MW_phase) · composition = feed moles
    for i in range(len(NAMES)):
        n_i = 0.0
        for p in ("gas", "oil", "water"):
            if split.fractions[p] <= 0.0:
                continue
            zs_p = split.stream({p})
            mw_p = sum(z * m for z, m in zip(zs_p, mws))
            n_i += split.fractions[p] * mw_feed / mw_p * zs_p[i]
        assert n_i == pytest.approx(ZS[i], rel=1e-9, abs=1e-12)
    # independent path: β, MW_g, MW_l of the HC flash + immiscible water
    hc = NAMES[:-1]
    s_hc = sum(ZS[:-1])
    beta, _rg, _rl, mw_g, mw_l, *_ = _hc_phase_split(hc, P, T, tuple(z / s_hc for z in ZS[:-1]), "PR")
    n_wl, n_wv = free_water_split(ZS[-1], beta, s_hc, T + 273.15, P)
    m_gas = s_hc * beta * mw_g + n_wv * mws[-1]
    m_oil = s_hc * (1.0 - beta) * mw_l
    assert split.fractions["gas"] == pytest.approx(m_gas / mw_feed, rel=1e-9)
    assert split.fractions["oil"] == pytest.approx(m_oil / mw_feed, rel=1e-9)
    assert split.fractions["water"] == pytest.approx(n_wl * mws[-1] / mw_feed, rel=1e-9)


@pytest.mark.parametrize("P,T", [(10e5, 20.0), (45.88e5, 43.83), (80e5, 60.0)])
def test_compositional_split_conserves_enthalpy(P, T):
    fl = _fluid()
    split = compositional_split(fl, ZS, P, T)

    def h(zs):
        return _mixture_enthalpy(NAMES, P, T, zs, "PR") - _reference_enthalpy(NAMES, zs, "PR")

    h_out = sum(split.fractions[p] * h(split.stream({p})) for p in ("gas", "oil", "water")
                if split.fractions[p] > 0.0)
    assert abs(h_out - h(ZS)) / 2_500.0 < 1e-5          # K·cp


# ── networks ──────────────────────────────────────────────────────────────────

def _black_oil_case(phases=("liquid",), p_gas=2.0e6, p_liq=0.3e6, per_inlet=True) -> NetworkCase:
    comps = [
        Pipe(component_id="feed", start_node=1, end_node=2, diameter_m=0.1524, length_m=3000.0,
             absolute_roughness_m=4.6e-5, heat_transfer_coefficient_w_per_m2k=2.0, ambient_temperature_c=15.0),
        Pipe(component_id="gas", start_node=2, end_node=4, diameter_m=0.1524, length_m=2000.0,
             absolute_roughness_m=4.6e-5),
    ]
    outs = [PressureBoundary(4, p_gas)]
    for k, ph in enumerate(phases):
        comps.append(SeparatorOutlet(component_id=f"lcv_{ph}", start_node=2, end_node=10 + k,
                                     diameter_m=0.0762, phase=ph))
        comps.append(Pipe(component_id=f"{ph}_line", start_node=10 + k, end_node=20 + k, diameter_m=0.1524,
                          length_m=500.0, absolute_roughness_m=4.6e-5))
        outs.append(PressureBoundary(20 + k, p_liq))
    return NetworkCase(
        name="separator", fluid_model=BlackOilFluid(**CRUDE),
        pressure_inlets=(PressureBoundary(1, 6.0e6),), pressure_outlets=tuple(outs),
        thermal_inlets=(ThermalBoundary(1, 60.0),),
        inlet_fluid_bcs=(InletFluidBC(1, **CRUDE),) if per_inlet else (),
        components=tuple(comps))


def _check_network(res, phases):
    assert res.converged
    flows = {c.label.split(":")[1]: c.mass_flow_kg_per_s for c in res.component_flows}
    sep = res.separators[0]
    feed = flows["feed"]
    assert sep.feed_mass_flow_kg_per_s == pytest.approx(feed, rel=1e-12)
    liquid = sum(flows[f"lcv_{ph}"] for ph in phases)
    assert flows["gas"] + liquid == pytest.approx(feed, rel=NODAL)
    assert flows["gas"] == pytest.approx(sep.gas_mass_flow_kg_per_s
                                         + (sep.oil_mass_flow_kg_per_s if phases == ("water",) else 0.0),
                                         rel=NODAL)
    gb = res.global_balance
    assert abs(gb.mass_inlet_kg_per_s - gb.mass_outlet_kg_per_s) < NODAL * feed
    geb = res.global_energy_balance
    # ṁ·cp·0.01 K, the solver's temperature tolerance (obtained 0.003–0.008 kW)
    assert abs(geb.energy_error_kw) < feed * 2500.0 * 0.01 / 1e3
    return flows, sep


@pytest.mark.parametrize("phases", [("liquid",), ("oil", "water"), ("water",)])
def test_black_oil_separator_network(phases):
    res = SteadyBlackOilSolver().solve(_black_oil_case(phases))
    flows, sep = _check_network(res, phases)
    for ph in phases:
        expected = {"liquid": sep.oil_mass_flow_kg_per_s + sep.water_mass_flow_kg_per_s,
                    "oil": sep.oil_mass_flow_kg_per_s, "water": sep.water_mass_flow_kg_per_s}[ph]
        # the outlet flow is set from the feed of the previous hydraulic
        # iteration: equal within the nodal tolerance (obtained 7e-7)
        assert flows[f"lcv_{ph}"] == pytest.approx(expected, rel=NODAL)
    # the gas fraction is the in-situ free gas at the separator's converged state
    pvt = compute_pvt(sep.pressure_pa, sep.temperature_c, *_params(BlackOilComposition(**CRUDE)))
    w_gas = pvt.holdup_gas * pvt.density_gas_kg_per_m3 / pvt.mixture_density_kg_per_m3
    assert sep.gas_mass_flow_kg_per_s / sep.feed_mass_flow_kg_per_s == pytest.approx(w_gas, rel=1e-5)


def test_black_oil_separator_without_inlet_compositions():
    a = SteadyBlackOilSolver().solve(_black_oil_case(per_inlet=False))
    b = SteadyBlackOilSolver().solve(_black_oil_case(per_inlet=True))
    assert a.converged
    assert a.separators[0].gas_mass_flow_kg_per_s == pytest.approx(b.separators[0].gas_mass_flow_kg_per_s, rel=1e-6)


def test_level_valve_pressure_drop_and_pump_message():
    ok = SteadyBlackOilSolver().solve(_black_oil_case(p_liq=0.3e6)).separators[0]
    assert ok.outlets[0].valve_pressure_drop_pa > 0.0 and ok.message == "ok"
    # liquid delivered against a pressure above the separator's
    high = SteadyBlackOilSolver().solve(_black_oil_case(p_gas=1.0e6, p_liq=1.5e6)).separators[0]
    assert high.outlets[0].valve_pressure_drop_pa < 0.0
    assert "needs a pump" in high.message


def test_separator_at_a_pressure_boundary():
    """The gas leaves into the separator node's own boundary (no gas line)."""
    case = NetworkCase(
        name="sep at boundary", fluid_model=BlackOilFluid(**CRUDE),
        pressure_inlets=(PressureBoundary(1, 6.0e6),),
        pressure_outlets=(PressureBoundary(2, 2.0e6), PressureBoundary(6, 0.3e6)),
        thermal_inlets=(ThermalBoundary(1, 60.0),), inlet_fluid_bcs=(InletFluidBC(1, **CRUDE),),
        components=(Pipe(component_id="feed", start_node=1, end_node=2, diameter_m=0.1524, length_m=3000.0,
                         absolute_roughness_m=4.6e-5),
                    SeparatorOutlet(component_id="lcv", start_node=2, end_node=5, diameter_m=0.0762),
                    Pipe(component_id="oil", start_node=5, end_node=6, diameter_m=0.1524, length_m=500.0,
                         absolute_roughness_m=4.6e-5)))
    res = SteadyBlackOilSolver().solve(case)
    assert res.converged
    sep = res.separators[0]
    assert res.component_flows[1].mass_flow_kg_per_s == pytest.approx(
        sep.oil_mass_flow_kg_per_s + sep.water_mass_flow_kg_per_s, rel=NODAL)
    gb = res.global_balance
    assert abs(gb.mass_inlet_kg_per_s - gb.mass_outlet_kg_per_s) < NODAL * sep.feed_mass_flow_kg_per_s


def _compositional_case(phases=("liquid",)) -> NetworkCase:
    comps = [
        Pipe(component_id="feed", start_node=1, end_node=2, diameter_m=0.1524, length_m=3000.0,
             absolute_roughness_m=4.6e-5, heat_transfer_coefficient_w_per_m2k=2.0, ambient_temperature_c=15.0),
        Pipe(component_id="gas", start_node=2, end_node=4, diameter_m=0.1524, length_m=2000.0,
             absolute_roughness_m=4.6e-5),
    ]
    outs = [PressureBoundary(4, 40e5)]
    for k, ph in enumerate(phases):
        comps.append(SeparatorOutlet(component_id=f"lcv_{ph}", start_node=2, end_node=10 + k,
                                     diameter_m=0.0508, phase=ph))
        comps.append(Pipe(component_id=f"{ph}_line", start_node=10 + k, end_node=20 + k, diameter_m=0.0762,
                          length_m=500.0, absolute_roughness_m=4.6e-5))
        outs.append(PressureBoundary(20 + k, 10e5))
    return NetworkCase(name="separator", fluid_model=_fluid(), pressure_inlets=(PressureBoundary(1, 60e5),),
                       pressure_outlets=tuple(outs), thermal_inlets=(ThermalBoundary(1, 50.0),),
                       inlet_composition_bcs=(InletCompositionBC(1, ZS),), components=tuple(comps))


@pytest.mark.parametrize("phases", [("liquid",), ("oil", "water")])
def test_compositional_separator_network(phases):
    res = SteadyCompositionalSolver().solve(_compositional_case(phases))
    flows, sep = _check_network(res, phases)
    links = {c.label.split(":")[1]: c for c in res.component_flows}
    split = compositional_split(_fluid(), ZS, sep.pressure_pa, sep.temperature_c)
    # The outlet compositions come from the last outer iteration, which stops
    # when temperatures change by less than 0.01 K: tolerance = change of the
    # gas composition over ±0.01 K (obtained 9e-7).
    y_hi = compositional_split(_fluid(), ZS, sep.pressure_pa, sep.temperature_c + 0.01).stream({"gas"})
    y_lo = compositional_split(_fluid(), ZS, sep.pressure_pa, sep.temperature_c - 0.01).stream({"gas"})
    tol = max(abs(a - b) for a, b in zip(y_hi, y_lo))
    assert links["gas"].zs == pytest.approx(split.stream({"gas"}), abs=tol)
    if phases == ("oil", "water"):
        assert links["lcv_water"].zs[-1] == 1.0
        assert links["lcv_oil"].zs[-1] == 0.0
    # the separator temperature is the feed's: separation at (P, T) conserves enthalpy
    assert sep.temperature_c == pytest.approx(res.node_temperatures_c[2])


def test_layout_errors():
    with pytest.raises(ValueError, match="phase"):
        SeparatorOutlet(component_id="x", start_node=1, end_node=2, diameter_m=0.05, phase="gas")
    with pytest.raises(ValueError, match="'water' outlet"):
        SteadyBlackOilSolver().solve(_black_oil_case(("oil",)))
    with pytest.raises(ValueError, match="already takes"):
        SteadyBlackOilSolver().solve(_black_oil_case(("liquid", "water")))
    no_gas = NetworkCase(
        name="no gas outlet", fluid_model=BlackOilFluid(**CRUDE),
        pressure_inlets=(PressureBoundary(1, 6.0e6),), pressure_outlets=(PressureBoundary(6, 0.3e6),),
        thermal_inlets=(ThermalBoundary(1, 60.0),),
        components=(Pipe(component_id="feed", start_node=1, end_node=2, diameter_m=0.15, length_m=100.0),
                    SeparatorOutlet(component_id="lcv", start_node=2, end_node=5, diameter_m=0.05),
                    Pipe(component_id="oil", start_node=5, end_node=6, diameter_m=0.15, length_m=100.0)))
    with pytest.raises(ValueError, match="gas has no outlet"):
        SteadyBlackOilSolver().solve(no_gas)


def test_separator_outlet_rejected_by_single_phase_solvers():
    from angelica.properties.single_component import SingleComponentFluid
    case = NetworkCase(
        name="water", fluid_model=SingleComponentFluid(density_kg_per_m3=998.25, viscosity_pa_s=0.001),
        pressure_inlets=(PressureBoundary(1, 3.0e5),),
        pressure_outlets=(PressureBoundary(3, 1.0e5), PressureBoundary(4, 1.0e5)),
        components=(Pipe(component_id="a", start_node=1, end_node=2, diameter_m=0.1, length_m=10.0),
                    Pipe(component_id="b", start_node=2, end_node=3, diameter_m=0.1, length_m=10.0),
                    SeparatorOutlet(component_id="s", start_node=2, end_node=4, diameter_m=0.05)))
    with pytest.raises(ValueError, match="SeparatorOutlet"):
        SteadyIsothermalIncompressibleSolver().solve(case)


def test_boundary_node_with_throughflow_is_at_the_mixing_temperature():
    """A separator at a pressure boundary (gas leaves through the boundary,
    liquid flows on) must sit at the temperature of the arriving stream, as
    the same separator one 1-m gas line away from the boundary does.  Before
    v1.26.0 the boundary node divided the inflowing enthalpy by the outflow
    only (−256 °C in compositional tutorial 10).  Tolerance 0.01 K, the
    solver's temperature tolerance (obtained 2e-4 K)."""
    def case(at_boundary: bool) -> NetworkCase:
        comps = [Pipe(component_id="feed", start_node=1, end_node=2, diameter_m=0.1524, length_m=3000.0,
                      absolute_roughness_m=4.6e-5, heat_transfer_coefficient_w_per_m2k=2.0,
                      ambient_temperature_c=15.0),
                 SeparatorOutlet(component_id="lcv", start_node=2, end_node=5, diameter_m=0.0762),
                 Pipe(component_id="oil", start_node=5, end_node=6, diameter_m=0.1524, length_m=500.0,
                      absolute_roughness_m=4.6e-5)]
        outs = [PressureBoundary(6, 0.3e6)]
        if at_boundary:
            outs.append(PressureBoundary(2, 2.0e6))
        else:
            comps.append(Pipe(component_id="gas", start_node=2, end_node=4, diameter_m=0.5, length_m=1.0))
            outs.append(PressureBoundary(4, 2.0e6))
        return NetworkCase(name="boundary separator", fluid_model=BlackOilFluid(**CRUDE),
                           pressure_inlets=(PressureBoundary(1, 6.0e6),), pressure_outlets=tuple(outs),
                           thermal_inlets=(ThermalBoundary(1, 60.0),),
                           inlet_fluid_bcs=(InletFluidBC(1, **CRUDE),), components=tuple(comps))

    a = SteadyBlackOilSolver().solve(case(True))
    b = SteadyBlackOilSolver().solve(case(False))
    assert a.converged and b.converged
    assert a.node_temperatures_c[2] == pytest.approx(b.node_temperatures_c[2], abs=0.01)
    assert a.node_temperatures_c[6] == pytest.approx(b.node_temperatures_c[6], abs=0.01)


def test_two_stage_train_with_interior_separators():
    """HP separator → oil LCV → LP separator (an ordinary node with a flare
    line) → liquid LCV → stock tank.  The LP feed is nearly all liquid
    (outlet fraction ≈ 1), the case that diverged while the outlet flow was
    lagged in the pressure correction.  Each outlet carries its fraction of
    its separator's feed (nodal tolerance) and the balances close."""
    def pipe(cid, a, b, D, L):
        return Pipe(component_id=cid, start_node=a, end_node=b, diameter_m=D, length_m=L,
                    absolute_roughness_m=4.6e-5)
    case = NetworkCase(
        name="train", fluid_model=BlackOilFluid(**CRUDE),
        pressure_inlets=(PressureBoundary(1, 6.0e6),),
        pressure_outlets=(PressureBoundary(4, 2.0e6), PressureBoundary(8, 0.3e6), PressureBoundary(10, 0.12e6),
                          PressureBoundary(12, 0.3e6)),
        thermal_inlets=(ThermalBoundary(1, 60.0),), inlet_fluid_bcs=(InletFluidBC(1, **CRUDE),),
        components=(pipe("feed", 1, 2, 0.1524, 3000.0), pipe("gas", 2, 4, 0.1524, 2000.0),
                    SeparatorOutlet(component_id="hp_oil", start_node=2, end_node=5, diameter_m=0.0762, phase="oil"),
                    pipe("oil", 5, 6, 0.1524, 500.0),
                    SeparatorOutlet(component_id="hp_water", start_node=2, end_node=11, diameter_m=0.05,
                                    phase="water"),
                    pipe("water", 11, 12, 0.1016, 300.0),
                    pipe("flare", 6, 8, 0.2032, 100.0),
                    SeparatorOutlet(component_id="lp_liquid", start_node=6, end_node=9, diameter_m=0.0762),
                    pipe("tank", 9, 10, 0.1524, 100.0)))
    res = SteadyBlackOilSolver().solve(case)
    assert res.converged
    flows = {c.label.split(":")[1]: c.mass_flow_kg_per_s for c in res.component_flows}
    seps = {s.node_id: s for s in res.separators}
    hp, lp = seps[2], seps[6]
    assert lp.feed_mass_flow_kg_per_s == pytest.approx(flows["hp_oil"], rel=NODAL)
    assert flows["lp_liquid"] / lp.feed_mass_flow_kg_per_s > 0.95          # the stiff case
    assert flows["lp_liquid"] == pytest.approx(lp.oil_mass_flow_kg_per_s + lp.water_mass_flow_kg_per_s, rel=NODAL)
    assert flows["flare"] == pytest.approx(lp.gas_mass_flow_kg_per_s, rel=NODAL)
    gb = res.global_balance
    assert abs(gb.mass_inlet_kg_per_s - gb.mass_outlet_kg_per_s) < NODAL * hp.feed_mass_flow_kg_per_s
    assert abs(res.global_energy_balance.energy_error_kw) < hp.feed_mass_flow_kg_per_s * 2500.0 * 0.01 / 1e3
