"""
Tutorial 12 — Wells with production tubing (nodal analysis in a network)
=======================================================================

Tutorial 11 referred each well's inflow to the wellhead.  Here every well
carries its production tubing: 2 100 m of tubing rising 2 000 m from the
bottom hole to the wellhead.  The ``Well`` component combines the reservoir
inflow (IPR) and the tubing, so the network sees each well through its
deliverability curve P_wh(q).

Network
-------

  Reservoir 1 ─ Well (IPR + tubing) ─ WH1 ─ flowline 3 km ─┐
  Reservoir 2 ─ Well (IPR + tubing) ─ WH2 ─ flowline 2 km ─┼─ Manifold ─ trunk 8 km ─ Separator
  Reservoir 3 ─ Well (IPR + tubing) ─ WH3 ─ flowline 4 km ─┘                         1.5 MPa

    Well 1: 40 °API, GOR 150, WOR 0.1, P_r = 25 MPa, 85 °C, J = 6 m³/d/bar
    Well 2: 32 °API, GOR  40, WOR 0.3, P_r = 22 MPa, 80 °C, J = 4 m³/d/bar
    Well 3: 25 °API, GOR  15, WOR 0.8, P_r = 17 MPa, 75 °C, J = 3 m³/d/bar
    Tubing 3.0 in ID (76.2 mm) unless stated; flowlines 4 in, Beggs-Brill.

How the Well component works
----------------------------
For an oil rate q the IPR (composite Vogel) gives the flowing bottom-hole
pressure P_wf.  The tubing pressure drop is then integrated upward with
Beggs-Brill and the local black-oil PVT, giving the wellhead pressure
P_wh(q).  With little flow the liquid holdup is high and the column is heavy;
a little gas lightens it, so P_wh(q) first rises and then falls — the
J-shaped curve of nodal analysis.  The network uses the stable, falling
branch; if the wellhead pressure exceeds the curve's maximum the well is dead.

Solving the tubing as an ordinary network pipe does not work: one pressure
drop can correspond to several flow rates on the J-curve and the pressure
correction jumps between them.

Temperature
-----------
The fluid enters the tubing at its reservoir enthalpy and loses heat to the
formation on the way up (Ramey's model: transient conduction in the rock,
which is at the reservoir temperature at the bottom and 0.03 K/m cooler
towards the surface — here 25 °C at the wellhead of an 85 °C reservoir).
The faster a well flows, the less time the fluid has to cool, so the
wellhead temperature rises with the rate.  The network receives each well at
its flowing wellhead temperature (Part E).  With the default
``tubing_u_w_per_m2k = inf`` the annulus, casing and cement offer no
resistance: the heat loss is an upper bound and the wellhead temperature a
lower bound.

Part F raises the separator pressure until one well dies against the
manifold.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures import BeggsBrillCorrelation
from angelica.closures.well_model import bottomhole_pressure, well_curve, well_tubing_state
from angelica.core.case import InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe, Well
from angelica.properties.black_oil import BlackOilFluid
from angelica.solvers import SteadyBlackOilSolver

# name, P_r (MPa), T (°C), API, gas gravity, GOR, WOR, J (m³/d/bar), flowline length (m)
WELLS = (
    ("Well 1", 25.0, 85.0, 40.0, 0.70, 150.0, 0.1, 6.0, 3000.0),
    ("Well 2", 22.0, 80.0, 32.0, 0.65,  40.0, 0.3, 4.0, 2000.0),
    ("Well 3", 17.0, 75.0, 25.0, 0.60,  15.0, 0.8, 3.0, 4000.0),
)
TUBING_LENGTH, TUBING_RISE = 2100.0, 2000.0
P_SEP = 1.5e6
MANIFOLD, SEPARATOR = 100, 101


def well_fluid(api, gg, gor, wor) -> BlackOilFluid:
    return BlackOilFluid(api_gravity=api, gas_gravity=gg, gor_sc_m3_per_m3=gor, wor_sc_m3_per_m3=wor)


def make_well(i: int, tubing_id: float, correlation: str = "beggs_brill") -> Well:
    name, *_rest = WELLS[i]
    J = WELLS[i][7]
    return Well(component_id=f"well_{name.replace(' ', '')}", start_node=10 * (i + 1),
                end_node=10 * (i + 1) + 1, diameter_m=tubing_id,
                productivity_index_m3_per_day_per_bar=J, ipr_model="vogel",
                tubing_length_m=TUBING_LENGTH, tubing_height_change_m=TUBING_RISE,
                tubing_correlation=correlation)


def build_case(tubing_id: float = 0.0762, correlation: str = "beggs_brill",
               p_sep: float = P_SEP) -> NetworkCase:
    components, inlets, thermal, fluids = [], [], [], []
    for i, (name, P_r, T, api, gg, gor, wor, _J, L) in enumerate(WELLS):
        res, wh = 10 * (i + 1), 10 * (i + 1) + 1
        inlets.append(PressureBoundary(node_id=res, pressure_pa=P_r * 1e6))
        thermal.append(ThermalBoundary(node_id=res, temperature_c=T, bc_type="fixed_temperature"))
        fluids.append(InletFluidBC(node_id=res, api_gravity=api, gas_gravity=gg,
                                   gor_sc_m3_per_m3=gor, wor_sc_m3_per_m3=wor))
        components.append(make_well(i, tubing_id, correlation))
        components.append(Pipe(component_id=f"flowline_{name.replace(' ', '')}", start_node=wh,
                               end_node=MANIFOLD, diameter_m=0.1016, length_m=L,
                               absolute_roughness_m=46e-6, heat_transfer_coefficient_w_per_m2k=5.0,
                               ambient_temperature_c=15.0))
    components.append(Pipe(component_id="trunk", start_node=MANIFOLD, end_node=SEPARATOR,
                           diameter_m=0.1524, length_m=8000.0, absolute_roughness_m=46e-6,
                           heat_transfer_coefficient_w_per_m2k=3.0, ambient_temperature_c=15.0))
    return NetworkCase(
        name=f"Wells with tubing ({tubing_id * 1000:.1f} mm ID)",
        fluid_model=well_fluid(*WELLS[1][3:7]),
        pressure_inlets=tuple(inlets),
        pressure_outlets=(PressureBoundary(node_id=SEPARATOR, pressure_pa=p_sep),),
        thermal_inlets=tuple(thermal),
        inlet_fluid_bcs=tuple(fluids),
        components=tuple(components),
    )


def well_table(result, tubing_id: float) -> float:
    """Print each well's operating point; return the total oil rate (m³/d)."""
    P = result.node_pressures_pa
    T_node = result.node_temperatures_c
    flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
    states = {ws.component_id: ws for ws in result.well_states}   # why a well is dead
    print(f"  {'Well':7} {'P_wh':>7} {'P_wf':>7} {'ṁ':>7} {'q_oil':>7} {'T_wh':>6}  State")
    print(f"  {'':7} {'(MPa)':>7} {'(MPa)':>7} {'(kg/s)':>7} {'(m³/d)':>7} {'(°C)':>6}")
    total = 0.0
    for i, (name, P_r, T, api, gg, gor, wor, J, _L) in enumerate(WELLS):
        fluid = well_fluid(api, gg, gor, wor)
        curve = well_curve(make_well(i, tubing_id), fluid, P_r * 1e6, T)
        m = flows[f"well_{name.replace(' ', '')}"]
        q = m * 86_400.0 / curve.mass_per_oil_sc_m3
        if q > 1e-6:
            p_b = fluid.pvt(P_r * 1e6, T).bubble_point_pa
            p_wf = f"{bottomhole_pressure(q, P_r * 1e6, J, 'vogel', p_b) / 1e6:7.2f}"
            state = "producing"
            t_wh = f"{T_node[10 * (i + 1) + 1]:6.1f}"
        else:
            q, m, p_wf, t_wh = 0.0, 0.0, "      —", "     —"
            state = states[f"well_{name.replace(' ', '')}"].message
        print(f"  {name:7} {P[10 * (i + 1) + 1] / 1e6:>7.3f} {p_wf} {m:>7.3f} {q:>7.1f} {t_wh}  {state}")
        total += q
    print(f"  Manifold {P[MANIFOLD] / 1e6:.3f} MPa  ·  oil {total:.1f} m³/d  ·  converged {result.converged}")
    return total


def run():
    solver = SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation())

    print("=" * 72)
    print("Tutorial 12 — Wells with production tubing")
    print("=" * 72)

    # ── Part A: deliverability curve of Well 1 ────────────────────────────────
    name, P_r, T, api, gg, gor, wor, J, _L = WELLS[0]
    fluid = well_fluid(api, gg, gor, wor)
    curve = well_curve(make_well(0, 0.0762), fluid, P_r * 1e6, T)
    p_b = fluid.pvt(P_r * 1e6, T).bubble_point_pa
    print(f"\nPart A — {name} deliverability (3.0 in tubing, P_r = {P_r:.0f} MPa, P_b = {p_b / 1e6:.1f} MPa)")
    print(f"  {'q_oil (m³/d)':>13} {'P_wf (MPa)':>11} {'P_wh (MPa)':>11}")
    q_all, p_all = curve.oil_rates_sc_m3_per_day, curve.wellhead_pressures_pa
    for k in sorted({0, curve.peak_index // 2, curve.peak_index, *range(curve.peak_index, len(q_all), 6), len(q_all) - 1}):
        mark = "  ← maximum P_wh (curve turns over)" if k == curve.peak_index else ""
        print(f"  {q_all[k]:>13.0f} {bottomhole_pressure(q_all[k], P_r * 1e6, J, 'vogel', p_b) / 1e6:>11.2f} "
              f"{p_all[k] / 1e6:>11.2f}{mark}")

    # ── Part B: the network ───────────────────────────────────────────────────
    print("\nPart B — gathering network, 3.0 in tubing")
    base = solver.solve(build_case(0.0762))
    well_table(base, 0.0762)

    # ── Part C: tubing size ───────────────────────────────────────────────────
    print("\nPart C — tubing size (all wells)")
    print(f"  {'Tubing ID':>10} {'oil (m³/d)':>11}")
    for tid in (0.0503, 0.0620, 0.0762):
        res = solver.solve(build_case(tid))
        q_tot = 0.0
        for i, (nm, P_r, T, api, gg, gor, wor, _J, _L) in enumerate(WELLS):
            m = next(cf.mass_flow_kg_per_s for cf in res.component_flows
                     if cf.label.endswith(f"well_{nm.replace(' ', '')}"))
            q_tot += max(m, 0.0) * 86_400.0 / well_curve(make_well(i, tid), well_fluid(api, gg, gor, wor),
                                                        P_r * 1e6, T).mass_per_oil_sc_m3
        print(f"  {tid * 1000:>7.1f} mm {q_tot:>11.1f}   converged {res.converged}")

    # ── Part D: tubing correlation ────────────────────────────────────────────
    print("\nPart D — tubing correlation (3.0 in tubing)")
    print(f"  {'Correlation':16} {'Well 1 max P_wh':>16} {'oil (m³/d)':>11}")
    for corr in ("beggs_brill", "hagedorn_brown", "gray"):
        res = solver.solve(build_case(0.0762, corr))
        q_tot = 0.0
        for i, (nm, P_r, T, api, gg, gor, wor, _J, _L) in enumerate(WELLS):
            m = next(cf.mass_flow_kg_per_s for cf in res.component_flows
                     if cf.label.endswith(f"well_{nm.replace(' ', '')}"))
            q_tot += max(m, 0.0) * 86_400.0 / well_curve(make_well(i, 0.0762, corr), well_fluid(api, gg, gor, wor),
                                                        P_r * 1e6, T).mass_per_oil_sc_m3
        c1 = well_curve(make_well(0, 0.0762, corr), well_fluid(*WELLS[0][3:7]), WELLS[0][1] * 1e6, WELLS[0][2])
        print(f"  {corr:16} {c1.max_wellhead_pressure_pa / 1e6:>12.2f} MPa {q_tot:>11.1f}")

    # ── Part E: tubing temperature ────────────────────────────────────────────
    print("\nPart E — tubing temperature of Well 1 at its operating rate")
    m1 = next(cf.mass_flow_kg_per_s for cf in base.component_flows if cf.label.endswith("well_Well1"))
    name, P_r, T, api, gg, gor, wor, J, _L = WELLS[0]
    state = well_tubing_state(make_well(0, 0.0762), well_fluid(api, gg, gor, wor), P_r * 1e6, T, m1,
                              keep_profile=True)
    print(f"  {'Height above':>13} {'Pressure':>9} {'Fluid':>7} {'Formation':>10}")
    print(f"  {'bottom (m)':>13} {'(MPa)':>9} {'(°C)':>7} {'(°C)':>10}")
    for s, P_s, T_s, T_e in state.profile[::5]:
        print(f"  {s * TUBING_RISE / TUBING_LENGTH:>13.0f} {P_s / 1e6:>9.2f} {T_s:>7.1f} {T_e:>10.1f}")
    print(f"  Heat given to the formation: {state.formation_heat_w / 1e3:.0f} kW")
    print("  Wellhead temperature against rate (same well):")
    for m in (1.0, 2.0, 4.0, m1):
        r = well_tubing_state(make_well(0, 0.0762), well_fluid(api, gg, gor, wor), P_r * 1e6, T, m)
        print(f"    {m:5.2f} kg/s → {r.wellhead_temperature_c:5.1f} °C")
    geb = base.global_energy_balance
    print(f"  Network energy balance: wells and pipes lose {geb.heat_wall_loss_kw:.0f} kW, "
          f"lifting takes {-geb.potential_energy_kw:.0f} kW, residual {geb.energy_error_kw:+.3f} kW")

    # ── Part F: back-pressure ─────────────────────────────────────────────────
    print("\nPart F — separator pressure (back-pressure on the wells)")
    print(f"  {'Separator':>10} {'Manifold':>9} {'Well 1':>8} {'Well 2':>8}   (oil, m³/d)")
    last = None
    for p_sep in (1.5e6, 4.0e6, 6.0e6):
        res = solver.solve(build_case(0.0762, p_sep=p_sep))
        rates = []
        for i in (0, 1):
            nm, P_r, T, api, gg, gor, wor, _J, _L = WELLS[i]
            m = next(cf.mass_flow_kg_per_s for cf in res.component_flows
                     if cf.label.endswith(f"well_{nm.replace(' ', '')}"))
            rates.append(max(m, 0.0) * 86_400.0 / well_curve(make_well(i, 0.0762), well_fluid(api, gg, gor, wor),
                                                            P_r * 1e6, T).mass_per_oil_sc_m3)
        print(f"  {p_sep / 1e6:>6.1f} MPa {res.node_pressures_pa[MANIFOLD] / 1e6:>5.2f} MPa "
              f"{rates[0]:>8.1f} {rates[1]:>8.1f}   converged {res.converged}")
        last = res
    print("  At 6 MPa:")
    well_table(last, 0.0762)

    print("\nPart G — nodal analysis of Well 1 at its wellhead (angelica.studies)")
    from angelica.studies import nodal_analysis
    case = build_case()
    na = nodal_analysis(case, lambda: SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation()),
                        "well_Well1", n_rates=7)
    print(f"  {'ṁ (kg/s)':>9} {'inflow P_wh':>12} {'outflow P_wh':>13}")
    print(f"  {'':>9} {'(MPa)':>12} {'(MPa)':>13}")
    peak = max(range(len(na.inflow)), key=lambda k: na.inflow[k][1])
    stable = na.inflow[peak:]
    for m, p_out in na.outflow:
        p_in = next((p0 + (p1 - p0) * (m - m0) / (m1 - m0) for (m0, p0), (m1, p1) in zip(stable, stable[1:])
                     if m0 <= m <= m1 and m1 > m0), float("nan"))
        print(f"  {m:9.3f} {p_in / 1e6:12.3f} {p_out / 1e6 if p_out is not None else float('nan'):13.3f}")
    m_op, p_op = na.operating_point
    net = SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation()).solve(case)
    ws = next(w for w in net.well_states if w.component_id == "well_Well1")
    print(f"  Crossing: {m_op:.4f} kg/s at {p_op / 1e6:.4f} MPa;  network solve: "
          f"{ws.mass_flow_kg_per_s:.4f} kg/s at {ws.wellhead_pressure_pa / 1e6:.4f} MPa")

    print()
    print("  • Well 3's static oil-and-water column (25 °API, WOR 0.8, 2 000 m) weighs")
    print("    more than its reservoir pressure: it cannot flow without artificial lift.")
    print("  • Larger tubing reduces friction in the high-rate wells; the Well model")
    print("    solves each well's tubing on its own, so the network converges in a")
    print("    few outer iterations even with 2 km two-phase risers.")
    print("  • The reservoir fluid warms slightly at the sandface (liquid expanding at")
    print("    constant enthalpy) and then cools along the tubing; the slower well")
    print("    arrives colder.  Wellhead temperature matters downstream: viscosity,")
    print("    wax and hydrate margins all start from it.")
    print("  • Hagedorn-Brown suits oil wells, Gray gas and gas-condensate wells; the")
    print("    spread between correlations is part of the uncertainty of any well model.")
    print("  • Raising the separator pressure lifts the manifold above Well 2's maximum")
    print("    wellhead pressure: Well 2 dies against the manifold while Well 1, with")
    print("    more reservoir pressure and gas, keeps producing.")
    print("  • A well whose curve peaks at a finite rate (a gassy well like Well 1) cannot")
    print("    flow below that rate.  Between its static wellhead pressure and that peak")
    print("    it is bistable: dead or producing.  If producing would lift the manifold")
    print("    above the peak, the well stays dead — shut in, it cannot start against")
    print("    the manifold (it would need unloading or gas lift).  result.well_states")
    print("    gives each well's state and the reason.")
    print("  • The nodal analysis cuts the network at Well 1's wellhead: the well's own")
    print("    curve (inflow) against the pressure the rest of the network needs to take")
    print("    each rate (outflow).  They cross where the network solve puts the well.")


if __name__ == "__main__":
    run()
