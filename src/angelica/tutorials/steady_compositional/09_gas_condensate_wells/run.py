"""
Tutorial 09 — Gas and gas-condensate wells in a gathering network
=================================================================

Two wells produce into a manifold and a trunk line to a separator:

  Reservoir G ─ gas well (3 000 m tubing) ─ WH-G ─ flowline 3 km, 6" ─┐
                                                                      ├─ Manifold ─ trunk 10 km, 8" ─ Separator
  Reservoir C ─ condensate well (3 200 m) ─ WH-C ─ flowline 4 km, 6" ─┘                              60 bar

    Gas well:        lean gas (86 % methane),  250 bar, 100 °C
    Condensate well: rich gas (74 % methane, 4 % n-decane), 280 bar, 110 °C
    Both: gas back-pressure IPR  q_g = C·(P_r² − P_wf²)ⁿ  (bar, standard m³/d),
          3-in. tubing, Gray correlation (gas and gas-condensate wells).

The ``Well`` component works as in black-oil tutorial 12: IPR + tubing give
the deliverability curve P_wh(q); the network uses its stable branch; the
fluid enters the tubing at its reservoir enthalpy and loses heat to the
formation (Ramey).  Here the properties come from the Peng-Robinson
equation of state: each well's fluid is flashed once on a pressure–
temperature table and interpolated along the tubing.

What to look at
---------------
  A. Deliverability curves.  The standard volumes come from a flash of the
     stream to standard conditions: the condensate well also gives a stock-
     tank liquid rate.
  B. The network: each well's rate, wellhead temperature and state.  The gas
     cools on the drawdown (Joule-Thomson) instead of warming like a liquid.
  C. The condensate well's tubing: pressure, temperature and liquid content.
     The liquid mass fraction grows as the fluid climbs and cools
     (condensation), while the liquid volume fraction falls (the gas
     expands faster).
  D. Separator pressure: the condensate well dies against the manifold
     first.  Its column, loaded with liquid, is heavier: its shut-in
     wellhead pressure is lower than the gas well's although its reservoir
     pressure is higher.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures.well_fluid import well_fluid_for
from angelica.closures.well_model import well_curve, well_tubing_state
from angelica.core.case import InletCompositionBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe, Well
from angelica.properties.compositional_fluid import CompositionalFluid, _hc_phase_split
from angelica.solvers import SteadyCompositionalSolver

COMPONENTS = ("methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane", "n-decane",
              "nitrogen", "carbon dioxide")
# name, mole fractions, P_r (bar), T_r (°C), C, n, tubing length = rise (m), flowline (m)
WELLS = (
    ("Gas", (0.86, 0.06, 0.03, 0.015, 0.005, 0.005, 0.005, 0.01, 0.01), 250.0, 100.0, 150.0, 0.80, 3000.0, 3000.0),
    ("Condensate", (0.74, 0.08, 0.05, 0.03, 0.02, 0.02, 0.04, 0.01, 0.01), 280.0, 110.0, 120.0, 0.75, 3200.0, 4000.0),
)
BAR = 1e5
P_SEP = 60.0
MANIFOLD, SEPARATOR = 100, 101


def make_well(i: int) -> Well:
    name, _zs, _p, _t, C, n, L, _fl = WELLS[i]
    return Well(component_id=f"well_{name}", start_node=10 * (i + 1), end_node=10 * (i + 1) + 1,
                diameter_m=0.0762, ipr_model="backpressure", ipr_coefficient=C, ipr_exponent=n,
                tubing_length_m=L, tubing_height_change_m=L, tubing_correlation="gray")


def fluid_of(i: int) -> CompositionalFluid:
    return CompositionalFluid(COMPONENTS, WELLS[i][1])


def build_case(p_sep_bar: float = P_SEP) -> NetworkCase:
    components, inlets, thermal, comps = [], [], [], []
    for i, (name, zs, p_r, t_r, _C, _n, _L, fl) in enumerate(WELLS):
        res, wh = 10 * (i + 1), 10 * (i + 1) + 1
        inlets.append(PressureBoundary(node_id=res, pressure_pa=p_r * BAR))
        thermal.append(ThermalBoundary(node_id=res, temperature_c=t_r, bc_type="fixed_temperature"))
        comps.append(InletCompositionBC(node_id=res, zs=zs))
        components.append(make_well(i))
        components.append(Pipe(component_id=f"flowline_{name}", start_node=wh, end_node=MANIFOLD,
                               diameter_m=0.1524, length_m=fl, absolute_roughness_m=46e-6,
                               heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0))
    components.append(Pipe(component_id="trunk", start_node=MANIFOLD, end_node=SEPARATOR, diameter_m=0.2032,
                           length_m=10_000.0, absolute_roughness_m=46e-6,
                           heat_transfer_coefficient_w_per_m2k=3.0, ambient_temperature_c=15.0))
    return NetworkCase(
        name=f"Gas and condensate wells, separator {p_sep_bar:.0f} bar",
        fluid_model=fluid_of(0),
        pressure_inlets=tuple(inlets),
        pressure_outlets=(PressureBoundary(node_id=SEPARATOR, pressure_pa=p_sep_bar * BAR),),
        thermal_inlets=tuple(thermal),
        inlet_composition_bcs=tuple(comps),
        components=tuple(components),
    )


def stream(i: int):
    """Standard-conditions view of a well's stream (gas and condensate volumes)."""
    return well_fluid_for(fluid_of(i))


def well_table(result) -> None:
    flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
    states = {ws.component_id: ws for ws in result.well_states}
    print(f"  {'Well':11} {'P_wh':>6} {'T_wh':>6} {'ṁ':>6} {'gas':>9} {'condensate':>11}  State")
    print(f"  {'':11} {'(bar)':>6} {'(°C)':>6} {'(kg/s)':>6} {'(kSm³/d)':>9} {'(Sm³/d)':>11}")
    for i, (name, *_rest) in enumerate(WELLS):
        wh = 10 * (i + 1) + 1
        m = flows[f"well_{name}"]
        mass, v_gas, v_oil = stream(i)._standard_volumes
        kmol_per_s = m / mass / 1000.0 if m > 0 else 0.0          # (mass is kg per mole of feed)
        gas = kmol_per_s * 1000.0 * v_gas * 86_400.0
        cond = kmol_per_s * 1000.0 * v_oil * 86_400.0
        t_wh = f"{result.node_temperatures_c[wh]:6.1f}" if m > 0 else "     —"
        print(f"  {name:11} {result.node_pressures_pa[wh] / BAR:>6.1f} {t_wh} {m:>6.2f} {gas / 1e3:>9.0f} "
              f"{cond:>11.1f}  {states[f'well_{name}'].message}")
    geb = result.global_energy_balance
    print(f"  Manifold {result.node_pressures_pa[MANIFOLD] / BAR:.1f} bar  ·  converged {result.converged}  ·  "
          f"energy residual {geb.energy_error_kw:+.3f} kW")


def run():
    solver = SteadyCompositionalSolver()
    print("=" * 76)
    print("Tutorial 09 — Gas and gas-condensate wells (Peng-Robinson, Gray tubing)")
    print("=" * 76)

    print("\nPart A — deliverability curves (rate in standard m³/d of gas)")
    for i, (name, zs, p_r, t_r, C, n, L, _fl) in enumerate(WELLS):
        s = stream(i)
        mass, v_gas, v_oil = s._standard_volumes
        yield_ = v_oil / v_gas * 1e6 if v_gas > 0 else 0.0
        curve = well_curve(make_well(i), fluid_of(i), p_r * BAR, t_r)
        q, p = curve.oil_rates_sc_m3_per_day, curve.wellhead_pressures_pa
        print(f"  {name} well: P_r {p_r:.0f} bar, condensate yield {yield_:.0f} Sm³ per million Sm³ of gas, "
              f"shut-in P_wh {curve.static_wellhead_pressure_pa / BAR:.1f} bar, "
              f"maximum {curve.max_wellhead_pressure_pa / BAR:.1f} bar")
        for k in (0, curve.peak_index, len(q) // 3, 2 * len(q) // 3, len(q) - 1):
            print(f"    {q[k] / 1e3:>8.0f} kSm³/d → P_wh {p[k] / BAR:6.1f} bar")

    print(f"\nPart B — network, separator {P_SEP:.0f} bar")
    base = solver.solve(build_case())
    well_table(base)

    print("\nPart C — condensate well tubing at its operating rate")
    m = next(cf.mass_flow_kg_per_s for cf in base.component_flows if cf.label.endswith("well_Condensate"))
    name, zs, p_r, t_r, *_ = WELLS[1]
    state = well_tubing_state(make_well(1), fluid_of(1), p_r * BAR, t_r, m, keep_profile=True)
    wf = well_fluid_for(fluid_of(1), None, make_well(1), round(p_r * BAR), round(t_r, 3))
    print(f"  {'Height':>8} {'P':>7} {'T':>7} {'Formation':>10} {'Liquid':>7} {'Liquid':>7}")
    print(f"  {'(m)':>8} {'(bar)':>7} {'(°C)':>7} {'(°C)':>10} {'mass':>7} {'volume':>7}")
    for s_, P_s, T_s, T_e in state.profile[::5]:
        C_L = wf.phase_properties(P_s, T_s)[0]
        beta, _rg, _rl, mw_g, mw_l, *_ = _hc_phase_split(COMPONENTS, round(P_s, -1), round(T_s, 3), zs, "PR")
        x_l = (1 - beta) * mw_l / ((1 - beta) * mw_l + beta * mw_g)
        print(f"  {s_:>8.0f} {P_s / BAR:>7.1f} {T_s:>7.1f} {T_e:>10.1f} {100 * x_l:>6.1f}% {100 * C_L:>6.1f}%")
    print(f"  Bottom hole {state.bottomhole_temperature_c:.1f} °C (reservoir {t_r:.0f} °C); "
          f"heat to the formation {state.formation_heat_w / 1e3:.0f} kW")

    print("\nPart D — separator pressure")
    print(f"  {'Separator':>10} {'Gas well':>16} {'Condensate well':>16}   (kg/s)")
    last = None
    for p_sep in (60.0, 150.0, 186.0):
        res = solver.solve(build_case(p_sep))
        st = {ws.component_id: ws for ws in res.well_states}
        print(f"  {p_sep:>6.0f} bar {st['well_Gas'].mass_flow_kg_per_s:>16.2f} "
              f"{st['well_Condensate'].mass_flow_kg_per_s:>16.2f}   converged {res.converged}")
        last = res
    print("  At 186 bar:")
    well_table(last)

    print()
    print("  • The standard gas and condensate rates come from a flash of each stream to")
    print("    standard conditions (Peneloux-corrected liquid volume).")
    print("  • Gas cools as it expands into the well (Joule-Thomson) and keeps cooling")
    print("    up the tubing; liquid condenses on the way.")
    print("  • The condensate well's liquid-laden column is the heavier one: it dies")
    print("    against the manifold first; result.well_states says why.")


if __name__ == "__main__":
    run()
