"""
Tutorial 11 — Wells with inflow performance (IPR) and check valves
==================================================================

Tutorial 09 fixed each well's pressure at the wellhead, so the weakest well
took back-flow from the manifold.  Here each well is a reservoir that
delivers oil according to its inflow performance relationship (IPR), and a
check valve at each wellhead stops back-flow.

Network
-------

  Reservoir 1 ─IPR─ WH1 ─check─ ─ flowline 4 km, 4"  ─┐
  Reservoir 2 ─IPR─ WH2 ─check─ ─ flowline 2 km, 4"  ─┼─ Manifold ─ trunk 8 km, 6" ─ Separator
  Reservoir 3 ─IPR─ WH3 ─check─ ─ flowline 6 km, 3"  ─┘                              1.5 MPa

  Each reservoir node is a pressure boundary at the well's shut-in pressure
  referred to the wellhead (the tubing is folded into the deliverability
  curve, as is common in gathering-network models), with its own crude:

    Well 1: 40 °API, GOR 150, WOR 0.1, P = 8.0 MPa, 80 °C, J = 25 m³/d/bar
    Well 2: 32 °API, GOR  40, WOR 0.3, P = 6.5 MPa, 65 °C, J = 20 m³/d/bar
    Well 3: 25 °API, GOR  15, WOR 0.8, P = 3.2 MPa, 50 °C, J = 15 m³/d/bar

IPR (``WellInflow``)
--------------------
  linear:  q_o = J·(P_r − P_wf)
  Vogel:   below the bubble point P_b the composite Vogel curve
           q_o = J·(P_r − P_b) + J·P_b/1.8·(1 − 0.2·x − 0.8·x²),  x = P_wf/P_b
           (liberated gas reduces the inflow; all three crudes are saturated
           at these pressures, so P_b is capped at P_r)

  The mass rate adds the produced gas and water,
  ṁ = q_o·(ρ_o,sc + GOR·ρ_g,sc + WOR·ρ_w,sc).

What this tutorial shows
------------------------
  Part A — no check valves: the manifold sits above Well 3's reservoir
           pressure, so oil is injected into Well 3.
  Part B — check valves: Well 3 is shut in, the other two wells carry the
           whole production.
  Part C — Vogel vs linear IPR for the same productivity index.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures import BeggsBrillCorrelation
from angelica.closures.well_inflow import mass_per_oil_sc_m3, oil_rate_sc_m3_per_day
from angelica.core.case import InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import CheckValve, Pipe, WellInflow
from angelica.properties.black_oil import BlackOilFluid
from angelica.solvers import SteadyBlackOilSolver

# name, P_r (MPa), T (°C), API, gas gravity, GOR, WOR, J (m³/d/bar), flowline L (m), D (m)
WELLS = (
    ("Well 1", 8.0, 80.0, 40.0, 0.70, 150.0, 0.1, 25.0, 4000.0, 0.1016),
    ("Well 2", 6.5, 65.0, 32.0, 0.65,  40.0, 0.3, 20.0, 2000.0, 0.1016),
    ("Well 3", 3.2, 50.0, 25.0, 0.60,  15.0, 0.8, 15.0, 6000.0, 0.0762),
)
P_SEP = 1.5e6
MANIFOLD, SEPARATOR = 100, 101


def _nodes(i: int) -> tuple[int, int, int]:
    """(reservoir, wellhead, downstream of check valve) node ids of well i."""
    return 10 * (i + 1), 10 * (i + 1) + 1, 10 * (i + 1) + 2


def well_fluid(api, gg, gor, wor) -> BlackOilFluid:
    return BlackOilFluid(api_gravity=api, gas_gravity=gg, gor_sc_m3_per_m3=gor,
                         wor_sc_m3_per_m3=wor, reference_pressure_pa=5e6,
                         reference_temperature_c=65.0)


def build_case(check_valves: bool = True, ipr_model: str = "vogel") -> NetworkCase:
    components, inlets, thermal, fluids = [], [], [], []
    for i, (name, P_r, T, api, gg, gor, wor, J, L, D) in enumerate(WELLS):
        res, wh, down = _nodes(i)
        tag = name.replace(" ", "")
        inlets.append(PressureBoundary(node_id=res, pressure_pa=P_r * 1e6))
        thermal.append(ThermalBoundary(node_id=res, temperature_c=T, bc_type="fixed_temperature"))
        fluids.append(InletFluidBC(node_id=res, api_gravity=api, gas_gravity=gg,
                                   gor_sc_m3_per_m3=gor, wor_sc_m3_per_m3=wor))
        components.append(WellInflow(component_id=f"ipr_{tag}", start_node=res, end_node=wh,
                                     diameter_m=0.1016, productivity_index_m3_per_day_per_bar=J,
                                     ipr_model=ipr_model))
        start = wh
        if check_valves:
            components.append(CheckValve(component_id=f"check_{tag}", start_node=wh, end_node=down,
                                         diameter_m=D))
            start = down
        components.append(Pipe(component_id=f"flowline_{tag}", start_node=start, end_node=MANIFOLD,
                               diameter_m=D, length_m=L, absolute_roughness_m=46e-6,
                               heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0))
    components.append(Pipe(component_id="trunk", start_node=MANIFOLD, end_node=SEPARATOR,
                           diameter_m=0.1524, length_m=8000.0, absolute_roughness_m=46e-6,
                           heat_transfer_coefficient_w_per_m2k=3.0, ambient_temperature_c=15.0))
    return NetworkCase(
        name=f"Wells with IPR ({ipr_model}){' and check valves' if check_valves else ''}",
        fluid_model=well_fluid(*WELLS[1][3:7]),
        pressure_inlets=tuple(inlets),
        pressure_outlets=(PressureBoundary(node_id=SEPARATOR, pressure_pa=P_SEP),),
        thermal_inlets=tuple(thermal),
        inlet_fluid_bcs=tuple(fluids),
        components=tuple(components),
    )


def well_table(result) -> float:
    """Print the operating point of each well; return the total oil rate."""
    P = result.node_pressures_pa
    flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
    print(f"  {'Well':7} {'P_r':>6} {'P_wh':>7} {'ṁ':>8} {'q_oil':>9}  State")
    print(f"  {'':7} {'(MPa)':>6} {'(MPa)':>7} {'(kg/s)':>8} {'(m³/d)':>9}")
    total = 0.0
    for i, (name, P_r, _T, api, gg, gor, wor, *_rest) in enumerate(WELLS):
        _res, wh, _down = _nodes(i)
        m = flows[f"ipr_{name.replace(' ', '')}"]
        m = 0.0 if abs(m) < 1e-6 else m          # shut-in well: round-off level
        q = m * 86_400.0 / mass_per_oil_sc_m3(well_fluid(api, gg, gor, wor))
        state = "producing" if m > 1e-6 else ("INJECTED (back-flow)" if m < -1e-6 else "shut in by check valve")
        print(f"  {name:7} {P_r:>6.2f} {P[wh] / 1e6:>7.3f} {m:>8.3f} {q:>9.1f}  {state}")
        total += max(q, 0.0)
    print(f"  Manifold {P[MANIFOLD] / 1e6:.3f} MPa  ·  oil to separator {total:.1f} m³/d  ·  "
          f"converged {result.converged}")
    return total


def run():
    solver = SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation())

    print("=" * 74)
    print("Tutorial 11 — Wells with IPR and check valves (Beggs-Brill flowlines)")
    print("=" * 74)

    print("\nIPR of Well 1 (Vogel vs linear, J = 25 m³/d/bar, P_r = 8 MPa):")
    _n, P_r, T, api, gg, gor, wor, J, *_ = WELLS[0]
    p_b = well_fluid(api, gg, gor, wor).pvt(P_r * 1e6, T).bubble_point_pa
    print(f"  bubble point at {T:.0f} °C: {p_b / 1e6:.1f} MPa (> P_r: saturated reservoir)")
    print(f"  {'P_wf (MPa)':>11} {'Vogel (m³/d)':>13} {'linear (m³/d)':>14}")
    for p_wf in (8.0, 6.0, 4.0, 2.0, 0.0):
        qv = oil_rate_sc_m3_per_day(P_r * 1e6, p_wf * 1e6, J, "vogel", p_b)
        ql = oil_rate_sc_m3_per_day(P_r * 1e6, p_wf * 1e6, J, "linear", p_b)
        print(f"  {p_wf:>11.1f} {qv:>13.1f} {ql:>14.1f}")

    print("\nPart A — no check valves (Vogel IPR)")
    well_table(solver.solve(build_case(check_valves=False)))

    print("\nPart B — check valve at every wellhead (Vogel IPR)")
    q_vogel = well_table(solver.solve(build_case(check_valves=True)))

    print("\nPart C — same network with a linear IPR")
    q_linear = well_table(solver.solve(build_case(check_valves=True, ipr_model="linear")))

    print()
    print("  • Without check valves the manifold pressure exceeds Well 3's reservoir")
    print("    pressure and oil flows back into the reservoir.  The check valve shuts")
    print("    the well in instead; its wellhead then sits at the reservoir pressure.")
    print("  • The linear IPR ignores the gas liberated near the wellbore and")
    print(f"    over-predicts production: {q_linear:.0f} vs {q_vogel:.0f} m³/d of oil "
          f"({100 * (q_linear / q_vogel - 1):+.0f} %).")


if __name__ == "__main__":
    run()
