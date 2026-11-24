"""
Tutorial 13 — Gas lift
======================

Well 3 of tutorial 12 — a heavy, wet crude (25 °API, GOR 15, WOR 0.8) — is
dead: its static column of oil and water outweighs the reservoir pressure.
Gas lift injects gas into the tubing through a valve; above the valve the
column is lighter and the well flows.

  Reservoir ─ Well (IPR + 2 100 m tubing, gas injected at depth) ─ WH ─ flowline 2 km, 4" ─ Separator
  17 MPa, 75 °C                                                                              1.5 MPa

How gas lift is modelled (``Well(gas_lift_rate_sm3_per_day=…, gas_lift_depth_m=…)``)
-------------------------------------------------------------------------------------
  - Below the valve the tubing carries the reservoir fluid; above it, the
    reservoir fluid plus the lift gas.  The lift gas reaches the valve at the
    formation temperature of that depth and stays gas (no mass transfer with
    the oil).
  - The lift gas leaves the well at the wellhead node, where the network
    receives it as an internal source: the flowline carries production +
    lift gas, and its gas-oil ratio includes the lift gas.

What to look at
---------------
  A. The gas-lift performance curve: oil rate against injection rate rises
     to a maximum and then falls, when the friction of the extra gas costs
     more than its lightening of the column.
  B. Deeper injection lightens more of the column.
  C. The tubing profile with the injection point.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures import BeggsBrillCorrelation
from angelica.closures.well_model import well_curve, well_tubing_state
from angelica.core.case import InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe, Well
from angelica.properties.black_oil import BlackOilFluid
from angelica.solvers import SteadyBlackOilSolver

CRUDE = {"api_gravity": 25.0, "gas_gravity": 0.60, "gor_sc_m3_per_m3": 15.0, "wor_sc_m3_per_m3": 0.8}
P_RES, T_RES, J = 17.0e6, 75.0, 3.0
TUBING_LENGTH, TUBING_RISE = 2100.0, 2000.0
P_SEP = 1.5e6


def make_well(q_gl: float = 0.0, depth: float = 1800.0) -> Well:
    lift = {"gas_lift_rate_sm3_per_day": q_gl, "gas_lift_depth_m": depth} if q_gl > 0.0 else {}
    return Well(component_id="well", start_node=1, end_node=2, diameter_m=0.0762,
                productivity_index_m3_per_day_per_bar=J, ipr_model="vogel",
                tubing_length_m=TUBING_LENGTH, tubing_height_change_m=TUBING_RISE, **lift)


def build_case(q_gl: float = 0.0, depth: float = 1800.0) -> NetworkCase:
    return NetworkCase(
        name=f"Gas lift {q_gl:.0f} Sm3/d at {depth:.0f} m",
        fluid_model=BlackOilFluid(**CRUDE),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=P_RES),),
        pressure_outlets=(PressureBoundary(node_id=3, pressure_pa=P_SEP),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=T_RES, bc_type="fixed_temperature"),),
        inlet_fluid_bcs=(InletFluidBC(node_id=1, **CRUDE),),
        components=(make_well(q_gl, depth),
                    Pipe(component_id="flowline", start_node=2, end_node=3, diameter_m=0.1016, length_m=2000.0,
                         absolute_roughness_m=46e-6, heat_transfer_coefficient_w_per_m2k=5.0,
                         ambient_temperature_c=15.0)),
    )


def oil_rate(result, q_gl: float, depth: float) -> float:
    curve = well_curve(make_well(q_gl, depth), BlackOilFluid(**CRUDE), P_RES, T_RES)
    return result.component_flows[0].mass_flow_kg_per_s * 86_400.0 / curve.mass_per_oil_sc_m3


def run():
    solver = SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation())
    print("=" * 74)
    print("Tutorial 13 — Gas lift (heavy, wet crude; Beggs-Brill tubing and flowline)")
    print("=" * 74)

    print("\nPart A — gas-lift performance curve (valve at 1 800 m)")
    print(f"  {'Lift gas':>10} {'Oil':>8} {'P_wh':>7} {'T_wh':>6} {'Valve P':>8}  State")
    print(f"  {'(Sm³/d)':>10} {'(m³/d)':>8} {'(MPa)':>7} {'(°C)':>6} {'(MPa)':>8}")
    best = (0.0, 0.0)
    for q_gl in (0.0, 5_000.0, 10_000.0, 20_000.0, 40_000.0, 70_000.0, 100_000.0, 150_000.0, 250_000.0):
        res = solver.solve(build_case(q_gl))
        q_o = oil_rate(res, q_gl, 1800.0)
        ws = res.well_states[0]
        if q_o > best[1]:
            best = (q_gl, q_o)
        p_valve = "      —"
        if ws.state == "producing":
            st = well_tubing_state(make_well(q_gl), BlackOilFluid(**CRUDE), P_RES, T_RES, ws.mass_flow_kg_per_s)
            if q_gl > 0.0:
                p_valve = f"{st.injection_pressure_pa / 1e6:>8.2f}"
        t_wh = f"{res.node_temperatures_c[2]:6.1f}" if ws.state == "producing" else "     —"
        print(f"  {q_gl:>10.0f} {q_o:>8.1f} {res.node_pressures_pa[2] / 1e6:>7.3f} {t_wh} {p_valve}  {ws.message}")
    print(f"  Maximum oil {best[1]:.0f} m³/d at {best[0]:.0f} Sm³/d of lift gas")

    print("\nPart B — injection depth (60 000 Sm³/d)")
    print(f"  {'Valve depth':>12} {'Oil (m³/d)':>11} {'Valve P (MPa)':>14}")
    for depth in (600.0, 1200.0, 1800.0, 2100.0):
        res = solver.solve(build_case(60_000.0, depth))
        st = well_tubing_state(make_well(60_000.0, depth), BlackOilFluid(**CRUDE), P_RES, T_RES,
                               res.component_flows[0].mass_flow_kg_per_s)
        print(f"  {depth:>10.0f} m {oil_rate(res, 60_000.0, depth):>11.1f} {st.injection_pressure_pa / 1e6:>14.2f}")

    print("\nPart C — tubing profile at 60 000 Sm³/d, valve at 1 800 m")
    res = solver.solve(build_case(60_000.0))
    m = res.component_flows[0].mass_flow_kg_per_s
    st = well_tubing_state(make_well(60_000.0), BlackOilFluid(**CRUDE), P_RES, T_RES, m, keep_profile=True)
    s_valve = TUBING_LENGTH - 1800.0
    print(f"  {'Depth (m)':>10} {'P (MPa)':>8} {'T (°C)':>7} {'Formation':>10}")
    for s_, P_s, T_s, T_e in st.profile:
        depth = TUBING_LENGTH - s_
        if abs(s_ - s_valve) < 1.0 or round(depth) % 300 == 0:
            mark = "  ← gas injected" if abs(s_ - s_valve) < 1.0 else ""
            print(f"  {depth:>10.0f} {P_s / 1e6:>8.2f} {T_s:>7.1f} {T_e:>10.1f}{mark}")
    m_line = res.component_flows[1].mass_flow_kg_per_s
    geb = res.global_energy_balance
    print(f"  Reservoir {m:.3f} kg/s + lift gas {res.well_states[0].gas_lift_kg_per_s:.3f} kg/s = "
          f"flowline {m_line:.3f} kg/s")
    print(f"  Energy balance residual {geb.energy_error_kw:+.3f} kW")

    print()
    print("  • Without lift the well is dead; a little gas brings it to life.")
    print("  • Past the optimum the extra gas adds more friction than it removes")
    print("    weight: the oil rate falls although more gas is injected.")
    print("  • A deeper valve lightens a longer part of the column, but needs a")
    print("    higher injection pressure at the valve (Part B).")


if __name__ == "__main__":
    run()
