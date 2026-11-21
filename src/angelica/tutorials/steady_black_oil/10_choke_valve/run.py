"""
Tutorial 10 — Critical two-phase flow through a production choke
=================================================================

A choke (bean) is the primary flow-control device on a production well.
When the pressure ratio across it falls below a critical value, the
gas-liquid mixture reaches sonic velocity in the throat and the flow becomes
"choked": the rate no longer depends on what happens downstream.  That is
precisely why chokes are used — they isolate the well from separator and
flowline pressure fluctuations.

Network
-------

  Wellhead (N=1)  P=5 MPa  T=60°C   (below bubble point: free gas present)
        │
      Choke (D_pipe = 0.1016 m, bean size varies, Cd = 0.75)
        │
  Downstream of choke (N=2)
        │──── Flowline (4 km, D = 0.1016 m, Beggs-Brill) ────│
                                                    Separator (N=3)

Physics
-------
The choke uses the Sachdeva et al. (1986) model (``Choke(flow_model=
"sachdeva")``, the default).  At the wellhead the fluid carries ~7 % gas
by mass.  The critical pressure ratio is y_c ≈ 0.55; whenever P2/P1 < y_c
the flow is critical.

Part A — bean-size sweep (separator at 0.5 MPa)
  Small beans: P2/P1 ≪ y_c → critical; rate set by the bean alone.
  Large beans: flowline back-pressure pushes P2/P1 above y_c → subcritical.
  The incompressible orifice equation (``flow_model="orifice"``) is shown
  for comparison: it uses the average-pressure mixture density and has no
  sonic limit.

Part B — separator-pressure sweep (24/64" bean)
  The rate stays constant while the choke is critical, then falls once
  P2/P1 exceeds y_c — the classic choke performance curve.

Reference: Sachdeva, R., Schmidt, Z., Brill, J.P. and Blais, R.M. (1986).
Two-Phase Flow Through Chokes. SPE 15657.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures.beggs_brill import BeggsBrillCorrelation
from angelica.closures.choke import black_oil_choke_inputs, choke_flow
from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Choke, Pipe
from angelica.properties.black_oil import BlackOilFluid
from angelica.solvers import SteadyBlackOilSolver

FLUID = BlackOilFluid(
    api_gravity=35.0,
    gas_gravity=0.65,
    gor_sc_m3_per_m3=120.0,
    wor_sc_m3_per_m3=0.2,
    reference_pressure_pa=3e6,
    reference_temperature_c=60.0,
)

P_WELLHEAD = 5.0e6    # Pa
T_WELLHEAD = 60.0     # °C
D_PIPE     = 0.1016   # m (4 inch)
L_FLOWLINE = 4_000.0  # m
CD         = 0.75     # Sachdeva et al. recommendation for liquid-dominated flow


def build_case(orifice_diameter_m: float, p_separator_pa: float = 0.5e6,
               flow_model: str = "sachdeva") -> NetworkCase:
    """Well → choke → flowline → separator for a given bean and model."""
    return NetworkCase(
        name=f"Choke D={orifice_diameter_m*1000:.1f}mm ({flow_model})",
        fluid_model=FLUID,
        pressure_inlets= (PressureBoundary(node_id=1, pressure_pa=P_WELLHEAD),),
        pressure_outlets=(PressureBoundary(node_id=3, pressure_pa=p_separator_pa),),
        thermal_inlets  =(ThermalBoundary(node_id=1, temperature_c=T_WELLHEAD,
                                          bc_type="fixed_temperature"),),
        components=(
            Choke(component_id="wellhead_choke", start_node=1, end_node=2,
                  diameter_m=D_PIPE, orifice_diameter_m=orifice_diameter_m,
                  discharge_coefficient=CD, flow_model=flow_model),
            Pipe(component_id="flowline", start_node=2, end_node=3,
                 diameter_m=D_PIPE, length_m=L_FLOWLINE, absolute_roughness_m=46e-6,
                 heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=20.0),
        ),
    )


def choke_state(result, orifice_diameter_m: float):
    """Re-evaluate the Sachdeva operating point at the solved pressures."""
    P1 = result.node_pressures_pa[1]
    P2 = result.node_pressures_pa[2]
    x1, rho_g, rho_l, k, n = black_oil_choke_inputs(FLUID, P1, T_WELLHEAD)
    area = math.pi * orifice_diameter_m ** 2 / 4.0
    return choke_flow(P1, P2, area, CD, x1, rho_g, rho_l, k, n)


def bean_to_m(bean_64ths: int) -> float:
    return bean_64ths * 0.0254 / 64.0


def run():
    solver = SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation())

    print("=" * 78)
    print("Tutorial 10 — Critical two-phase flow through a production choke")
    print("=" * 78)
    pvt = FLUID.pvt(P_WELLHEAD, T_WELLHEAD)
    x1, _, _, k, n = black_oil_choke_inputs(FLUID, P_WELLHEAD, T_WELLHEAD)
    print(f"\nWellhead: P = {P_WELLHEAD/1e6:.1f} MPa, T = {T_WELLHEAD:.0f} °C, "
          f"Pb = {pvt.bubble_point_pa/1e6:.1f} MPa")
    print(f"Gas mass fraction x1 = {x1:.3f},  k = Cp/Cv = {k:.3f},  polytropic n = {n:.3f}")
    print(f"Flowline: D = {D_PIPE*1000:.0f} mm, L = {L_FLOWLINE/1000:.0f} km (Beggs-Brill)")

    # ── Part A: bean sweep ────────────────────────────────────────────────────
    print("\nPart A — bean-size sweep, separator at 0.5 MPa")
    print(f"{'Bean':>7}  {'ṁ Sachdeva':>11}  {'ṁ orifice':>10}  {'P2 (MPa)':>9}  "
          f"{'P2/P1':>6}  {'y_c':>6}  {'Regime':>12}")
    print(f"{'':>7}  {'(kg/s)':>11}  {'(kg/s)':>10}")
    print("-" * 78)
    for b in (8, 12, 16, 20, 24, 32, 40, 48):
        d = bean_to_m(b)
        res_s = solver.solve(build_case(d))
        res_o = solver.solve(build_case(d, flow_model="orifice"))
        op = choke_state(res_s, d)
        m_s = abs(res_s.component_flows[0].mass_flow_kg_per_s)
        m_o = abs(res_o.component_flows[0].mass_flow_kg_per_s)
        regime = "critical" if op.critical else "subcritical"
        conv = "✓" if (res_s.converged and res_o.converged) else "✗"
        print(f"{b:>4}/64  {m_s:>11.3f}  {m_o:>10.3f}  {res_s.node_pressures_pa[2]/1e6:>9.3f}  "
              f"{op.pressure_ratio:>6.3f}  {op.critical_ratio:>6.3f}  {regime:>12}  {conv}")

    # ── Part B: separator-pressure sweep ─────────────────────────────────────
    d = bean_to_m(24)
    print("\nPart B — separator-pressure sweep, 24/64\" bean")
    print(f"{'P_sep (MPa)':>11}  {'ṁ (kg/s)':>9}  {'P2 (MPa)':>9}  {'P2/P1':>6}  {'Regime':>12}")
    print("-" * 56)
    for p_sep in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5):
        res = solver.solve(build_case(d, p_separator_pa=p_sep * 1e6))
        op = choke_state(res, d)
        regime = "critical" if op.critical else "subcritical"
        print(f"{p_sep:>11.1f}  {abs(res.component_flows[0].mass_flow_kg_per_s):>9.3f}  "
              f"{res.node_pressures_pa[2]/1e6:>9.3f}  {op.pressure_ratio:>6.3f}  {regime:>12}")

    print()
    print("  • While critical, the rate is fixed by the bean: separator changes do")
    print("    not reach the well.  Once P2/P1 > y_c the flow becomes subcritical.")
    print("  • The incompressible orifice equation has no sonic limit and uses an")
    print("    average-pressure density: here it overestimates the rate by ~8 % at")
    print("    48/64\" and by ~40 % for the smallest, deeply critical beans.")


if __name__ == "__main__":
    run()
