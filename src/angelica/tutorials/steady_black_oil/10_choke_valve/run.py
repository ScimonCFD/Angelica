"""
Tutorial 10 — Choke valve in a production well flowline
========================================================

A choke (bean) is the primary flow-control device on a production well.
By restricting the cross-section, it throttles the wellhead pressure down
to separator pressure while limiting sand production and protecting surface
equipment.

This tutorial models a well → choke → flowline → separator system and
shows how the bean size controls production rate.

Network
-------

  Wellhead (N=1)  P=20 MPa  T=90°C
        │
      Choke (D_pipe=0.1016 m, various bean sizes)
        │
  Downstream of choke (N=2)
        │──── Flowline (8 km, D=0.1016 m) ────│
                                             Separator (N=3)  P=0.5 MPa

Physics
-------
The Choke component is modeled as a sharp-edged orifice:

    ṁ = Cd · A_orifice · √(2 · ρ · ΔP_choke)

with default Cd = 0.62 (ISO 5167 sharp-edged orifice plate).
Pressure drop across the choke is ΔP = K · ρ · v²/2 where
K = (D_pipe / D_orifice)⁴ / Cd².

The flowline uses the Beggs-Brill correlation for two-phase pressure drop.

Bean-size sweep
---------------
Eight bean sizes from 8/64" to 48/64" are compared.  The results show:

  - Production rate increases strongly with bean size (roughly ∝ D²)
  - The choke pressure drop falls as bean size increases
  - At large bean sizes the flowline friction becomes the dominant
    restriction and gains per additional bean size diminish

Reference: Choke sizing methodology from Guo & Ghalambor (2005),
"Natural Gas Engineering Handbook", 2nd ed.  Orifice coefficient from
ISO 5167-1:2003.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures.beggs_brill import BeggsBrillCorrelation
from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Choke, Pipe
from angelica.properties.black_oil import BlackOilFluid
from angelica.solvers import SteadyBlackOilSolver

# ── Fluid ──────────────────────────────────────────────────────────────────────
FLUID = BlackOilFluid(
    api_gravity=35.0,
    gas_gravity=0.65,
    gor_sc_m3_per_m3=60.0,
    wor_sc_m3_per_m3=0.2,
    reference_pressure_pa=10e6,
    reference_temperature_c=80.0,
)

D_PIPE     = 0.1016   # m  (4 inch pipe)
D_FLOWLINE = 0.1016   # m
L_FLOWLINE = 8_000.0  # m


def build_case(orifice_diameter_m: float) -> NetworkCase:
    """Build the well-choke-flowline-separator case for a given bean size."""
    return NetworkCase(
        name=f"Choke D={orifice_diameter_m*1000:.0f}mm",
        fluid_model=FLUID,
        pressure_inlets= (PressureBoundary(node_id=1, pressure_pa=20e6),),
        pressure_outlets=(PressureBoundary(node_id=3, pressure_pa=0.5e6),),
        thermal_inlets  =(ThermalBoundary(node_id=1, temperature_c=90.0,
                                          bc_type="fixed_temperature"),),
        components=(
            Choke(
                component_id="wellhead_choke",
                start_node=1, end_node=2,
                diameter_m=D_PIPE,
                orifice_diameter_m=orifice_diameter_m,
                discharge_coefficient=0.62,
            ),
            Pipe(
                component_id="flowline",
                start_node=2, end_node=3,
                diameter_m=D_FLOWLINE, length_m=L_FLOWLINE,
                absolute_roughness_m=46e-6,
                heat_transfer_coefficient_w_per_m2k=5.0,
                ambient_temperature_c=20.0,
            ),
        ),
    )


def run():
    # Bean sizes in 64ths of an inch (oilfield convention)
    beans_64ths = [8, 12, 16, 20, 24, 32, 40, 48]

    print("=" * 70)
    print("Tutorial 10 — Choke valve: bean-size sensitivity")
    print("=" * 70)
    print("\nWell P = 20 MPa, T = 90 °C  →  Separator P = 0.5 MPa")
    print(f"Flowline: D={D_FLOWLINE*1000:.0f} mm, L={L_FLOWLINE/1000:.0f} km (Beggs-Brill)")
    print(f"\n{'Bean':>8}  {'D_orifice':>10}  {'K_choke':>9}  "
          f"{'q (m³/h)':>10}  {'ΔP_choke':>10}  {'ΔP_line':>10}  {'P_2 (MPa)':>10}")
    print("-" * 75)

    solver = SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation())

    q_max = 0.0
    for b in beans_64ths:
        D_orifice = b * 0.0254 / 64.0      # convert 64ths of inch → metres
        K = (D_PIPE / D_orifice) ** 4 / 0.62 ** 2

        case   = build_case(D_orifice)
        result = solver.solve(case)

        q   = abs(result.component_flows[0].volumetric_flow_m3_per_h)
        P1  = result.node_pressures_pa[1]
        P2  = result.node_pressures_pa[2]
        P3  = result.node_pressures_pa[3]
        dP_choke = (P1 - P2) / 1e6    # MPa
        dP_line  = (P2 - P3) / 1e6    # MPa

        q_max = max(q_max, q)
        conv  = "✓" if result.converged else "✗"
        print(f"{b:>6}/64  {D_orifice*1000:>7.1f} mm  {K:>9.1f}  "
              f"{q:>10.2f}  {dP_choke:>8.2f} MPa  {dP_line:>8.2f} MPa  "
              f"{P2/1e6:>8.3f} MPa  {conv}")

    print()
    print(f"  Note: flow rate spans from small-bean to large-bean: {q_max:.1f} m³/h max.")
    print("  At large beans the flowline friction dominates → diminishing returns.")
    print()
    print("  Physical insight:")
    print("  • ΔP_choke decreases rapidly as bean size grows (∝ D⁻⁴)")
    print("  • ΔP_flowline increases as flow increases (∝ q²)")
    print("  • The operating point is where both sum to the total available ΔP")


if __name__ == "__main__":
    run()
