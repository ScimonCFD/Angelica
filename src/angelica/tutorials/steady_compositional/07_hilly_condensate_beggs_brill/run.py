"""
Tutorial 07 — Gas-condensate line over hilly terrain: homogeneous vs Beggs-Brill
================================================================================

A wet gas (with condensate and produced water) flows 20 km over hilly
terrain.  The same network is solved twice:

  • homogeneous no-slip flow (the compositional default): gas and liquid
    move at the same velocity, so the liquid in the pipe equals the
    no-slip liquid fraction C_L;
  • Beggs-Brill: the gas slips past the liquid.  Uphill the liquid
    accumulates (holdup H_L > C_L), the mixture column is heavier and the
    two-phase friction factor rises; downhill the liquid drains faster
    and the holdup is lower than uphill.

Network (pressure at both ends, flow is the result)
---------------------------------------------------

  Inlet ── A: 5 km, +80 m ── B: 5 km, −60 m ── C: 5 km, +120 m ── D: 5 km, −40 m ── Outlet
  80 bar, 45 °C                                                                   55 bar

  10 in (254 mm) line, U = 4 W/m²K to soil at 10 °C.

Fluid (mol %): methane 78, ethane 8, propane 5, n-butane 3, n-pentane 2,
n-hexane 2, water 2 — Peng-Robinson; condensate and free water form the
liquid, water vapour travels with the gas.

What to look at
---------------
  1. Beggs-Brill delivers less gas for the same pressures: the line holds
     more liquid, which is heavier to lift and raises the friction factor.
  2. Section table (``result.component_flows[i].multiphase``): flow
     pattern, no-slip C_L, in-situ H_L, superficial gas velocity, and the
     split of the pressure drop between elevation and friction.
  3. Liquid inventory: the line holds far more liquid than the homogeneous
     model predicts — what a pig would push out, or a slug catcher see.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures import BeggsBrillCorrelation
from angelica.core.case import InletCompositionBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.compositional_fluid import CompositionalFluid
from angelica.solvers import SteadyCompositionalSolver

COMPONENTS = ["methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane", "water"]
ZS = [0.78, 0.08, 0.05, 0.03, 0.02, 0.02, 0.02]
P_IN, P_OUT, T_IN, T_AMB = 80e5, 55e5, 45.0, 10.0
D, ROUGH, U = 0.254, 46e-6, 4.0
SECTIONS = (("A", 5000.0, 80.0), ("B", 5000.0, -60.0), ("C", 5000.0, 120.0), ("D", 5000.0, -40.0))


def build_case() -> NetworkCase:
    pipes = tuple(
        Pipe(component_id=f"section_{name}", start_node=i + 1, end_node=i + 2, diameter_m=D,
             length_m=L, height_change_m=dz, absolute_roughness_m=ROUGH,
             heat_transfer_coefficient_w_per_m2k=U, ambient_temperature_c=T_AMB)
        for i, (name, L, dz) in enumerate(SECTIONS)
    )
    return NetworkCase(
        name="Hilly gas-condensate line",
        fluid_model=CompositionalFluid(components=COMPONENTS, default_zs=ZS),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=P_IN),),
        pressure_outlets=(PressureBoundary(node_id=len(SECTIONS) + 1, pressure_pa=P_OUT),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=T_IN, bc_type="fixed_temperature"),),
        inlet_composition_bcs=(InletCompositionBC(node_id=1, zs=tuple(ZS)),),
        components=pipes,
    )


def section_table(res) -> None:
    """Per-section two-phase state from ``ComponentFlowResult.multiphase``."""
    P, T = res.node_pressures_pa, res.node_temperatures_c
    print(f"  {'Sec':3} {'dz':>5} {'P_out':>6} {'T_out':>6} {'pattern':>13} {'C_L':>6} {'H_L':>6}"
          f" {'v_sg':>6} {'ΔP elev':>8} {'ΔP fric':>8} {'liquid':>7}")
    print(f"  {'':3} {'(m)':>5} {'(bar)':>6} {'(°C)':>6} {'':>13} {'':>6} {'':>6}"
          f" {'(m/s)':>6} {'(bar)':>8} {'(bar)':>8} {'(m³)':>7}")
    inv_ns = inv = 0.0
    for i, ((name, L, dz), cf) in enumerate(zip(SECTIONS, res.component_flows)):
        mp = cf.multiphase
        inv += mp.liquid_volume_m3
        inv_ns += mp.no_slip_liquid_fraction * mp.liquid_volume_m3 / mp.liquid_holdup
        print(f"  {name:3} {dz:>+5.0f} {P[i + 2] / 1e5:>6.2f} {T[i + 2]:>6.1f} {mp.flow_pattern:>13}"
              f" {mp.no_slip_liquid_fraction:>6.3f} {mp.liquid_holdup:>6.3f}"
              f" {mp.superficial_gas_velocity_m_per_s:>6.2f}"
              f" {mp.elevation_gradient_pa_per_m * L / 1e5:>8.2f} {mp.friction_gradient_pa_per_m * L / 1e5:>8.2f}"
              f" {mp.liquid_volume_m3:>7.1f}")
    print(f"  Liquid in the line: {inv:.0f} m³ with slip  ·  {inv_ns:.0f} m³ at the no-slip fraction")


def run() -> None:
    print("=" * 76)
    print("Tutorial 07 — Hilly gas-condensate line: homogeneous vs Beggs-Brill")
    print("=" * 76)

    hom = SteadyCompositionalSolver().solve(build_case())
    bb = SteadyCompositionalSolver(turbulent_pipe_correlation=BeggsBrillCorrelation()).solve(build_case())

    m_h = hom.component_flows[0].mass_flow_kg_per_s
    m_b = bb.component_flows[0].mass_flow_kg_per_s
    print(f"\n  {'Model':14} {'ṁ (kg/s)':>9} {'converged':>10} {'mass bal.':>10} {'energy bal.':>12}")
    for label, r, m in (("Homogeneous", hom, m_h), ("Beggs-Brill", bb, m_b)):
        print(f"  {label:14} {m:>9.2f} {str(r.converged):>10} {r.global_balance.mass_error_pct:>9.4f}%"
              f" {r.global_energy_balance.energy_error_pct:>11.4f}%")
    print(f"  Beggs-Brill delivers {100 * (1 - m_b / m_h):.1f} % less for the same 25 bar.")

    print("\nSections with Beggs-Brill")
    section_table(bb)
    print()
    print("  • At these low mixture velocities the gas slips past the liquid everywhere:")
    print("    H_L is five to seven times the no-slip fraction C_L.")
    print("  • Uphill sections (A, C) hold the most liquid; downhill (B, D) the liquid")
    print("    drains faster and H_L is lower, so the descent recovers less elevation")
    print("    pressure than the climb cost.")
    print("  • The extra liquid also raises the two-phase friction factor: friction,")
    print("    not elevation, takes most of the 25 bar.")


if __name__ == "__main__":
    run()
