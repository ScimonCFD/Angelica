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
  2. Section table: flow pattern, no-slip C_L, in-situ H_L, and the split
     of the pressure drop between elevation and friction.
  3. Liquid inventory: the line holds far more liquid than the homogeneous
     model predicts — what a pig would push out, or a slug catcher see.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures import BeggsBrillCorrelation
from angelica.closures.beggs_brill import beggs_brill_gradient
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


class _LinkView:
    """Mean conditions of one solved section, in the form the fluid model reads."""

    class _Node:
        def __init__(self, p: float) -> None:
            self.pressure_pa = p

    def __init__(self, p_in: float, p_out: float, t_mean: float) -> None:
        self.start_node, self.end_node = self._Node(p_in), self._Node(p_out)
        self.temperature_c = t_mean
        self.zs = tuple(ZS)


def section_table(case: NetworkCase, res) -> None:
    fluid = case.fluid_model
    area = math.pi * D * D / 4.0
    P, T = res.node_pressures_pa, res.node_temperatures_c
    m = res.component_flows[0].mass_flow_kg_per_s
    print(f"  {'Sec':3} {'dz':>5} {'P_out':>6} {'T_out':>6} {'pattern':>13} {'C_L':>6} {'H_L':>6}"
          f" {'ΔP elev':>8} {'ΔP fric':>8} {'liquid':>7}")
    print(f"  {'':3} {'(m)':>5} {'(bar)':>6} {'(°C)':>6} {'':>13} {'':>6} {'':>6}"
          f" {'(bar)':>8} {'(bar)':>8} {'(m³)':>7}")
    inv_ns = inv_bb = 0.0
    for i, (name, L, dz) in enumerate(SECTIONS):
        link = _LinkView(P[i + 1], P[i + 2], 0.5 * (T[i + 1] + T[i + 2]))
        CL, rl, rg, ml, mg, sig = fluid.beggs_brill_phase_properties(link)
        v_m = m / ((CL * rl + (1 - CL) * rg) * area)
        st = beggs_brill_gradient(v_m, CL, rl, rg, ml, mg, D, ROUGH, math.asin(dz / L), sig)
        inv_ns += CL * area * L
        inv_bb += st.holdup * area * L
        print(f"  {name:3} {dz:>+5.0f} {P[i + 2] / 1e5:>6.2f} {T[i + 2]:>6.1f} {st.regime:>13} {CL:>6.3f}"
              f" {st.holdup:>6.3f} {st.dP_dz_gravity * L / 1e5:>8.2f} {st.dP_dz_friction * L / 1e5:>8.2f}"
              f" {st.holdup * area * L:>7.1f}")
    print(f"  Liquid in the line: {inv_bb:.0f} m³ with slip  ·  {inv_ns:.0f} m³ at the no-slip fraction")


def run() -> None:
    print("=" * 76)
    print("Tutorial 07 — Hilly gas-condensate line: homogeneous vs Beggs-Brill")
    print("=" * 76)

    hom = SteadyCompositionalSolver().solve(build_case())
    case_bb = build_case()
    bb = SteadyCompositionalSolver(turbulent_pipe_correlation=BeggsBrillCorrelation()).solve(case_bb)

    m_h = hom.component_flows[0].mass_flow_kg_per_s
    m_b = bb.component_flows[0].mass_flow_kg_per_s
    print(f"\n  {'Model':14} {'ṁ (kg/s)':>9} {'converged':>10} {'mass bal.':>10} {'energy bal.':>12}")
    for label, r, m in (("Homogeneous", hom, m_h), ("Beggs-Brill", bb, m_b)):
        print(f"  {label:14} {m:>9.2f} {str(r.converged):>10} {r.global_balance.mass_error_pct:>9.4f}%"
              f" {r.global_energy_balance.energy_error_pct:>11.4f}%")
    print(f"  Beggs-Brill delivers {100 * (1 - m_b / m_h):.1f} % less for the same 25 bar.")

    print("\nSections with Beggs-Brill")
    section_table(case_bb, bb)
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
