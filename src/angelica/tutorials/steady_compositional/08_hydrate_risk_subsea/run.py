"""
Tutorial 08 — Hydrate risk in a subsea wet-gas line
===================================================

A wet natural gas with produced water leaves a platform at 90 bar and 40 °C
and flows 30 km along the seabed (4 °C) to shore at 60 bar.  Gas hydrates —
ice-like solids of water and light hydrocarbons — form wherever the fluid is
colder than the hydrate formation temperature at the local pressure and free
water is present.  They can block the line.

The same line is solved bare (U = 15 W/m²K) and insulated (U = 1.5 W/m²K).

Hydrate model
-------------
van der Waals–Platteeuw with Munck et al. (1988) Langmuir constants, PR
fugacities, liquid-water branch (≥ 0 °C), pure water (no inhibitor).  Every
compositional pipe reports ``hydrate_temperature_c`` and
``hydrate_subcooling_k = T_hydrate − T_fluid`` (> 0: inside the hydrate
region) at its colder end.

What to look at
---------------
  1. The hydrate curve of this gas: with 3 % propane it forms structure II
     hydrate at much lower pressure than pure methane.
  2. Bare line: the gas cools towards the seabed temperature within a few
     kilometres and enters the hydrate region.
  3. Insulated line: the gas stays warmer than the hydrate temperature all
     the way to shore.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import InletCompositionBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.compositional_fluid import CompositionalFluid
from angelica.properties.hydrate import hydrate_curve, hydrate_formation_pressure
from angelica.solvers import SteadyCompositionalSolver

COMPONENTS = ["methane", "ethane", "propane", "n-butane", "water"]
ZS = [0.81, 0.06, 0.03, 0.01, 0.09]
P_IN, P_OUT, T_IN, T_SEA = 90e5, 60e5, 40.0, 4.0
N_SECTIONS, L_SECTION, D = 6, 5000.0, 0.30


def build_case(U: float) -> NetworkCase:
    pipes = tuple(
        Pipe(component_id=f"km{5 * i}-{5 * (i + 1)}", start_node=i + 1, end_node=i + 2, diameter_m=D,
             length_m=L_SECTION, absolute_roughness_m=46e-6,
             heat_transfer_coefficient_w_per_m2k=U, ambient_temperature_c=T_SEA)
        for i in range(N_SECTIONS)
    )
    return NetworkCase(
        name=f"Subsea wet-gas line, U = {U} W/m2K",
        fluid_model=CompositionalFluid(components=COMPONENTS, default_zs=ZS),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=P_IN),),
        pressure_outlets=(PressureBoundary(node_id=N_SECTIONS + 1, pressure_pa=P_OUT),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=T_IN, bc_type="fixed_temperature"),),
        inlet_composition_bcs=(InletCompositionBC(node_id=1, zs=tuple(ZS)),),
        components=pipes,
    )


def line_table(res) -> None:
    P, T = res.node_pressures_pa, res.node_temperatures_c
    print(f"  {'Section':10} {'P_out':>6} {'T_out':>6} {'T_hyd':>6} {'subcool':>8} {'free water':>11}")
    print(f"  {'':10} {'(bar)':>6} {'(°C)':>6} {'(°C)':>6} {'(K)':>8} {'(mol frac)':>11}")
    first = None
    for i, cf in enumerate(res.component_flows):
        sub = cf.hydrate_subcooling_k
        flag = "  ← hydrate risk" if sub is not None and sub > 0 and cf.free_water_fraction > 0 else ""
        if flag and first is None:
            first = cf.label.split(":")[-1]
        print(f"  {cf.label.split(':')[-1]:10} {P[i + 2] / 1e5:>6.1f} {T[i + 2]:>6.1f} "
              f"{cf.hydrate_temperature_c:>6.1f} {sub:>+8.1f} {cf.free_water_fraction:>11.3f}{flag}")
    print(f"  ṁ = {res.component_flows[0].mass_flow_kg_per_s:.2f} kg/s  ·  converged {res.converged}"
          f"  ·  energy balance {res.global_energy_balance.energy_error_pct:.4f} %")
    print(f"  First section at risk: {first or 'none'}")


def run() -> None:
    print("=" * 72)
    print("Tutorial 08 — Hydrate risk in a subsea wet-gas line")
    print("=" * 72)

    print("\nHydrate formation curve of the gas (dry basis, free water present)")
    print(f"  {'T (°C)':>7} {'P (bar)':>9} {'structure':>10}   pure methane (bar)")
    for pt in hydrate_curve(COMPONENTS, ZS, 0.0, 20.0, 6):
        c1 = hydrate_formation_pressure(["methane"], [1.0], pt.temperature_c)
        print(f"  {pt.temperature_c:>7.1f} {pt.pressure_pa / 1e5:>9.1f} {pt.structure:>10}   {c1.pressure_pa / 1e5:>9.1f}")

    solver = SteadyCompositionalSolver()
    for U, label in ((15.0, "Bare line (U = 15 W/m²K)"), (1.5, "Insulated line (U = 1.5 W/m²K)")):
        print(f"\n{label}")
        line_table(solver.solve(build_case(U)))

    print()
    print("  • The bare line reaches seabed temperature long before shore and runs")
    print("    inside the hydrate region: it would need inhibitor injection")
    print("    (methanol, MEG) or a heated / insulated design.")
    print("  • Insulation keeps the gas above the hydrate temperature to shore.")


if __name__ == "__main__":
    run()
