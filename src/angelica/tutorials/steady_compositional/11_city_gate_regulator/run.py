"""
Tutorial 11 — A pressure-reducing station (city gate)
=====================================================

A dry natural gas arrives through a 20 km transmission line and is let
down to a 20 bar distribution network by a pressure-reducing regulator.

  Supply 70 bar, 20 °C ─ transmission 20 km, 8" ─ [heater] ─ Regulator (20 bar) ─ distribution 2 km, 6" ─ Consumers 5 bar

How the regulator is modelled (``PressureRegulator``, mode "downstream")
-------------------------------------------------------------------------
  - While it controls, its outlet node is held at the set pressure and it
    passes whatever the distribution side takes.
  - If the supply cannot keep its inlet above the set point it is wide
    open (a fitting with K = 5); if the outlet side is already at or
    above the set point it is closed.
  - The gas crosses it at constant enthalpy, so it cools (Joule-Thomson).

What to look at
---------------
  A. The let-down: the gas leaves the regulator far below 0 °C.
  B. Preheating: the heater duty that keeps the regulator outlet above
     0 °C.
  C. The supply pressure: the regulator holds 20 bar until the supply can
     no longer feed the demand, then it is wide open.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import InletCompositionBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import HeatSource, Pipe, PressureRegulator
from angelica.properties.compositional_fluid import CompositionalFluid
from angelica.solvers import SteadyCompositionalSolver

COMPONENTS = ("methane", "ethane", "propane", "n-butane", "nitrogen", "carbon dioxide")
GAS = (0.90, 0.05, 0.02, 0.01, 0.01, 0.01)
SUPPLY, HEATER_IN, REG_IN, REG_OUT, CONSUMERS = 1, 2, 3, 4, 5
P_SET = 20e5


def pipe(cid, a, b, D, L) -> Pipe:
    return Pipe(component_id=cid, start_node=a, end_node=b, diameter_m=D, length_m=L,
                absolute_roughness_m=46e-6, heat_transfer_coefficient_w_per_m2k=2.0, ambient_temperature_c=10.0)


def build_case(p_supply: float = 70e5, heater_w: float = 0.0, p_consumers: float = 5e5) -> NetworkCase:
    return NetworkCase(
        name=f"City gate, supply {p_supply / 1e5:.0f} bar, heater {heater_w / 1e3:.0f} kW",
        fluid_model=CompositionalFluid(COMPONENTS, GAS),
        pressure_inlets=(PressureBoundary(SUPPLY, p_supply),),
        pressure_outlets=(PressureBoundary(CONSUMERS, p_consumers),),
        thermal_inlets=(ThermalBoundary(SUPPLY, 20.0),),
        inlet_composition_bcs=(InletCompositionBC(SUPPLY, GAS),),
        components=(
            pipe("transmission", SUPPLY, HEATER_IN, 0.2032, 20000.0),
            HeatSource(component_id="heater", start_node=HEATER_IN, end_node=REG_IN, diameter_m=0.2032,
                       power_w=heater_w),
            PressureRegulator(component_id="regulator", start_node=REG_IN, end_node=REG_OUT, diameter_m=0.2032,
                              set_pressure_pa=P_SET),
            pipe("distribution", REG_OUT, CONSUMERS, 0.1524, 2000.0),
        ),
    )


def solve(**kw):
    return SteadyCompositionalSolver().solve(build_case(**kw))


def run():
    print("=" * 78)
    print("Tutorial 11 — A pressure-reducing station (compositional, Peng-Robinson)")
    print("=" * 78)

    print("\nPart A — the let-down, no preheating")
    res = solve()
    reg = res.regulators[0]
    T = res.node_temperatures_c
    print(f"  Converged: {res.converged}")
    print(f"  Regulator {reg.state}: {reg.inlet_pressure_pa / 1e5:.2f} → {reg.outlet_pressure_pa / 1e5:.2f} bar, "
          f"{reg.mass_flow_kg_per_s:.3f} kg/s")
    print(f"  Temperature {T[REG_IN]:.2f} → {T[REG_OUT]:.2f} °C "
          f"({(T[REG_OUT] - T[REG_IN]) / (reg.pressure_drop_pa / 1e5):+.3f} K/bar)")
    print(f"  Consumers receive the gas at {T[CONSUMERS]:.2f} °C after 2 km of distribution line")
    geb = res.global_energy_balance
    print(f"  Energy balance residual {geb.energy_error_kw:+.3f} kW")

    print("\nPart B — preheating before the regulator")
    print(f"  {'Heater (kW)':>11} {'T in (°C)':>10} {'T out (°C)':>11} {'Flow (kg/s)':>12}")
    for kw in (0.0, 200.0, 400.0, 600.0):
        r = solve(heater_w=kw * 1e3)
        Tn = r.node_temperatures_c
        print(f"  {kw:11.0f} {Tn[REG_IN]:10.2f} {Tn[REG_OUT]:11.2f} {r.regulators[0].mass_flow_kg_per_s:12.3f}")

    print("\nPart C — supply pressure")
    print(f"  {'Supply (bar)':>12} {'State':>12} {'P in (bar)':>11} {'P out (bar)':>12} {'Flow (kg/s)':>12}")
    for p in (25e5, 30e5, 35e5, 40e5, 50e5, 70e5):
        r = solve(p_supply=p)
        g = r.regulators[0]
        print(f"  {p / 1e5:12.0f} {g.state:>12} {g.inlet_pressure_pa / 1e5:11.2f} {g.outlet_pressure_pa / 1e5:12.2f} "
              f"{g.mass_flow_kg_per_s:12.3f}" + ("" if r.converged else "  (not converged)"))

    print()
    print("  • While the regulator controls, the consumers see 20 bar.  The flow still")
    print("    rises a little with the supply pressure: a larger let-down cools the gas")
    print("    more, and the colder, denser gas passes the distribution line faster.")
    print("  • Below about 36 bar the supply line limits the flow: the regulator is")
    print("    wide open and the consumers get less than 20 bar.")


if __name__ == "__main__":
    run()
