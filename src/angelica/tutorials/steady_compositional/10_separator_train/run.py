"""
Tutorial 10 — A separator train (two-stage separation)
======================================================

A volatile oil with free water arrives from a production manifold through a
5 km flowline at a high-pressure (HP) separator.  The HP gas goes to export;
the water to disposal; the oil, through its level-control valve (LCV), to a
low-pressure (LP) separator, whose gas goes to a flare and whose liquid goes
to the stock tank.

  Manifold ─ flowline 5 km, 6" ─ HP separator ─┬─ gas line 1 km, 10" ─ Gas export (P_exp)
  120 bar, 70 °C                               ├─ water LCV ─ line ─ Disposal 5 bar
                                               └─ oil LCV ─ oil line 500 m, 6" ─ LP separator ─┬─ flare line 100 m, 8" ─ Flare (P_flare)
                                                                                                 └─ liquid LCV ─ tank line ─ Stock tank 1.2 bar

Both separators are ordinary nodes: their pressures follow from the export
and flare lines.

Each separator flashes the arriving streams at its pressure and temperature
(Peng-Robinson; free water by the immiscible-water model).  The separators'
"oil" rates are stock-tank volumes: the oil phase flashed once to standard
conditions (1 atm, 15.56 °C).  So the HP separator's oil rate is what the
stock tank would receive with the HP separator alone, and the LP
separator's what it receives after the second stage.

What to look at
---------------
  A. The train at an export pressure of 30 bar.
  B. The HP separator pressure: the stock-tank liquid per tonne of feed has
     a maximum — near 18 bar with the HP separator alone, near 30 bar with
     the LP stage.  The usual explanation: at high pressure much light gas
     stays in the oil and, flashing off downstream, carries intermediates
     with it; at low pressure the intermediates leave with the HP gas.
  C. The LP stage pressure: a second stage recovers more stock-tank liquid
     than one flash from the HP separator to the stock tank; the gain is
     largest at an intermediate LP pressure.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import InletCompositionBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe, SeparatorOutlet
from angelica.properties.compositional_fluid import CompositionalFluid
from angelica.solvers import SteadyCompositionalSolver

COMPONENTS = ("methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane", "n-decane", "water")
FEED = (0.62, 0.09, 0.07, 0.05, 0.03, 0.03, 0.08, 0.03)
MANIFOLD, HP, EXPORT, WATER_OUT, DISPOSAL, OIL_OUT, LP, LIQ_OUT, STOCK_TANK, FLARE = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10


def pipe(cid, a, b, D, L, U=0.0) -> Pipe:
    return Pipe(component_id=cid, start_node=a, end_node=b, diameter_m=D, length_m=L,
                absolute_roughness_m=46e-6, heat_transfer_coefficient_w_per_m2k=U, ambient_temperature_c=15.0)


def build_case(p_export: float = 30e5, p_flare: float = 5e5) -> NetworkCase:
    fluid = CompositionalFluid(COMPONENTS, FEED)
    return NetworkCase(
        name=f"Separator train, export {p_export / 1e5:.0f} bar, flare {p_flare / 1e5:.0f} bar",
        fluid_model=fluid,
        pressure_inlets=(PressureBoundary(MANIFOLD, 120e5),),
        pressure_outlets=(PressureBoundary(EXPORT, p_export), PressureBoundary(DISPOSAL, 5e5),
                          PressureBoundary(FLARE, p_flare), PressureBoundary(STOCK_TANK, 1.2e5)),
        thermal_inlets=(ThermalBoundary(MANIFOLD, 70.0),),
        inlet_composition_bcs=(InletCompositionBC(MANIFOLD, FEED),),
        components=(
            pipe("flowline", MANIFOLD, HP, 0.1524, 5000.0, U=3.0),
            pipe("gas_line", HP, EXPORT, 0.2540, 1000.0),
            SeparatorOutlet(component_id="HP_water_LCV", start_node=HP, end_node=WATER_OUT,
                            diameter_m=0.0508, phase="water"),
            pipe("water_line", WATER_OUT, DISPOSAL, 0.0508, 300.0),
            SeparatorOutlet(component_id="HP_oil_LCV", start_node=HP, end_node=OIL_OUT,
                            diameter_m=0.0508, phase="oil"),
            pipe("oil_line", OIL_OUT, LP, 0.1524, 500.0),
            pipe("flare_line", LP, FLARE, 0.2032, 100.0),
            SeparatorOutlet(component_id="LP_liquid_LCV", start_node=LP, end_node=LIQ_OUT,
                            diameter_m=0.0508, phase="liquid"),
            pipe("tank_line", LIQ_OUT, STOCK_TANK, 0.1524, 100.0),
        ),
    )


def solve(p_export: float = 30e5, p_flare: float = 5e5):
    res = SteadyCompositionalSolver().solve(build_case(p_export, p_flare))
    seps = {s.node_id: s for s in res.separators}
    return res, seps[HP], seps[LP]


def run():
    print("=" * 74)
    print("Tutorial 10 — A separator train (Peng-Robinson)")
    print("=" * 74)

    print("\nPart A — export 30 bar, flare 5 bar")
    res, hp, lp = solve()
    print(f"  Converged: {res.converged}")
    print(f"  {'':14} {'P':>6} {'T':>6} {'Feed':>6} {'Gas':>6} {'Oil':>6} {'Water':>6} "
          f"{'Gas':>8} {'Oil':>6}")
    print(f"  {'':14} {'(bar)':>6} {'(°C)':>6} {'(kg/s)':>6} {'(kg/s)':>6} {'(kg/s)':>6} {'(kg/s)':>6} "
          f"{'(kSm³/d)':>8} {'(m³/d)':>6}")
    for name, s in (("HP separator", hp), ("LP separator", lp)):
        print(f"  {name:14} {s.pressure_pa / 1e5:6.2f} {s.temperature_c:6.1f} {s.feed_mass_flow_kg_per_s:6.2f} "
              f"{s.gas_mass_flow_kg_per_s:6.2f} {s.oil_mass_flow_kg_per_s:6.2f} {s.water_mass_flow_kg_per_s:6.3f} "
              f"{s.gas_rate_sm3_per_day / 1e3:8.1f} {s.oil_rate_sm3_per_day:6.1f}")
    for s in (hp, lp):
        for o in s.outlets:
            print(f"  {o.component_id:14}: {o.mass_flow_kg_per_s:6.3f} kg/s, LCV drop {o.valve_pressure_drop_pa / 1e5:+7.2f} bar")
    print(f"  Stock tank: {lp.oil_rate_sm3_per_day:.1f} m³/d after two stages; "
          f"{hp.oil_rate_sm3_per_day:.1f} m³/d if the HP oil went straight to the tank")
    gb, geb = res.global_balance, res.global_energy_balance
    print(f"  Mass balance {gb.mass_inlet_kg_per_s - gb.mass_outlet_kg_per_s:+.4f} kg/s, "
          f"energy balance {geb.energy_error_kw:+.3f} kW")

    print("\nPart B — HP separator pressure (flare 5 bar); liquid per tonne of feed")
    print(f"  {'Export':>7} {'P_HP':>6} {'T_HP':>6} {'Feed':>7} {'One stage':>10} {'Two stages':>11}")
    print(f"  {'(bar)':>7} {'(bar)':>6} {'(°C)':>6} {'(kg/s)':>7} {'(m³/t)':>10} {'(m³/t)':>11}")
    for p_exp in (5e5, 15e5, 30e5, 50e5, 80e5):
        _r, hp, lp = solve(p_exp)
        t_per_day = hp.feed_mass_flow_kg_per_s * 86.4
        print(f"  {p_exp / 1e5:7.0f} {hp.pressure_pa / 1e5:6.1f} {hp.temperature_c:6.1f} "
              f"{hp.feed_mass_flow_kg_per_s:7.2f} {hp.oil_rate_sm3_per_day / t_per_day:10.4f} "
              f"{lp.oil_rate_sm3_per_day / t_per_day:11.4f}")

    print("\nPart C — LP stage pressure, set by the flare pressure (export 30 bar)")
    print(f"  {'P_LP':>6} {'T_LP':>6} {'Two stages':>11} {'vs one stage':>13}")
    print(f"  {'(bar)':>6} {'(°C)':>6} {'(m³/t)':>11}")
    for p_flare in (2e5, 5e5, 10e5):
        _r, hp, lp = solve(30e5, p_flare)
        t_per_day = hp.feed_mass_flow_kg_per_s * 86.4
        two, one = lp.oil_rate_sm3_per_day / t_per_day, hp.oil_rate_sm3_per_day / t_per_day
        print(f"  {lp.pressure_pa / 1e5:6.2f} {lp.temperature_c:6.1f} {two:11.4f} {100.0 * (two / one - 1.0):+12.1f}%")


if __name__ == "__main__":
    run()
