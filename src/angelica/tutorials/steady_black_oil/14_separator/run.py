"""
Tutorial 14 — A three-phase separator in the network
====================================================

Wells 1 and 2 of tutorial 12 deliver through a trunk line to a production
separator.  The separator is not a boundary: its pressure follows from the
gas export line downstream, and its gas, oil and water leave through their
own lines.

  Reservoir 1 ─ Well ─ WH1 ─ flowline 3 km ─┐                    ┌─ gas line 4 km, 6" ─ Gas export (P_exp)
                                            ├─ Manifold ─ trunk ─ Separator
  Reservoir 2 ─ Well ─ WH2 ─ flowline 2 km ─┘    8 km, 6"        ├─ oil LCV ─ oil line 3 km, 4"  ─ Stock tank 0.2 MPa
                                                                 └─ water LCV ─ water line 1 km, 3" ─ Disposal 0.3 MPa

How the separator is modelled (``SeparatorOutlet``)
---------------------------------------------------
  - The separator node mixes the streams arriving there and splits them at
    its pressure and temperature: the oil leaves with the gas still
    dissolved at the separator, Rs(P, T); the free gas, GOR − Rs, leaves
    through the gas line; the water through its own outlet.
  - Each ``SeparatorOutlet`` link carries its phases whatever its pressure
    difference: it stands for the level-control valve (LCV), whose pressure
    drop is reported.  A negative drop means the liquid could not reach its
    outlet without a pump.

What to look at
---------------
  A. The split at the base export pressure.
  B. The export pressure sets the separator pressure: at lower pressure
     less gas stays dissolved in the oil (Rs falls), and the wells produce
     more.
  C. A liquid line too small for its flow: the LCV would need a pump.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures import BeggsBrillCorrelation
from angelica.closures.well_model import well_curve
from angelica.core.case import InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe, SeparatorOutlet, Well
from angelica.properties.black_oil import BlackOilFluid
from angelica.solvers import SteadyBlackOilSolver

# name, P_r (MPa), T (°C), API, gas gravity, GOR, WOR, J (m³/d/bar), flowline length (m)
WELLS = (
    ("Well 1", 25.0, 85.0, 40.0, 0.70, 150.0, 0.1, 6.0, 3000.0),
    ("Well 2", 22.0, 80.0, 32.0, 0.65,  40.0, 0.3, 4.0, 2000.0),
)
TUBING_LENGTH, TUBING_RISE = 2100.0, 2000.0
MANIFOLD, SEPARATOR = 100, 101
GAS_EXPORT, OIL_LCV_OUT, STOCK_TANK, WATER_LCV_OUT, DISPOSAL = 102, 103, 104, 105, 106
P_EXPORT = 1.5e6


def well_fluid(api, gg, gor, wor) -> BlackOilFluid:
    return BlackOilFluid(api_gravity=api, gas_gravity=gg, gor_sc_m3_per_m3=gor, wor_sc_m3_per_m3=wor)


def make_well(i: int) -> Well:
    name = WELLS[i][0].replace(" ", "")
    return Well(component_id=f"well_{name}", start_node=10 * (i + 1), end_node=10 * (i + 1) + 1,
                diameter_m=0.0762, productivity_index_m3_per_day_per_bar=WELLS[i][7], ipr_model="vogel",
                tubing_length_m=TUBING_LENGTH, tubing_height_change_m=TUBING_RISE)


def pipe(cid, a, b, D, L, U=3.0) -> Pipe:
    return Pipe(component_id=cid, start_node=a, end_node=b, diameter_m=D, length_m=L,
                absolute_roughness_m=46e-6, heat_transfer_coefficient_w_per_m2k=U, ambient_temperature_c=15.0)


def build_case(p_export: float = P_EXPORT, oil_line_d: float = 0.1016) -> NetworkCase:
    components, inlets, thermal, fluids = [], [], [], []
    for i, (name, P_r, T, api, gg, gor, wor, _J, L) in enumerate(WELLS):
        res, wh = 10 * (i + 1), 10 * (i + 1) + 1
        inlets.append(PressureBoundary(node_id=res, pressure_pa=P_r * 1e6))
        thermal.append(ThermalBoundary(node_id=res, temperature_c=T, bc_type="fixed_temperature"))
        fluids.append(InletFluidBC(node_id=res, api_gravity=api, gas_gravity=gg,
                                   gor_sc_m3_per_m3=gor, wor_sc_m3_per_m3=wor))
        components.append(make_well(i))
        components.append(pipe(f"flowline_{name.replace(' ', '')}", wh, MANIFOLD, 0.1016, L, U=5.0))
    components += [
        pipe("trunk", MANIFOLD, SEPARATOR, 0.1524, 8000.0),
        pipe("gas_line", SEPARATOR, GAS_EXPORT, 0.1524, 4000.0),
        SeparatorOutlet(component_id="oil_LCV", start_node=SEPARATOR, end_node=OIL_LCV_OUT,
                        diameter_m=0.0762, phase="oil"),
        pipe("oil_line", OIL_LCV_OUT, STOCK_TANK, oil_line_d, 3000.0),
        SeparatorOutlet(component_id="water_LCV", start_node=SEPARATOR, end_node=WATER_LCV_OUT,
                        diameter_m=0.0508, phase="water"),
        pipe("water_line", WATER_LCV_OUT, DISPOSAL, 0.0762, 1000.0),
    ]
    return NetworkCase(
        name=f"Separator, export {p_export / 1e6:.1f} MPa",
        fluid_model=well_fluid(*WELLS[1][3:7]),
        pressure_inlets=tuple(inlets),
        pressure_outlets=(PressureBoundary(node_id=GAS_EXPORT, pressure_pa=p_export),
                          PressureBoundary(node_id=STOCK_TANK, pressure_pa=0.2e6),
                          PressureBoundary(node_id=DISPOSAL, pressure_pa=0.3e6)),
        thermal_inlets=tuple(thermal),
        inlet_fluid_bcs=tuple(fluids),
        components=tuple(components),
    )


def well_oil_rates(result) -> list[float]:
    """Stock-tank oil of each well (m³/d)."""
    flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
    rates = []
    for i, (name, P_r, T, api, gg, gor, wor, *_r) in enumerate(WELLS):
        curve = well_curve(make_well(i), well_fluid(api, gg, gor, wor), P_r * 1e6, T)
        rates.append(flows[f"well_{name.replace(' ', '')}"] * 86_400.0 / curve.mass_per_oil_sc_m3)
    return rates


def wells_oil(result) -> float:
    return sum(well_oil_rates(result))


def feed_gor(result) -> float:
    """Gas-oil ratio of the separator feed (Sm³/m³): the wells' GORs weighted by their oil."""
    q = well_oil_rates(result)
    return sum(qi * w[5] for qi, w in zip(q, WELLS)) / sum(q)


def run():
    solver = SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation())
    print("=" * 78)
    print("Tutorial 14 — A three-phase separator in the network (black-oil, Beggs-Brill)")
    print("=" * 78)

    print(f"\nPart A — the split at an export pressure of {P_EXPORT / 1e6:.1f} MPa")
    res = solver.solve(build_case())
    sep = res.separators[0]
    P = res.node_pressures_pa
    print(f"  Converged: {res.converged}")
    print(f"  Separator: {sep.pressure_pa / 1e6:.3f} MPa, {sep.temperature_c:.1f} °C "
          f"(export {P[GAS_EXPORT] / 1e6:.2f} MPa + gas-line friction)")
    print(f"  {'Phase':8} {'kg/s':>8} {'Standard rate':>16}")
    print(f"  {'feed':8} {sep.feed_mass_flow_kg_per_s:8.2f}")
    print(f"  {'gas':8} {sep.gas_mass_flow_kg_per_s:8.2f} {sep.gas_rate_sm3_per_day:>11.0f} Sm³/d")
    print(f"  {'oil':8} {sep.oil_mass_flow_kg_per_s:8.2f} {sep.oil_rate_sm3_per_day:>11.1f} Sm³/d")
    print(f"  {'water':8} {sep.water_mass_flow_kg_per_s:8.2f} {sep.water_rate_m3_per_day:>11.1f} m³/d")
    print(f"  Wells' stock-tank oil {wells_oil(res):.1f} m³/d (the separator's oil, as it must be)")
    for o in sep.outlets:
        print(f"  {o.component_id:9}: {o.mass_flow_kg_per_s:6.2f} kg/s, LCV pressure drop "
              f"{o.valve_pressure_drop_pa / 1e6:+.3f} MPa")
    gb, geb = res.global_balance, res.global_energy_balance
    print(f"  Mass balance {gb.mass_inlet_kg_per_s - gb.mass_outlet_kg_per_s:+.4f} kg/s, "
          f"energy balance {geb.energy_error_kw:+.3f} kW")

    print("\nPart B — export pressure")
    print(f"  {'Export':>7} {'P_sep':>7} {'T_sep':>6} {'Oil':>7} {'Free gas':>9} {'GOR':>5} {'Rs':>5} {'Free':>5} {'Well 1':>7} {'Well 2':>7}")
    print(f"  {'(MPa)':>7} {'(MPa)':>7} {'(°C)':>6} {'(m³/d)':>7} {'(kSm³/d)':>9} {'feed':>5} {'sep.':>5} {'gas':>5} {'(m³/d)':>7} {'(m³/d)':>7}")
    for p_exp in (0.5e6, 1.0e6, 1.5e6, 3.0e6, 5.0e6, 8.0e6, 10.0e6):
        r = solver.solve(build_case(p_exp))
        s = r.separators[0]
        if not r.converged:
            print(f"  {p_exp / 1e6:7.1f}  not converged (see the README: the oil line has no solution)")
            continue
        q_o = s.oil_rate_sm3_per_day
        gor = feed_gor(r)
        free = s.gas_rate_sm3_per_day / q_o          # free gas per m³ of oil
        print(f"  {p_exp / 1e6:7.1f} {s.pressure_pa / 1e6:7.3f} {s.temperature_c:6.1f} {q_o:7.1f} "
              f"{s.gas_rate_sm3_per_day / 1e3:9.1f} {gor:5.0f} {gor - free:5.1f} {free:5.0f} "
              + " ".join(f"{q:7.1f}" for q in well_oil_rates(r)))
    print("  (GOR and Rs in Sm³/m³.  At 5 MPa Well 2, the low-GOR well, is close to")
    print("  dying and the feed's GOR rises to Well 1's; above it Well 2 is dead.)")

    print("\nPart C — a 2-inch oil line")
    r = solver.solve(build_case(oil_line_d=0.0508))
    s = r.separators[0]
    oil = next(o for o in s.outlets if o.phase == "oil")
    print(f"  Separator {s.pressure_pa / 1e6:.3f} MPa; oil LCV outlet "
          f"{r.node_pressures_pa[OIL_LCV_OUT] / 1e6:.3f} MPa: LCV drop {oil.valve_pressure_drop_pa / 1e6:+.3f} MPa")
    print(f"  → {s.message}")

    print()
    print("  • The separator's pressure is set by what is downstream of its gas;")
    print("    its liquid outlets only take what the flash leaves as liquid.")
    print("  • In black-oil the stock-tank oil is the wells' oil whatever the")
    print("    separator pressure: the model has no stage-separation effect")
    print("    (see compositional tutorial 10).")


if __name__ == "__main__":
    run()
