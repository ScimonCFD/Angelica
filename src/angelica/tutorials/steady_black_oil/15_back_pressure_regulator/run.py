"""
Tutorial 15 — A back-pressure regulator holds the separator pressure
====================================================================

The network of tutorial 14 with a back-pressure regulator at the start of
the separator's gas line.  The regulator holds the separator at its set
pressure whatever the export pressure downstream — as long as it can.

  Reservoir 1 ─ Well ─ WH1 ─ flowline 3 km ─┐                    ┌─ BPR (P_set) ─ gas line 4 km, 6" ─ Gas export (P_exp)
                                            ├─ Manifold ─ trunk ─ Separator
  Reservoir 2 ─ Well ─ WH2 ─ flowline 2 km ─┘    8 km, 6"        ├─ oil LCV ─ oil line 3 km, 4"  ─ Stock tank 0.2 MPa
                                                                 └─ water LCV ─ water line 1 km, 3" ─ Disposal 0.3 MPa

How the regulator is modelled (``PressureRegulator``, mode "upstream")
-----------------------------------------------------------------------
  - While it controls, its inlet node — here the separator — is held at
    the set pressure and it passes the separator's gas, whatever that is.
  - If the export is above the set point, it cannot hold the separator
    down: it is wide open (a fitting with K = 5).
  - If nothing arrives to hold the separator up (the wells dead), it is
    closed.

What to look at
---------------
  A. The separator held at 3 MPa with the export at 1.5 MPa, against
     tutorial 14 without the regulator.
  B. The export pressure: the separator stays at the set point until the
     export reaches it.
  C. The set point: what a higher separator pressure costs in production.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures import BeggsBrillCorrelation
from angelica.core.case import InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe, PressureRegulator, SeparatorOutlet, Well
from angelica.properties.black_oil import BlackOilFluid
from angelica.solvers import SteadyBlackOilSolver

# name, P_r (MPa), T (°C), API, gas gravity, GOR, WOR, J (m³/d/bar), flowline length (m) — as tutorial 14
WELLS = (
    ("Well 1", 25.0, 85.0, 40.0, 0.70, 150.0, 0.1, 6.0, 3000.0),
    ("Well 2", 22.0, 80.0, 32.0, 0.65,  40.0, 0.3, 4.0, 2000.0),
)
TUBING_LENGTH, TUBING_RISE = 2100.0, 2000.0
MANIFOLD, SEPARATOR, BPR_OUT = 100, 101, 107
GAS_EXPORT, OIL_LCV_OUT, STOCK_TANK, WATER_LCV_OUT, DISPOSAL = 102, 103, 104, 105, 106
P_EXPORT, P_SET = 1.5e6, 3.0e6


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


def build_case(p_export: float = P_EXPORT, p_set: float | None = P_SET) -> NetworkCase:
    """``p_set=None``: no regulator (tutorial 14)."""
    components, inlets, thermal, fluids = [], [], [], []
    for i, (name, P_r, T, api, gg, gor, wor, _J, L) in enumerate(WELLS):
        res, wh = 10 * (i + 1), 10 * (i + 1) + 1
        inlets.append(PressureBoundary(node_id=res, pressure_pa=P_r * 1e6))
        thermal.append(ThermalBoundary(node_id=res, temperature_c=T, bc_type="fixed_temperature"))
        fluids.append(InletFluidBC(node_id=res, api_gravity=api, gas_gravity=gg,
                                   gor_sc_m3_per_m3=gor, wor_sc_m3_per_m3=wor))
        components.append(make_well(i))
        components.append(pipe(f"flowline_{name.replace(' ', '')}", wh, MANIFOLD, 0.1016, L, U=5.0))
    gas_start = SEPARATOR if p_set is None else BPR_OUT
    if p_set is not None:
        components.append(PressureRegulator(component_id="BPR", start_node=SEPARATOR, end_node=BPR_OUT,
                                            diameter_m=0.1524, set_pressure_pa=p_set, mode="upstream"))
    components += [
        pipe("trunk", MANIFOLD, SEPARATOR, 0.1524, 8000.0),
        pipe("gas_line", gas_start, GAS_EXPORT, 0.1524, 4000.0),
        SeparatorOutlet(component_id="oil_LCV", start_node=SEPARATOR, end_node=OIL_LCV_OUT,
                        diameter_m=0.0762, phase="oil"),
        pipe("oil_line", OIL_LCV_OUT, STOCK_TANK, 0.1016, 3000.0),
        SeparatorOutlet(component_id="water_LCV", start_node=SEPARATOR, end_node=WATER_LCV_OUT,
                        diameter_m=0.0508, phase="water"),
        pipe("water_line", WATER_LCV_OUT, DISPOSAL, 0.0762, 1000.0),
    ]
    return NetworkCase(
        name=f"Back-pressure regulator, export {p_export / 1e6:.1f} MPa",
        fluid_model=well_fluid(*WELLS[1][3:7]),
        pressure_inlets=tuple(inlets),
        pressure_outlets=(PressureBoundary(node_id=GAS_EXPORT, pressure_pa=p_export),
                          PressureBoundary(node_id=STOCK_TANK, pressure_pa=0.2e6),
                          PressureBoundary(node_id=DISPOSAL, pressure_pa=0.3e6)),
        thermal_inlets=tuple(thermal),
        inlet_fluid_bcs=tuple(fluids),
        components=tuple(components),
    )


def solve(**kw):
    return SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation()).solve(build_case(**kw))


def run():
    print("=" * 78)
    print("Tutorial 15 — A back-pressure regulator holds the separator pressure (black-oil)")
    print("=" * 78)

    print(f"\nPart A — export {P_EXPORT / 1e6:.1f} MPa, with and without the regulator ({P_SET / 1e6:.1f} MPa)")
    print(f"  {'':16} {'P_sep':>7} {'Oil':>7} {'Gas':>9} {'Regulator':>22}")
    print(f"  {'':16} {'(MPa)':>7} {'(m³/d)':>7} {'(kSm³/d)':>9}")
    for label, p_set in (("no regulator", None), ("regulator", P_SET)):
        r = solve(p_set=p_set)
        s = r.separators[0]
        reg = (f"{r.regulators[0].state}, ΔP {r.regulators[0].pressure_drop_pa / 1e6:.2f} MPa"
               if r.regulators else "—")
        print(f"  {label:16} {s.pressure_pa / 1e6:7.3f} {s.oil_rate_sm3_per_day:7.1f} "
              f"{s.gas_rate_sm3_per_day / 1e3:9.1f} {reg:>22}" + ("" if r.converged else "  (not converged)"))

    print("\nPart B — export pressure (set point 3.0 MPa)")
    print(f"  {'Export':>7} {'P_sep':>7} {'State':>12} {'Oil':>7} {'Gas':>9}")
    print(f"  {'(MPa)':>7} {'(MPa)':>7} {'':>12} {'(m³/d)':>7} {'(kSm³/d)':>9}")
    for p in (0.5e6, 1.5e6, 2.5e6, 3.0e6, 3.5e6, 5.0e6):
        r = solve(p_export=p)
        s, g = r.separators[0], r.regulators[0]
        print(f"  {p / 1e6:7.1f} {s.pressure_pa / 1e6:7.3f} {g.state:>12} {s.oil_rate_sm3_per_day:7.1f} "
              f"{s.gas_rate_sm3_per_day / 1e3:9.1f}" + ("" if r.converged else "  (not converged)"))

    print("\nPart C — set point (export 1.0 MPa)")
    print(f"  {'P_set':>7} {'P_sep':>7} {'Oil':>7} {'Gas':>9} {'Well 2':>12}")
    print(f"  {'(MPa)':>7} {'(MPa)':>7} {'(m³/d)':>7} {'(kSm³/d)':>9}")
    for p_set in (1.5e6, 2.0e6, 3.0e6, 4.0e6, 5.0e6):
        r = solve(p_export=1.0e6, p_set=p_set)
        s = r.separators[0]
        print(f"  {p_set / 1e6:7.1f} {s.pressure_pa / 1e6:7.3f} {s.oil_rate_sm3_per_day:7.1f} "
              f"{s.gas_rate_sm3_per_day / 1e3:9.1f} {r.well_states[1].state:>12}"
              + ("" if r.converged else "  (not converged)"))

    print()
    print("  • The regulator decouples the separator from the export pressure: the")
    print("    wells see the same back-pressure whatever happens downstream of it.")
    print("  • A higher separator pressure costs oil from these wells; in black-oil it")
    print("    gains nothing at the stock tank (see tutorial 14).")


if __name__ == "__main__":
    run()
