"""
Tutorial 16 — A heater station set by its outlet temperature
============================================================

A 50 km, 8 in. line carrying a heavy crude (22 °API) that leaves the field
at 80 °C, with a heater station half-way.  The station is specified the way
an operator specifies it: by the temperature the oil must leave it at, not
by its duty.  The solver finds the duty.

  Inlet 10 MPa, 80 °C ─ pipe 25 km ─ [heater: outlet 80 °C, max duty] ─ pipe 25 km ─ Outlet 2 MPa

How the heater is modelled (``HeatSource`` with ``outlet_temperature_c``)
---------------------------------------------------------------------------
  - Its duty is the enthalpy the stream needs to reach the set outlet
    temperature, ṁ·(h_out − h_in), with the black-oil enthalpy model.
  - With ``max_duty_w`` the duty is capped: beyond it the heater works at
    its maximum and the oil leaves colder than the set point ("duty
    limited").
  - Without ``outlet_temperature_c`` it is a fixed-duty heater
    (``power_w``), as before.

What to look at
---------------
  A. With and without the station: the flow the same 8 MPa drives.
  B. The ambient temperature: the duty the station needs to hold 80 °C.
  C. A station with a 2.5 MW limit: where it stops holding the set point.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import HeatSource, Pipe
from angelica.properties.black_oil import BlackOilFluid
from angelica.solvers import SteadyBlackOilSolver

CRUDE = {"api_gravity": 22.0, "gas_gravity": 0.65, "gor_sc_m3_per_m3": 25.0, "wor_sc_m3_per_m3": 0.5}
INLET, HEATER_IN, HEATER_OUT, OUTLET = 1, 2, 3, 4
P_IN, P_OUT, T_IN = 10.0e6, 2.0e6, 80.0


def pipe(cid, a, b, t_amb) -> Pipe:
    return Pipe(component_id=cid, start_node=a, end_node=b, diameter_m=0.2, length_m=25_000.0,
                absolute_roughness_m=46e-6, heat_transfer_coefficient_w_per_m2k=3.0,
                ambient_temperature_c=t_amb, n_thermal_segments=20)


def build_case(t_amb: float = 10.0, t_set: float | None = 80.0, max_duty_w: float = float("inf")) -> NetworkCase:
    """``t_set=None``: no heating (the station is a 0 W heater)."""
    heater = HeatSource(component_id="station", start_node=HEATER_IN, end_node=HEATER_OUT, diameter_m=0.2,
                        pressure_drop_mode="rated", pressure_drop_pa=0.5e5, rated_mass_flow_kg_per_s=50.0,
                        outlet_temperature_c=t_set, max_duty_w=max_duty_w)
    return NetworkCase(
        name=f"Heater station, ambient {t_amb:g} °C",
        fluid_model=BlackOilFluid(**CRUDE),
        pressure_inlets=(PressureBoundary(INLET, P_IN),),
        pressure_outlets=(PressureBoundary(OUTLET, P_OUT),),
        thermal_inlets=(ThermalBoundary(INLET, T_IN),),
        components=(pipe("first_half", INLET, HEATER_IN, t_amb), heater,
                    pipe("second_half", HEATER_OUT, OUTLET, t_amb)),
    )


def solve(**kw):
    return SteadyBlackOilSolver().solve(build_case(**kw))


def run():
    print("=" * 78)
    print("Tutorial 16 — A heater station set by its outlet temperature (black-oil)")
    print("=" * 78)

    print("\nPart A — the line with and without the station (ambient 10 °C)")
    print(f"  {'':16} {'Flow':>7} {'T before':>9} {'T after':>8} {'T outlet':>9} {'Duty':>8}")
    print(f"  {'':16} {'(kg/s)':>7} {'(°C)':>9} {'(°C)':>8} {'(°C)':>9} {'(MW)':>8}")
    for label, t_set in (("no heating", None), ("station 80 °C", 80.0)):
        r = solve(t_set=t_set)
        h, T = r.heaters[0], r.node_temperatures_c
        print(f"  {label:16} {h.mass_flow_kg_per_s:7.2f} {T[HEATER_IN]:9.2f} {T[HEATER_OUT]:8.2f} "
              f"{T[OUTLET]:9.2f} {h.duty_w / 1e6:8.3f}" + ("" if r.converged else "  (not converged)"))
    geb = r.global_energy_balance
    print(f"  Energy balance with the station: {geb.heat_sources_kw:.1f} kW in, residual "
          f"{geb.energy_error_kw:+.3f} kW ({geb.energy_error_pct:+.1e} %)")

    print("\nPart B — ambient temperature (station holds 80 °C)")
    print(f"  {'Ambient':>8} {'Flow':>7} {'T before':>9} {'Duty':>8} {'T outlet':>9}")
    print(f"  {'(°C)':>8} {'(kg/s)':>7} {'(°C)':>9} {'(MW)':>8} {'(°C)':>9}")
    for t_amb in (25.0, 20.0, 15.0, 10.0, 5.0):
        r = solve(t_amb=t_amb)
        h, T = r.heaters[0], r.node_temperatures_c
        print(f"  {t_amb:8.0f} {h.mass_flow_kg_per_s:7.2f} {T[HEATER_IN]:9.2f} {h.duty_w / 1e6:8.3f} "
              f"{T[OUTLET]:9.2f}" + ("" if r.converged else "  (not converged)"))

    print("\nPart C — a station limited to 2.5 MW")
    print(f"  {'Ambient':>8} {'State':>22} {'Duty':>7} {'T after':>8} {'T outlet':>9}")
    print(f"  {'(°C)':>8} {'':>22} {'(MW)':>7} {'(°C)':>8} {'(°C)':>9}")
    for t_amb in (25.0, 20.0, 15.0, 10.0, 5.0):
        r = solve(t_amb=t_amb, max_duty_w=2.5e6)
        h, T = r.heaters[0], r.node_temperatures_c
        print(f"  {t_amb:8.0f} {h.state:>22} {h.duty_w / 1e6:7.3f} {T[HEATER_OUT]:8.2f} {T[OUTLET]:9.2f}"
              + ("" if r.converged else "  (not converged)"))

    print()
    print("  • The station turns the second half into a copy of the first.  The oil in")
    print("    it is about 16 K warmer (mean 70 against 54 °C) and its viscosity about")
    print("    3.5 against 6 mPa·s, so the same pressure drop moves 2.3 % more of it.")
    print("  • The colder the ground, the more the first half loses and the more the")
    print("    station must give back: about 40 kW per kelvin of ambient here.")
    print("  • A station at its limit no longer holds the set point; the solution says")
    print("    so (\"duty limited\") instead of silently delivering colder oil.")


if __name__ == "__main__":
    run()
