"""
Tutorial: Non-Isothermal Black-Oil Hot-Oil Pipeline
====================================================
A heated crude-oil pipeline where temperature-dependent PVT properties
(viscosity, density, solution GOR) materially affect the hydraulic result.

Geometry
--------
  [Node 1] ──────────── Pipe (20 km, D = 0.25 m) ──────────── [Node 2]
   P = 8 MPa                                                    P = 2 MPa
   T = 80 °C                                                    T = ?
   Ambient = 10 °C, U = 2 W/m²·K

Physics
-------
The black-oil solver already couples hydraulics with the energy equation:
each outer iteration updates nodal temperatures, and the PVT correlations
(Standing bubble-point, Beggs-Robinson viscosity, Lee-Gonzalez-Eakin gas
viscosity) re-evaluate at the local (P, T).  This tutorial makes that
coupling explicit by comparing:

  (a) Non-isothermal: T_in = 80 °C, pipe cools toward ambient = 10 °C.
  (b) Isothermal cold: T_in = 10 °C (isothermal at ambient).

The hot case benefits from low viscosity at elevated temperature, giving
significantly higher flow for the same pressure drop.

Temperature profile verification
---------------------------------
For a single pipe at steady state, neglecting viscous dissipation and
assuming nearly constant Cp and mass flow:

    T(L) = T_amb + (T_in - T_amb) * exp(−U · π · D · L / (ṁ · Cp))

This analytical formula is used as the quantitative benchmark (tolerance
±2 °C to allow for variable Cp along the pipe).
"""
import math
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.black_oil import BlackOilFluid
from angelica.solvers import SteadyBlackOilSolver

# ── Parameters ────────────────────────────────────────────────────────────────
API     = 32.0
GAS_GR  = 0.65
GOR_SC  = 25.0   # m³/m³ at standard conditions
WOR_SC  = 0.5    # m³/m³ at standard conditions

D       = 0.25   # m
L       = 20_000.0  # m
U       = 2.0    # W/m²·K
T_AMB   = 10.0   # °C
T_HOT   = 80.0   # °C
P_IN    = 8.0e6  # Pa
P_OUT   = 2.0e6  # Pa

# ── Fluids ────────────────────────────────────────────────────────────────────
fluid_hot = BlackOilFluid(
    api_gravity             = API,
    gas_gravity             = GAS_GR,
    gor_sc_m3_per_m3        = GOR_SC,
    wor_sc_m3_per_m3        = WOR_SC,
    reference_pressure_pa   = 5.0e6,
    reference_temperature_c = T_HOT,
)

fluid_cold = BlackOilFluid(
    api_gravity             = API,
    gas_gravity             = GAS_GR,
    gor_sc_m3_per_m3        = GOR_SC,
    wor_sc_m3_per_m3        = WOR_SC,
    reference_pressure_pa   = 5.0e6,
    reference_temperature_c = T_AMB,
)

# ── PVT comparison at inlet vs outlet temperatures ────────────────────────────
pvt_80  = fluid_hot.pvt(P_IN,  T_HOT)
pvt_10  = fluid_hot.pvt(P_OUT, T_AMB)

print("=" * 65)
print("PVT comparison: same fluid at different (P, T) conditions")
print(f"{'':40s}  {'T=80°C,8MPa':>11s}  {'T=10°C,2MPa':>11s}")
print(f"{'Bubble point Pb (MPa)':40s}  {pvt_80.bubble_point_pa/1e6:>11.2f}  {pvt_10.bubble_point_pa/1e6:>11.2f}")
print(f"{'Solution GOR Rs (m³/m³)':40s}  {pvt_80.rs_m3_per_m3:>11.2f}  {pvt_10.rs_m3_per_m3:>11.2f}")
print(f"{'Mixture density ρ (kg/m³)':40s}  {pvt_80.mixture_density_kg_per_m3:>11.1f}  {pvt_10.mixture_density_kg_per_m3:>11.1f}")
print(f"{'Mixture viscosity μ (mPa·s)':40s}  {pvt_80.mixture_viscosity_pa_s*1e3:>11.3f}  {pvt_10.mixture_viscosity_pa_s*1e3:>11.3f}")
print(f"{'Mixture Cp (J/kg·K)':40s}  {pvt_80.mixture_specific_heat_j_per_kg_k:>11.0f}  {pvt_10.mixture_specific_heat_j_per_kg_k:>11.0f}")
print()
print(f"Viscosity ratio cold/hot = {pvt_10.mixture_viscosity_pa_s / pvt_80.mixture_viscosity_pa_s:.1f}×")
print("  → hot oil has far lower viscosity; same ΔP drives much higher flow")

# ── Network cases ─────────────────────────────────────────────────────────────
def _make_case(T_in, fluid, name):
    return NetworkCase(
        name             = name,
        fluid_model      = fluid,
        pressure_inlets  = (PressureBoundary(node_id=1, pressure_pa=P_IN),),
        pressure_outlets = (PressureBoundary(node_id=2, pressure_pa=P_OUT),),
        components       = (Pipe(
            component_id                        = "hot_oil_pipeline",
            start_node                          = 1,
            end_node                            = 2,
            diameter_m                          = D,
            length_m                            = L,
            absolute_roughness_m                = 46e-6,
            heat_transfer_coefficient_w_per_m2k = U,
            ambient_temperature_c               = T_AMB,
            n_thermal_segments                  = 20,
        ),),
        thermal_inlets   = (
            ThermalBoundary(node_id=1, temperature_c=T_in, bc_type="fixed_temperature"),
        ),
    )

# ── Solve both cases ──────────────────────────────────────────────────────────
solver     = SteadyBlackOilSolver()
res_hot    = solver.solve(_make_case(T_HOT, fluid_hot,  "hot-oil-20km"))
res_cold   = solver.solve(_make_case(T_AMB, fluid_cold, "cold-oil-20km"))

T_out_hot  = res_hot.node_temperatures_c[2]
mdot_hot   = res_hot.component_flows[0].mass_flow_kg_per_s
mdot_cold  = res_cold.component_flows[0].mass_flow_kg_per_s

# ── Analytical temperature verification ───────────────────────────────────────
pvt_mean   = fluid_hot.pvt(5.0e6, 0.5 * (T_HOT + T_out_hot))
Cp_mean    = pvt_mean.mixture_specific_heat_j_per_kg_k
NTU        = math.pi * D * U * L / (mdot_hot * Cp_mean)
T_out_ana  = T_AMB + (T_HOT - T_AMB) * math.exp(-NTU)
err_T      = abs(T_out_hot - T_out_ana)

print()
print("=" * 65)
print("HOT case (T_in = 80 °C):")
print(f"  Converged:       {res_hot.converged}")
print(f"  T_outlet:        {T_out_hot:.2f} °C")
print(f"  T_out analytical:{T_out_ana:.2f} °C  (NTU = {NTU:.4f})")
print(f"  Error:           {err_T:.3f} °C  (tolerance: ±2 °C)")
print(f"  Mass flow:       {mdot_hot:.2f} kg/s")

print()
print("COLD case (T_in = 10 °C, isothermal at ambient):")
print(f"  Converged:       {res_cold.converged}")
print(f"  Mass flow:       {mdot_cold:.2f} kg/s")

print()
flow_ratio = mdot_hot / mdot_cold
print(f"Flow ratio (hot / cold) = {flow_ratio:.3f}")
print(f"  → Hot oil delivers {(flow_ratio - 1.0) * 100:.1f}% more flow for the same ΔP")
print(f"    driven by viscosity: {pvt_80.mixture_viscosity_pa_s*1e3:.3f} mPa·s (inlet)")
print(f"    vs cold viscosity:   {pvt_10.mixture_viscosity_pa_s*1e3:.3f} mPa·s (10 °C, 2 MPa)")

# ── PVT change along the hot pipe ─────────────────────────────────────────────
pvt_outlet = fluid_hot.pvt(P_OUT, T_out_hot)
print()
print("PVT change along the hot pipe:")
print(f"  Inlet  (8 MPa, 80 °C): ρ = {pvt_80.mixture_density_kg_per_m3:.1f} kg/m³  "
      f"μ = {pvt_80.mixture_viscosity_pa_s*1e3:.3f} mPa·s  "
      f"Rs = {pvt_80.rs_m3_per_m3:.2f} m³/m³  (undersaturated: P > Pb)")
print(f"  Outlet (2 MPa, {T_out_hot:.1f} °C): ρ = {pvt_outlet.mixture_density_kg_per_m3:.1f} kg/m³  "
      f"μ = {pvt_outlet.mixture_viscosity_pa_s*1e3:.3f} mPa·s  "
      f"Rs = {pvt_outlet.rs_m3_per_m3:.2f} m³/m³  (two-phase: P < Pb)")
print(f"  Gas liberation: ΔRs = {pvt_80.rs_m3_per_m3 - pvt_outlet.rs_m3_per_m3:.2f} m³/m³  "
      f"({(pvt_80.rs_m3_per_m3 - pvt_outlet.rs_m3_per_m3)/pvt_80.rs_m3_per_m3*100:.0f}% of inlet GOR)")
