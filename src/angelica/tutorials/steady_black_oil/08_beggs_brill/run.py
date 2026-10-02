"""
Tutorial 08 — Beggs-Brill two-phase flow correlation
=====================================================

This tutorial compares two pressure-drop models on the same horizontal
black-oil pipeline:

  A. Default homogeneous (no-slip) Colebrook-White — treats the gas-oil-water
     mixture as a single fluid with mixture properties.
  B. Beggs-Brill (1973) — accounts for slip between gas and liquid, computes
     the actual in-situ liquid holdup H_L, and corrects both the gravity and
     friction gradients.

For a horizontal pipe the gravity term is zero and the main difference is
in the friction term (two-phase friction factor f_tp vs no-slip f_n).
Beggs-Brill typically predicts a higher friction drop for the same boundary
conditions because:
  • H_L > C_L (liquid accumulates due to slip) → two-phase correction y > 1
  • f_tp = f_n · exp(s) ≥ f_n

Expected results
----------------
Pipe:      D = 0.1016 m (4 in), L = 10 km, roughness = 46 µm, horizontal
Fluid:     32 API, GOR = 80 m³/m³ (high gas), WOR = 0.2
Boundary:  P_in = 6 MPa (fixed), P_out = 2 MPa (fixed)

Both solvers converge on the same pressure boundary conditions, so the flow
rate differs:
  • Colebrook:   higher flow rate (underestimates frictional resistance)
  • Beggs-Brill: lower flow rate (correctly accounts for slip-enhanced friction)

The ratio  q_colebrook / q_BB  illustrates the error from ignoring slip.

Reference: Beggs, H.D. and Brill, J.P. (1973). A Study of Two-Phase Flow in
Inclined Pipes. JPT, May, 607-617.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures.beggs_brill import BeggsBrillCorrelation, beggs_brill_gradient
from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.black_oil import BlackOilFluid
from angelica.solvers import SteadyBlackOilSolver

# ── Fluid ──────────────────────────────────────────────────────────────────
FLUID = BlackOilFluid(
    api_gravity=32.0,
    gas_gravity=0.65,
    gor_sc_m3_per_m3=80.0,     # high GOR → significant gas phase at line conditions
    wor_sc_m3_per_m3=0.2,
    reference_pressure_pa=4e6,
    reference_temperature_c=50.0,
)

# ── Pipe geometry ──────────────────────────────────────────────────────────
D        = 0.1016   # m  (4 inch)
L        = 10_000.0 # m
ROUGH    = 46e-6    # m  (commercial steel)

# ── Network case helper ────────────────────────────────────────────────────

def build_case() -> NetworkCase:
    return NetworkCase(
        name="beggs-brill-tutorial-08",
        fluid_model=FLUID,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=6e6),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=2e6),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=50.0,
                                        bc_type="fixed_temperature"),),
        components=(Pipe(
            start_node=1, end_node=2,
            diameter_m=D, length_m=L,
            absolute_roughness_m=ROUGH,
            height_change_m=0.0,    # horizontal
        ),),
    )


def run():
    case = build_case()

    # ── A. Colebrook (default, no-slip homogeneous) ────────────────────────
    res_ck = SteadyBlackOilSolver().solve(case)

    # ── B. Beggs-Brill ─────────────────────────────────────────────────────
    res_bb = SteadyBlackOilSolver(
        turbulent_pipe_correlation=BeggsBrillCorrelation()
    ).solve(case)

    # ── Extract results ────────────────────────────────────────────────────
    q_ck = abs(res_ck.component_flows[0].volumetric_flow_m3_per_h)
    q_bb = abs(res_bb.component_flows[0].volumetric_flow_m3_per_h)

    mdot_ck = res_ck.component_flows[0].mass_flow_kg_per_s
    mdot_bb = res_bb.component_flows[0].mass_flow_kg_per_s

    # Phase fractions at average conditions (midpoint pressure)
    P_avg = 4e6  # Pa (midpoint of 6 MPa–2 MPa)
    T_avg = FLUID.reference_temperature_c
    pvt = FLUID.pvt(P_avg, T_avg)
    alpha_l = pvt.holdup_oil + pvt.holdup_water
    alpha_g = pvt.holdup_gas
    C_L = alpha_l / max(alpha_l + alpha_g, 1e-9)

    print("=" * 60)
    print("Tutorial 08 — Beggs-Brill two-phase flow")
    print("=" * 60)
    print(f"\nPipe: D={D*1000:.0f} mm, L={L/1e3:.0f} km, horizontal")
    print(f"P_in = 6.0 MPa  →  P_out = 2.0 MPa  (ΔP = 4.0 MPa)")
    print(f"\nFluid at avg conditions ({P_avg/1e6:.1f} MPa, {T_avg:.0f} °C):")
    print(f"  No-slip liquid content  C_L = {C_L:.3f}")
    print(f"  No-slip gas content     C_G = {1-C_L:.3f}")
    print(f"  Mix density             ρ_m = {pvt.mixture_density_kg_per_m3:.1f} kg/m³")
    print(f"  Mix viscosity           µ_m = {pvt.mixture_viscosity_pa_s*1000:.2f} mPa·s")

    # Holdup at representative velocity (from BB solve)
    v_bb = abs(res_bb.component_flows[0].mass_flow_kg_per_s) / (
        pvt.mixture_density_kg_per_m3 * 3.14159 * (D / 2) ** 2
    )
    st = beggs_brill_gradient(
        v_bb, C_L,
        (pvt.holdup_oil * pvt.density_oil_kg_per_m3 + pvt.holdup_water * pvt.density_water_kg_per_m3) / max(alpha_l, 1e-9),
        pvt.density_gas_kg_per_m3,
        (pvt.holdup_oil * pvt.viscosity_oil_pa_s + pvt.holdup_water * pvt.viscosity_water_pa_s) / max(alpha_l, 1e-9),
        pvt.viscosity_gas_pa_s,
        D, ROUGH, 0.0,
    )

    print(f"\nBeggs-Brill at v_m = {v_bb:.2f} m/s:")
    print(f"  Flow regime        = {st.regime}")
    print(f"  In-situ holdup H_L = {st.holdup:.3f}  (vs no-slip C_L = {C_L:.3f})")
    print(f"  Froude number N_Fr = {st.N_Fr:.2f}")
    print(f"  No-slip f_n        = {st.f_n:.5f}")
    print(f"  Two-phase f_tp     = {st.f_tp:.5f}  (ratio = {st.f_tp/st.f_n:.2f}×)")
    print(f"  dP/dz friction     = {st.dP_dz_friction:.0f} Pa/m")

    print(f"\n{'':>22} {'Colebrook':>12} {'Beggs-Brill':>12}")
    print(f"  {'Converged':20} {str(res_ck.converged):>12} {str(res_bb.converged):>12}")
    print(f"  {'Mass flow (kg/s)':20} {mdot_ck:>12.3f} {mdot_bb:>12.3f}")
    print(f"  {'Volume flow (m³/h)':20} {q_ck:>12.2f} {q_bb:>12.2f}")
    print(f"\n  Colebrook / BB ratio: {q_ck/q_bb:.3f}  (BB predicts lower flow)")
    print(f"  → ignoring slip overestimates flow by {(q_ck/q_bb - 1)*100:.1f} %")


if __name__ == "__main__":
    run()
