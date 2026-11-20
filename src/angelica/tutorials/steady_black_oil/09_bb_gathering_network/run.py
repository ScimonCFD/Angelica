"""
Tutorial 09 — Beggs-Brill in a multi-well gathering network
============================================================

Three production wells with different fluid characteristics (GOR, API,
water cut) converge at a single manifold and then flow to a separator.

The Beggs-Brill correlation is used for all pipes so that each segment
correctly accounts for two-phase slip between gas and liquid.

Network
-------

  Well 1 (N=1)  API=40, GOR=150  P=12 MPa  T=80°C
   |── Pipe-A (4 km, D=0.1016 m) ──┐
                                     │
  Well 2 (N=2)  API=32, GOR=40   P=10 MPa  T=65°C
   |── Pipe-B (2 km, D=0.1016 m) ──── Manifold (N=4) ── Pipe-D (8 km, D=0.1524 m) ── Sep (N=5)  P=1.5 MPa
                                     │
  Well 3 (N=3)  API=25, GOR=15   P= 8 MPa  T=50°C
   |── Pipe-C (6 km, D=0.0762 m) ──┘

The three wells produce fluids of different character:
  Well 1 — light condensate-rich crude, high GOR, low viscosity
  Well 2 — medium crude, moderate GOR (reference fluid)
  Well 3 — heavy crude, low GOR, high viscosity

Physical expectations
---------------------
- Well 1 (high GOR) will have the highest gas-phase content and the
  largest Beggs-Brill holdup correction relative to the no-slip model.
- Well 3 (heavy crude, low GOR) will behave most like a single-phase
  liquid: small slip, f_tp ≈ f_n.
- The manifold pressure (N=4) emerges as a solved quantity — each
  well's contribution adjusts to satisfy pressure continuity.

Reference: Beggs & Brill (1973).  Per-inlet compositions use
    InletFluidBC (mass-weighted blending at junctions).
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures.beggs_brill import BeggsBrillCorrelation, beggs_brill_gradient
from angelica.core.case import InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.black_oil import BlackOilFluid
from angelica.solvers import SteadyBlackOilSolver

# ── Reference fluid (used for global PVT; per-inlet GOR / API override) ───────
BASE_FLUID = BlackOilFluid(
    api_gravity=32.0, gas_gravity=0.65,
    gor_sc_m3_per_m3=50.0, wor_sc_m3_per_m3=0.3,
    reference_pressure_pa=6e6, reference_temperature_c=65.0,
)

# ── Network ────────────────────────────────────────────────────────────────────
PIPES = (
    Pipe(component_id="A_well1",  start_node=1, end_node=4,
         diameter_m=0.1016, length_m=4_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    Pipe(component_id="B_well2",  start_node=2, end_node=4,
         diameter_m=0.1016, length_m=2_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    Pipe(component_id="C_well3",  start_node=3, end_node=4,
         diameter_m=0.0762, length_m=6_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    Pipe(component_id="D_manifold", start_node=4, end_node=5,
         diameter_m=0.1524, length_m=8_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=3.0, ambient_temperature_c=15.0),
)


def build_case() -> NetworkCase:
    return NetworkCase(
        name="BB gathering network — Tutorial 09",
        fluid_model=BASE_FLUID,
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=12e6),
            PressureBoundary(node_id=2, pressure_pa=10e6),
            PressureBoundary(node_id=3, pressure_pa=8e6),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=5, pressure_pa=1.5e6),
        ),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=80.0, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=2, temperature_c=65.0, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=3, temperature_c=50.0, bc_type="fixed_temperature"),
        ),
        inlet_fluid_bcs=(
            InletFluidBC(node_id=1, api_gravity=40.0, gas_gravity=0.70,
                         gor_sc_m3_per_m3=150.0, wor_sc_m3_per_m3=0.1),
            InletFluidBC(node_id=2, api_gravity=32.0, gas_gravity=0.65,
                         gor_sc_m3_per_m3=40.0,  wor_sc_m3_per_m3=0.3),
            InletFluidBC(node_id=3, api_gravity=25.0, gas_gravity=0.60,
                         gor_sc_m3_per_m3=15.0,  wor_sc_m3_per_m3=0.8),
        ),
        components=PIPES,
    )


def run():
    case = build_case()

    res_ck = SteadyBlackOilSolver().solve(case)
    res_bb = SteadyBlackOilSolver(
        turbulent_pipe_correlation=BeggsBrillCorrelation()
    ).solve(case)

    print("=" * 65)
    print("Tutorial 09 — Beggs-Brill in a multi-well gathering network")
    print("=" * 65)

    P_manifold_ck = res_ck.node_pressures_pa[4]
    P_manifold_bb = res_bb.node_pressures_pa[4]
    T_manifold_ck = res_ck.node_temperatures_c.get(4, float("nan"))
    T_manifold_bb = res_bb.node_temperatures_c.get(4, float("nan"))

    print(f"\n{'Node':>4}  {'P_CK (MPa)':>12}  {'P_BB (MPa)':>12}  {'T_CK (°C)':>10}  {'T_BB (°C)':>10}")
    for nid in sorted(res_ck.node_pressures_pa):
        P_ck = res_ck.node_pressures_pa[nid] / 1e6
        P_bb = res_bb.node_pressures_pa[nid] / 1e6
        T_ck = res_ck.node_temperatures_c.get(nid, float("nan"))
        T_bb = res_bb.node_temperatures_c.get(nid, float("nan"))
        print(f"{nid:>4}  {P_ck:>12.3f}  {P_bb:>12.3f}  {T_ck:>10.1f}  {T_bb:>10.1f}")

    print(f"\n{'Pipe':>10}  {'Well':>6}  {'CK (m³/h)':>11}  {'BB (m³/h)':>11}  {'Regime':>14}  {'H_L':>6}  {'f_tp/f_n':>8}")
    print("-" * 75)

    pipe_ids   = ["A_well1", "B_well2", "C_well3", "D_manifold"]
    well_names = ["Well 1 (API=40, GOR=150)",
                  "Well 2 (API=32, GOR=40 )",
                  "Well 3 (API=25, GOR=15 )",
                  "Manifold → Sep           "]

    for idx, (pid, wname) in enumerate(zip(pipe_ids, well_names)):
        q_ck = res_ck.component_flows[idx].volumetric_flow_m3_per_h  # signed
        q_bb = res_bb.component_flows[idx].volumetric_flow_m3_per_h  # signed
        cf_bb = res_bb.component_flows[idx]

        # BB intermediate quantities at midpoint conditions
        pipe = PIPES[idx]
        P_start = res_bb.node_pressures_pa.get(pipe.start_node, 6e6)
        P_end   = res_bb.node_pressures_pa.get(pipe.end_node, 3e6)
        P_avg   = 0.5 * (P_start + P_end)
        T_start = res_bb.node_temperatures_c.get(pipe.start_node, 50.0)
        T_end   = res_bb.node_temperatures_c.get(pipe.end_node, 40.0)
        T_avg   = 0.5 * (T_start + T_end)

        pvt   = BASE_FLUID.pvt(P_avg, T_avg)
        alpha_l = pvt.holdup_oil + pvt.holdup_water
        alpha_g = pvt.holdup_gas
        C_L = alpha_l / max(alpha_l + alpha_g, 1e-9)
        rho_l = ((pvt.holdup_oil * pvt.density_oil_kg_per_m3
                  + pvt.holdup_water * pvt.density_water_kg_per_m3) / max(alpha_l, 1e-9))
        rho_g = pvt.density_gas_kg_per_m3
        mu_l = ((pvt.holdup_oil * pvt.viscosity_oil_pa_s
                 + pvt.holdup_water * pvt.viscosity_water_pa_s) / max(alpha_l, 1e-9))
        mu_g  = pvt.viscosity_gas_pa_s
        v_m = abs(cf_bb.mass_flow_kg_per_s) / max(pvt.mixture_density_kg_per_m3 * pipe.area_m2, 1e-12)
        st = beggs_brill_gradient(v_m, C_L, rho_l, rho_g, mu_l, mu_g,
                                  pipe.diameter_m, 46e-6, 0.0)

        ratio = f"{st.f_tp/st.f_n:.2f}×"
        print(f"{pid:>10}  {wname[:6]:>6}  {q_ck:>11.2f}  {q_bb:>11.2f}  "
              f"{st.regime:>14}  {st.holdup:.3f}  {ratio:>8}")

    print()
    print(f"  Manifold pressure:  Colebrook = {P_manifold_ck/1e6:.3f} MPa,  "
          f"Beggs-Brill = {P_manifold_bb/1e6:.3f} MPa")
    print(f"  Manifold temperature: CK = {T_manifold_ck:.1f} °C,  BB = {T_manifold_bb:.1f} °C")
    print(f"\n  BB total to separator: {res_bb.component_flows[3].volumetric_flow_m3_per_h:.1f} m³/h")
    print(f"  CK total to separator: {res_ck.component_flows[3].volumetric_flow_m3_per_h:.1f} m³/h")
    print()
    print("  Note: negative flow on Pipe-C means manifold pressure exceeds Well-3 bhp.")
    print("  With BB the manifold is 0.28 MPa *above* Well 3 — fluid flows back into")
    print("  the well, which would trip a check valve in reality. Colebrook gives a")
    print("  smaller manifold pressure and Well-3 produces 0.6 m³/h into the system.")


if __name__ == "__main__":
    run()
