"""Tutorial 06 — Gas compressor station: shaft-work heating.

A natural-gas compressor station boosts methane from 30 bar to 60 bar.
Shaft work raises the gas temperature significantly.  This tutorial shows:

  1. Compression heats the gas by ΔT ≈ ΔP / (ρ · Cp) — much larger than
     for liquid pumps because gas density is ~100× lower.
  2. Real gas (Peng-Robinson EOS) heats up *more* than ideal gas at the same
     ΔP.  Although Z < 1 gives higher density and therefore less shaft work
     per unit mass (w = ΔP/ρ), the enthalpy departure at high pressure is
     more negative and stores energy in intermolecular forces, requiring a
     larger temperature rise to absorb that shaft work.  This is the
     complement of the JT effect (Tutorial 05): the same departure-enthalpy
     physics that causes real-gas cooling on expansion causes enhanced
     heating on compression.
  3. The hot compressed gas cools as it flows downstream.
  4. (Part B) A ``Compressor`` at the same pressure ratio: the work of the gas's
     own isentrope, ∫v·dP, and of a real machine with an efficiency.

Network layout (series):

    Source ──── Pipe 1 ──── Compressor ──── Pipe 2 ──── Sink
    (30 bar,       200 m,                     1000 m,   60 bar)
    20 °C       U=2 W/m²K                 U=2 W/m²K
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]  # …/angelica/angelica/src
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe, Pump
from angelica.properties.compressible_fluid import CompressibleFluid
from angelica.properties.eos import IdealGasEOS, PengRobinsonEOS
from angelica.solvers import SteadyCompressibleSolver

_M     = 0.016043   # kg/mol
_MU    = 1.1e-5     # Pa·s
_CP    = 2220.0     # J/(kg·K)
_K     = 0.033      # W/(m·K)
_ROUGH = 4.6e-5     # m

_D     = 0.30       # m
_U     = 2.0        # W/(m²·K)
_T_AMB = 10.0       # °C
_T_IN  = 20.0       # °C
_P_IN  = 30e5       # Pa
_P_OUT = 60e5       # Pa

# Compressor Q-H curve.  Q in m³/h (inlet conditions), H in m.
# Designed for ΔP ≈ 30 bar at Q ≈ 900 m³/h
# (ρ_methane_PR ≈ 21.6 kg/m³ at 30 bar, 20 °C → ṁ ≈ 5.4 kg/s).
_COMPRESSOR_CURVE = (
    (    0.0, 16000.0),
    (  400.0, 15500.0),
    (  800.0, 14500.0),
    ( 1200.0, 13000.0),
    ( 1600.0, 10000.0),
    ( 2000.0,  5000.0),
    ( 2200.0,     0.0),
)


def _build_case(eos) -> NetworkCase:
    fluid = CompressibleFluid.from_constants(
        eos=eos,
        viscosity_pa_s=_MU,
        specific_heat_j_per_kg_k=_CP,
        thermal_conductivity_w_per_m_k=_K,
        reference_pressure_pa=_P_IN,
        reference_temperature_c=_T_IN,
    )
    return NetworkCase(
        name="Gas compressor station — shaft-work heating",
        fluid_model=fluid,
        node_ids=(1, 2, 3, 4),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=_P_IN),),
        pressure_outlets=(PressureBoundary(node_id=4, pressure_pa=_P_OUT),),
        thermal_inlets=(
            ThermalBoundary(node_id=1, bc_type="fixed_temperature", temperature_c=_T_IN),
            ThermalBoundary(node_id=4, bc_type="zero_gradient"),
        ),
        components=(
            Pipe(
                start_node=1, end_node=2,
                diameter_m=_D, length_m=200.0,
                absolute_roughness_m=_ROUGH,
                heat_transfer_coefficient_w_per_m2k=_U,
                ambient_temperature_c=_T_AMB,
                n_thermal_segments=3,
                component_id="feed_pipe",
            ),
            Pump(
                start_node=2, end_node=3,
                diameter_m=_D,
                curve_points_q_head=_COMPRESSOR_CURVE,
                component_id="compressor",
            ),
            Pipe(
                start_node=3, end_node=4,
                diameter_m=_D, length_m=1000.0,
                absolute_roughness_m=_ROUGH,
                heat_transfer_coefficient_w_per_m2k=_U,
                ambient_temperature_c=_T_AMB,
                n_thermal_segments=10,
                component_id="export_pipe",
            ),
        ),
    )


def main() -> None:
    eos_pr    = PengRobinsonEOS(
        molecular_weight_kg_per_mol=_M,
        critical_temperature_k=190.6,
        critical_pressure_pa=4.60e6,
        acentric_factor=0.011,
    )
    eos_ideal = IdealGasEOS(molecular_weight_kg_per_mol=_M)

    result_pr   = SteadyCompressibleSolver().solve(_build_case(eos_pr))
    result_id   = SteadyCompressibleSolver().solve(_build_case(eos_ideal))
    case_pr     = _build_case(eos_pr)
    case_id     = _build_case(eos_ideal)

    print("Gas compressor station — compression heating")
    print(f"Source: {_P_IN/1e5:.0f} bar, {_T_IN:.0f} °C   Sink: {_P_OUT/1e5:.0f} bar\n")

    for label, result, case in (
        ("Peng-Robinson (real gas)", result_pr, case_pr),
        ("Ideal gas              ", result_id, case_id),
    ):
        T1 = result.node_temperatures_c[1]
        T2 = result.node_temperatures_c[2]
        T3 = result.node_temperatures_c[3]
        T4 = result.node_temperatures_c[4]
        P1 = result.node_pressures_pa[1]
        P2 = result.node_pressures_pa[2]
        P3 = result.node_pressures_pa[3]
        P4 = result.node_pressures_pa[4]

        dT_compr = T3 - T2
        dP_compr = P3 - P2

        fluid   = case.fluid_model
        # The solver uses average pressure for density, matching energy.py
        rho_avg = fluid.eos.density(0.5 * (P2 + P3), T2)
        h_in    = fluid.enthalpy_j_per_kg(P2, T2)
        h_out   = fluid.enthalpy_j_per_kg(P3, T3)
        dh      = h_out - h_in
        w_shaft = dP_compr / rho_avg

        mdot = next(
            f.mass_flow_kg_per_s
            for f in result.component_flows
            if f.label and "compressor" in f.label.lower()
        )

        print(f"── {label} ──")
        print(f"  Converged: {result.converged}")
        print(f"  Node 1 (source):           {P1/1e5:6.2f} bar   {T1:7.3f} °C")
        print(f"  Node 2 (upstream compr.):  {P2/1e5:6.2f} bar   {T2:7.3f} °C")
        print(f"  Node 3 (downstream compr.):{P3/1e5:6.2f} bar   {T3:7.3f} °C")
        print(f"  Node 4 (export sink):      {P4/1e5:6.2f} bar   {T4:7.3f} °C")
        print()
        print(f"  ṁ:                          {mdot:.3f} kg/s")
        print(f"  ΔP compressor:              {dP_compr/1e5:+.3f} bar")
        print(f"  ΔT compressor:              {dT_compr:+.3f} °C")
        print(f"  ρ_avg (solver):             {rho_avg:.3f} kg/m³")
        print(f"  Shaft work w = ΔP/ρ_avg:   {w_shaft:,.0f} J/kg")
        print(f"  Δh = h_out − h_in:         {dh:,.0f} J/kg  (≈ w by construction)")
        print()

    dT_pr = result_pr.node_temperatures_c[3] - result_pr.node_temperatures_c[2]
    dT_id = result_id.node_temperatures_c[3] - result_id.node_temperatures_c[2]
    print("Compression heating comparison:")
    print(f"  Real gas (PR):  ΔT = {dT_pr:+.2f} °C")
    print(f"  Ideal gas:      ΔT = {dT_id:+.2f} °C")
    print(f"  Difference:     {dT_pr - dT_id:+.2f} °C  (PR heats more — departure enthalpy at high P)")

    # ── Part B: the same pressure ratio through a Compressor ─────────────────
    # The head-curve "compressor" above gives the gas w = ΔP/ρ (one density):
    # the work of a reversible machine on a liquid.  A Compressor follows the
    # gas's own isentrope (dh = v·dP) and its efficiency.
    from dataclasses import replace

    from angelica.core.components import Compressor

    ratio = result_pr.node_pressures_pa[3] / result_pr.node_pressures_pa[2]
    print()
    print(f"Part B — a Compressor at the same pressure ratio ({ratio:.4f}), Peng-Robinson")
    print(f"  {'Machine':34} {'ṁ (kg/s)':>9} {'T out (°C)':>11} {'Power (MW)':>11}")
    base = _build_case(eos_pr)
    fl = base.fluid_model
    T2, T3 = result_pr.node_temperatures_c[2], result_pr.node_temperatures_c[3]
    w_curve = fl.enthalpy_j_per_kg(result_pr.node_pressures_pa[3], T3) - \
        fl.enthalpy_j_per_kg(result_pr.node_pressures_pa[2], T2)
    works = {}
    for label, eta in (("isentropic, η = 1 (reversible)", 1.0), ("isentropic, η = 0.75", 0.75),
                       ("polytropic, η = 0.75", 0.75)):
        kind = "polytropic" if label.startswith("polytropic") else "isentropic"
        comps = tuple(Compressor(start_node=2, end_node=3, diameter_m=_D, pressure_ratio=ratio, efficiency=eta,
                                 efficiency_type=kind, component_id="compressor")
                      if c.component_id == "compressor" else c for c in base.components)
        r = SteadyCompressibleSolver().solve(replace(base, components=comps))
        c = r.compressors[0]
        works[label] = (c.power_w / c.mass_flow_kg_per_s, c.outlet_temperature_c)
        print(f"  {label:34} {c.mass_flow_kg_per_s:9.3f} {c.outlet_temperature_c:11.2f} {c.power_w / 1e6:11.3f}"
              + ("" if r.converged else "  (not converged)"))
    w_rev, t_rev = works["isentropic, η = 1 (reversible)"]
    w_75 = works["isentropic, η = 0.75"][0]
    print(f"  • Per kg of gas the head curve gives w = ΔP/ρ = {w_curve / 1e3:.1f} kJ/kg, the")
    print(f"    reversible compressor ∫v·dP = {w_rev / 1e3:.1f} kJ/kg ({(w_rev / w_curve - 1) * 100:+.1f} %): ρ at the")
    print("    mean pressure but the inlet temperature misses the heating along the")
    print(f"    compression, and the gas leaves {t_rev - T3:.1f} K hotter than the curve says.")
    print(f"  • A real machine (η = 0.75) needs {(w_75 / w_rev - 1) * 100:.0f} % more work per kg; the same number")
    print("    as a polytropic efficiency means a less efficient machine, hotter still.")


if __name__ == "__main__":
    main()
