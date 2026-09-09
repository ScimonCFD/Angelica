from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Fitting, Pipe
from angelica.properties.compressible_fluid import CompressibleFluid
from angelica.properties.eos import IdealGasEOS, PengRobinsonEOS
from angelica.solvers import SteadyCompressibleSolver

# ── Methane properties ─────────────────────────────────────────────────────────
_M = 0.016043      # kg/mol
_MU = 1.1e-5       # Pa·s
_CP = 2_220.0      # J/(kg·K)
_K = 0.033         # W/(m·K)
_ROUGHNESS = 4.6e-5


def _build_case(eos) -> NetworkCase:
    fluid = CompressibleFluid.from_constants(
        eos=eos,
        viscosity_pa_s=_MU,
        specific_heat_j_per_kg_k=_CP,
        thermal_conductivity_w_per_m_k=_K,
        reference_pressure_pa=60e5,
        reference_temperature_c=40.0,
    )
    return NetworkCase(
        name="Valve JT effect",
        fluid_model=fluid,
        node_ids=(1, 2, 3, 4),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=70e5),),
        pressure_outlets=(PressureBoundary(node_id=4, pressure_pa=50e5),),
        thermal_inlets=(
            ThermalBoundary(node_id=1, bc_type="fixed_temperature", temperature_c=40.0),
            ThermalBoundary(node_id=4, bc_type="zero_gradient"),
        ),
        components=(
            Pipe(
                start_node=1, end_node=2, diameter_m=0.15, length_m=500.0,
                absolute_roughness_m=_ROUGHNESS,
                heat_transfer_coefficient_w_per_m2k=5.0,
                ambient_temperature_c=10.0, n_thermal_segments=5,
                component_id="upstream_pipe",
            ),
            Fitting(
                start_node=2, end_node=3, diameter_m=0.15,
                loss_coefficient=50.0, component_id="valve",
            ),
            Pipe(
                start_node=3, end_node=4, diameter_m=0.15, length_m=500.0,
                absolute_roughness_m=_ROUGHNESS,
                heat_transfer_coefficient_w_per_m2k=5.0,
                ambient_temperature_c=10.0, n_thermal_segments=5,
                component_id="downstream_pipe",
            ),
        ),
    )


def main() -> None:
    eos_pr = PengRobinsonEOS(
        molecular_weight_kg_per_mol=_M,
        critical_temperature_k=190.6,
        critical_pressure_pa=4.60e6,
        acentric_factor=0.011,
    )
    eos_ideal = IdealGasEOS(molecular_weight_kg_per_mol=_M)

    result_pr = SteadyCompressibleSolver().solve(_build_case(eos_pr))
    result_id = SteadyCompressibleSolver().solve(_build_case(eos_ideal))

    case_pr = _build_case(eos_pr)
    case_id = _build_case(eos_ideal)

    for label, result, case in (
        ("Peng-Robinson (real gas)", result_pr, case_pr),
        ("Ideal gas              ", result_id, case_id),
    ):
        P2 = result.node_pressures_pa[2]
        P3 = result.node_pressures_pa[3]
        T2 = result.node_temperatures_c[2]
        T3 = result.node_temperatures_c[3]

        dh_str = ""
        fluid = case.fluid_model
        if hasattr(fluid, "enthalpy_j_per_kg"):
            h_in = fluid.enthalpy_j_per_kg(P2, T2)
            h_out = fluid.enthalpy_j_per_kg(P3, T3)
            dh_str = f"   Δh = {h_out - h_in:+.3f} J/kg"

        print(f"\n── {label} ──")
        print(f"  Converged: {result.converged}")
        print(f"  Node 2 (upstream valve):   {P2/1e5:.3f} bar  {T2:.3f} °C")
        print(f"  Node 3 (downstream valve): {P3/1e5:.3f} bar  {T3:.3f} °C")
        print(f"  ΔP valve: {(P2-P3)/1e5:.3f} bar")
        print(f"  ΔT valve: {T3-T2:+.4f} °C{dh_str}")

    print()
    dT_pr = result_pr.node_temperatures_c[3] - result_pr.node_temperatures_c[2]
    dT_id = result_id.node_temperatures_c[3] - result_id.node_temperatures_c[2]
    print("Joule-Thomson cooling:")
    print(f"  Real gas (PR):  ΔT = {dT_pr:+.4f} °C")
    print(f"  Ideal gas:      ΔT = {dT_id:+.4f} °C  (zero by definition)")


if __name__ == "__main__":
    main()
