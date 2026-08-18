from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[3] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import FlowBoundary, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.compressible_fluid import CompressibleFluid
from angelica.properties.eos import IdealGasEOS
from angelica.solvers import SteadyCompressibleSolver

_M_METHANE = 0.016043   # kg/mol
_MU_METHANE = 1.1e-5    # Pa·s
_CP_METHANE = 2_220.0   # J/(kg·K)
_K_METHANE = 0.033      # W/(m·K)
_ROUGHNESS = 4.6e-5     # m — commercial steel
_T_C = 20.0             # °C


def build_case() -> NetworkCase:
    fluid = CompressibleFluid.from_constants(
        eos=IdealGasEOS(molecular_weight_kg_per_mol=_M_METHANE),
        viscosity_pa_s=_MU_METHANE,
        specific_heat_j_per_kg_k=_CP_METHANE,
        thermal_conductivity_w_per_m_k=_K_METHANE,
        reference_pressure_pa=800_000.0,
        reference_temperature_c=_T_C,
    )
    return NetworkCase(
        name="Natural Gas Pipeline — Flow BC Cross-Validation",
        fluid_model=fluid,
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=800_000.0),
        ),
        pressure_outlets=(),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=_T_C, bc_type="fixed_temperature"),
        ),
        flow_outlets=(
            FlowBoundary(node_id=3, mass_flow_kg_per_s=1.255479),
            FlowBoundary(node_id=4, mass_flow_kg_per_s=1.540697),
        ),
        components=(
            Pipe(1, 2, diameter_m=0.15, length_m=500.0,
                 absolute_roughness_m=_ROUGHNESS, component_id="pipe_1"),
            Pipe(2, 3, diameter_m=0.10, length_m=300.0,
                 absolute_roughness_m=_ROUGHNESS, component_id="pipe_2"),
            Pipe(2, 4, diameter_m=0.10, length_m=200.0,
                 absolute_roughness_m=_ROUGHNESS, component_id="pipe_3"),
        ),
        node_ids=(1, 2, 3, 4),
    )


def main() -> None:
    case = build_case()
    result = SteadyCompressibleSolver().solve(case)

    print(f"Case: {case.name}")
    print(f"Converged:            {result.converged}")
    print(f"Density iterations:   {len(result.density_history)}")
    print()

    eos = IdealGasEOS(molecular_weight_kg_per_mol=_M_METHANE)
    p = result.node_pressures_pa

    print("Node pressures:")
    for node_id in (1, 2, 3, 4):
        rho = eos.density(p[node_id], _T_C)
        print(f"  Node {node_id}: {p[node_id]/1e3:8.2f} kPa   rho = {rho:.3f} kg/m³")
    print()

    print("Component flows:")
    for flow in result.component_flows:
        print(f"  {flow.label:<20}  {flow.mass_flow_kg_per_s:8.4f} kg/s")
    print()

    mdot_in = result.component_flows[0].mass_flow_kg_per_s
    mdot_3 = result.component_flows[1].mass_flow_kg_per_s
    mdot_4 = result.component_flows[2].mass_flow_kg_per_s
    print("Mass balance:")
    print(f"  ṁ_pipe1 = {mdot_in:.4f} kg/s  (inlet)")
    print(f"  ṁ_pipe2 + ṁ_pipe3 = {mdot_3 + mdot_4:.4f} kg/s  (outlets)")
    print(f"  Imbalance: {abs(mdot_in - mdot_3 - mdot_4):.2e} kg/s")
    print()

    print("Cross-validation: outlet pressures with flow BCs")
    print(f"  Node 3: {p[3]/1e3:.2f} kPa  (should match pressure-BC case ≈ 500 kPa)")
    print(f"  Node 4: {p[4]/1e3:.2f} kPa  (should match pressure-BC case ≈ 500 kPa)")


if __name__ == "__main__":
    main()
