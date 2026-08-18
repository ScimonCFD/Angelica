from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[3] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.compressible_fluid import CompressibleFluid
from angelica.properties.eos import IdealGasEOS
from angelica.solvers import SteadyCompressibleSolver

_M_METHANE = 0.016043   # kg/mol
_MU_METHANE = 1.1e-5    # Pa·s
_CP_METHANE = 2_220.0   # J/(kg·K)
_K_METHANE = 0.033      # W/(m·K)
_ROUGHNESS = 4.6e-5     # m
_U = 5.0                # W/(m²·K) — overall heat transfer coefficient
_T_AMB = 10.0           # °C


def _pipe(start, end, *, d, length, label=None):
    return Pipe(
        start_node=start,
        end_node=end,
        diameter_m=d,
        length_m=length,
        absolute_roughness_m=_ROUGHNESS,
        heat_transfer_coefficient_w_per_m2k=_U,
        ambient_temperature_c=_T_AMB,
        component_id=label,
    )


def build_case() -> NetworkCase:
    fluid = CompressibleFluid.from_constants(
        eos=IdealGasEOS(molecular_weight_kg_per_mol=_M_METHANE),
        viscosity_pa_s=_MU_METHANE,
        specific_heat_j_per_kg_k=_CP_METHANE,
        thermal_conductivity_w_per_m_k=_K_METHANE,
        reference_pressure_pa=500_000.0,
        reference_temperature_c=40.0,
    )
    return NetworkCase(
        name="Looped Gas Pipeline with Heat Loss",
        fluid_model=fluid,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=700_000.0),),
        pressure_outlets=(PressureBoundary(node_id=6, pressure_pa=500_000.0),),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=40.0, bc_type="fixed_temperature"),
        ),
        components=(
            _pipe(1, 2, d=0.20, length=600.0, label="L1-feeder"),
            _pipe(2, 3, d=0.15, length=500.0, label="L2-upper"),
            _pipe(2, 4, d=0.12, length=400.0, label="L3-lower"),
            _pipe(3, 4, d=0.10, length=350.0, label="L4-cross"),
            _pipe(3, 5, d=0.12, length=400.0, label="L5-upper"),
            _pipe(4, 5, d=0.15, length=500.0, label="L6-lower"),
            _pipe(5, 6, d=0.20, length=600.0, label="L7-collector"),
        ),
        node_ids=(1, 2, 3, 4, 5, 6),
    )


def main() -> None:
    case = build_case()
    result = SteadyCompressibleSolver().solve(case)

    print(f"Case: {case.name}")
    print(f"Converged:           {result.converged}")
    print(f"Density iterations:  {len(result.density_history)}")
    print()

    print("Node pressures and temperatures:")
    labels = {1: "Source", 2: "Junction A", 3: "Upper jct",
              4: "Lower jct", 5: "Junction B", 6: "Sink"}
    for node_id in (1, 2, 3, 4, 5, 6):
        p = result.node_pressures_pa[node_id]
        t = result.node_temperatures_c[node_id]
        name = labels[node_id]
        print(f"  Node {node_id} ({name:<12}) {p/1e3:6.1f} kPa   {t:5.1f} °C")
    print()

    print("Component flows:")
    for flow in result.component_flows:
        print(f"  {flow.label:<22}  {flow.mass_flow_kg_per_s:6.3f} kg/s")
    print()

    mdot_feeder = result.component_flows[0].mass_flow_kg_per_s
    mdot_collector = result.component_flows[6].mass_flow_kg_per_s
    print("Mass balance:")
    print(f"  Feeder (L1):    {mdot_feeder:.4f} kg/s")
    print(f"  Collector (L7): {mdot_collector:.4f} kg/s")
    print(f"  Imbalance:      {abs(mdot_feeder - mdot_collector):.2e} kg/s")


if __name__ == "__main__":
    main()
