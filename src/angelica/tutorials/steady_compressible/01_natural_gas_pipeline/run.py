from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[3] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases import build_natural_gas_pipeline_case
from angelica.properties.eos import IdealGasEOS
from angelica.solvers import SteadyCompressibleSolver

_M_METHANE = 0.016043   # kg/mol
_R = 8.314              # J/(mol·K)
_T_C = 15.0             # °C — reference temperature


def rho_ideal(pressure_pa: float) -> float:
    return pressure_pa * _M_METHANE / (_R * (_T_C + 273.15))


def main() -> None:
    case = build_natural_gas_pipeline_case()
    solver = SteadyCompressibleSolver()
    result = solver.solve(case)

    print(f"Case: {case.name}")
    print(f"Converged:            {result.converged}")
    print(f"Density iterations:   {len(result.density_history)}")
    print()

    eos = IdealGasEOS(molecular_weight_kg_per_mol=_M_METHANE)
    p_in  = result.node_pressures_pa[1]
    p_jct = result.node_pressures_pa[2]
    p_out_a = result.node_pressures_pa[3]
    p_out_b = result.node_pressures_pa[4]

    print("Node pressures:")
    print(f"  Node 1 (inlet):      {p_in/1e3:8.2f} kPa   rho = {eos.density(p_in,  _T_C):.3f} kg/m³")
    print(f"  Node 2 (junction):   {p_jct/1e3:8.2f} kPa   rho = {eos.density(p_jct, _T_C):.3f} kg/m³")
    print(f"  Node 3 (outlet A):   {p_out_a/1e3:8.2f} kPa   rho = {eos.density(p_out_a, _T_C):.3f} kg/m³")
    print(f"  Node 4 (outlet B):   {p_out_b/1e3:8.2f} kPa   rho = {eos.density(p_out_b, _T_C):.3f} kg/m³")
    print()

    print("Component flows:")
    for flow in result.component_flows:
        print(f"  {flow.label:<20}  {flow.mass_flow_kg_per_s:8.4f} kg/s  ({flow.volumetric_flow_m3_per_h:.2f} m³/h)")
    print()

    print("Mass balance (inlet = outlet A + outlet B):")
    mdot_1 = result.component_flows[0].mass_flow_kg_per_s
    mdot_2 = result.component_flows[1].mass_flow_kg_per_s
    mdot_3 = result.component_flows[2].mass_flow_kg_per_s
    print(f"  ṁ_pipe1 = {mdot_1:.4f} kg/s")
    print(f"  ṁ_pipe2 + ṁ_pipe3 = {mdot_2 + mdot_3:.4f} kg/s")
    print(f"  Imbalance: {abs(mdot_1 - mdot_2 - mdot_3):.2e} kg/s")
    print()

    print("Compressibility effect:")
    rho_in_node  = eos.density(p_in, _T_C)
    rho_out_node = eos.density(p_out_a, _T_C)
    print(f"  Inlet density:  {rho_in_node:.3f} kg/m³")
    print(f"  Outlet density: {rho_out_node:.3f} kg/m³")
    print(f"  Ratio (rho_in / rho_out): {rho_in_node / rho_out_node:.3f}  (= P_in / P_out for ideal gas)")


if __name__ == "__main__":
    main()
