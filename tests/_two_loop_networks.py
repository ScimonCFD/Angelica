"""
A gathering network with two independent loops, built for every fluid mode
(``tests/test_two_loop_networks.py``).

  Inlet A (1) ── feed A ──► N2 ───── 2-3 ─────► N3 ── 3-5 ──► N5 ── export ──► Outlet (6)
                             │                   │            ▲
                            2-4                 3-4           │
                             ▼                   ▼            │
  Inlet B (7) ── feed B ──► N4 ◄─────────────────┘ ── 4-5 ────┘

Junction links 2-3, 2-4, 3-4, 3-5, 4-5 join four nodes: 5 − 4 + 1 = 2
independent loops, (2, 3, 4) and (3, 5, 4).  Each mode carries a fluid of
that kind; oil–gas–water where the mode allows it (black-oil with WOR > 0,
compositional with water and free-water formation).
"""
from __future__ import annotations

from dataclasses import replace

from angelica.core.case import (
    FlowBoundary,
    InletCompositionBC,
    NetworkCase,
    PressureBoundary,
    ThermalBoundary,
)
from angelica.core.components import Pipe
from angelica.properties.black_oil import BlackOilFluid
from angelica.properties.compositional_fluid import CompositionalFluid
from angelica.properties.compressible_fluid import CompressibleFluid
from angelica.properties.eos import PengRobinsonEOS
from angelica.properties.single_component import SingleComponentFluid
from angelica.properties.water_liquid import build_water_thermal_fluid
from angelica.solvers import (
    SteadyBlackOilSolver,
    SteadyCompositionalSolver,
    SteadyCompressibleSolver,
    SteadyIsothermalIncompressibleSolver,
    SteadyNonIsothermalIncompressibleSolver,
)
from angelica.solvers.steady_black_oil import BlackOilSolverSettings
from angelica.solvers.steady_compositional import CompositionalSolverSettings
from angelica.solvers.steady_compressible import CompressibleSolverSettings

INLET_A, INLET_B, OUTLET = 1, 7, 6
LINKS = (   # id, start, end, length scale (× the mode's base length), diameter scale
    ("feed_A", 1, 2, 1.0, 1.0),
    ("feed_B", 7, 4, 0.8, 1.0),
    ("2-3",    2, 3, 1.2, 0.8),
    ("2-4",    2, 4, 0.9, 0.6),
    ("3-4",    3, 4, 0.7, 0.6),
    ("3-5",    3, 5, 1.1, 0.8),
    ("4-5",    4, 5, 1.0, 0.8),
    ("export", 5, 6, 0.5, 1.2),
)
LOOPS = (("2-3", "3-4", "-2-4"), ("3-5", "-4-5", "-3-4"))   # signed link ids around each loop
T_AMBIENT = 15.0

# mode: (base length m, base diameter m, P_A, P_B, P_out (Pa), T_A, T_B (°C), U W/m²K)
MODES = {
    "liquid":        (2_000.0, 0.20, 6.0e5, 5.5e5, 1.0e5, 20.0, 20.0, 0.0),
    "liquid_thermal": (2_000.0, 0.20, 6.0e5, 5.5e5, 1.0e5, 80.0, 60.0, 5.0),
    "gas":           (10_000.0, 0.25, 70e5, 66e5, 40e5, 50.0, 30.0, 3.0),
    "black_oil":     (3_000.0, 0.15, 40e5, 37e5, 10e5, 70.0, 55.0, 3.0),
    "compositional": (10_000.0, 0.25, 80e5, 75e5, 40e5, 90.0, 50.0, 3.0),
}

COMPONENTS = ("methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane", "water")
ZS_A = (0.55, 0.10, 0.10, 0.08, 0.05, 0.05, 0.07)   # rich condensate with water
ZS_B = (0.88, 0.06, 0.03, 0.01, 0.00, 0.00, 0.02)   # lean gas with some water
CRUDE = {"api_gravity": 35.0, "gas_gravity": 0.70, "gor_sc_m3_per_m3": 80.0, "wor_sc_m3_per_m3": 0.4}


def fluid_for(mode: str):
    if mode == "liquid":
        return SingleComponentFluid(density_kg_per_m3=998.0, viscosity_pa_s=1.0e-3)
    if mode == "liquid_thermal":
        return build_water_thermal_fluid()
    if mode == "gas":
        return CompressibleFluid.from_constants(
            eos=PengRobinsonEOS(molecular_weight_kg_per_mol=0.016043, critical_temperature_k=190.56,
                                critical_pressure_pa=4.599e6, acentric_factor=0.011),   # methane
            viscosity_pa_s=1.2e-5, specific_heat_j_per_kg_k=2_200.0,
            thermal_conductivity_w_per_m_k=0.035, reference_pressure_pa=50e5,
            reference_temperature_c=30.0)
    if mode == "black_oil":
        return BlackOilFluid(**CRUDE)
    return CompositionalFluid(components=list(COMPONENTS), default_zs=list(ZS_A))


def pipe(mode: str, cid: str, a: int, b: int, length: float, diameter: float) -> Pipe:
    U = MODES[mode][7]
    return Pipe(component_id=cid, start_node=a, end_node=b, diameter_m=diameter, length_m=length,
                absolute_roughness_m=46e-6, heat_transfer_coefficient_w_per_m2k=U,
                ambient_temperature_c=T_AMBIENT)


def build_case(mode: str) -> NetworkCase:
    L0, D0, pa, pb, pout, ta, tb, _U = MODES[mode]
    comps = tuple(pipe(mode, cid, a, b, L0 * kl, D0 * kd) for cid, a, b, kl, kd in LINKS)
    kw = {}
    if mode == "compositional":
        kw["inlet_composition_bcs"] = (InletCompositionBC(node_id=INLET_A, zs=ZS_A),
                                       InletCompositionBC(node_id=INLET_B, zs=ZS_B))
    if mode != "liquid":
        kw["thermal_inlets"] = (ThermalBoundary(node_id=INLET_A, temperature_c=ta),
                                ThermalBoundary(node_id=INLET_B, temperature_c=tb))
    return NetworkCase(
        name=f"Two-loop network ({mode})",
        fluid_model=fluid_for(mode),
        pressure_inlets=(PressureBoundary(node_id=INLET_A, pressure_pa=pa),
                         PressureBoundary(node_id=INLET_B, pressure_pa=pb)),
        pressure_outlets=(PressureBoundary(node_id=OUTLET, pressure_pa=pout),),
        components=comps,
        **kw,
    )


def solver_for(mode: str, isothermal_c: float | None = None):
    if mode == "liquid":
        return SteadyIsothermalIncompressibleSolver()
    if mode == "liquid_thermal":
        return SteadyNonIsothermalIncompressibleSolver()
    if mode == "gas":
        return SteadyCompressibleSolver(
            compressible_settings=CompressibleSolverSettings(isothermal_temperature_c=isothermal_c))
    if mode == "black_oil":
        return SteadyBlackOilSolver(
            black_oil_settings=BlackOilSolverSettings(isothermal_temperature_c=isothermal_c))
    return SteadyCompositionalSolver(
        compositional_settings=CompositionalSolverSettings(isothermal_temperature_c=isothermal_c))


def single_pipe_case(mode: str, case: NetworkCase, result, cid: str) -> tuple[NetworkCase, int, int, float]:
    """The pipe ``cid`` alone: its upstream node's pressure, temperature and
    composition from the network solution, and the network's flow through it
    imposed at its downstream node.  Returns the case, the upstream and
    downstream node ids, and the network's pressure at the downstream node."""
    comps = {c.component_id: c for c in case.components}
    flows = {c.component_id: f.mass_flow_kg_per_s for c, f in zip(case.components, result.component_flows)}
    p = comps[cid]
    m = flows[cid]
    up, down = (p.start_node, p.end_node) if m >= 0 else (p.end_node, p.start_node)
    alone = replace(p, start_node=up, end_node=down)
    kw = {}
    if mode != "liquid":
        kw["thermal_inlets"] = (ThermalBoundary(node_id=up, temperature_c=result.node_temperatures_c[up]),)
    if mode == "compositional":
        zs = {INLET_A: ZS_A, INLET_B: ZS_B}.get(up) or tuple(result.node_compositions[up])
        kw["inlet_composition_bcs"] = (InletCompositionBC(node_id=up, zs=zs),)
    single = NetworkCase(
        name=f"{cid} alone",
        fluid_model=case.fluid_model,
        pressure_inlets=(PressureBoundary(node_id=up, pressure_pa=result.node_pressures_pa[up]),),
        pressure_outlets=(),
        flow_outlets=(FlowBoundary(node_id=down, mass_flow_kg_per_s=abs(m)),),
        components=(alone,),
        **kw,
    )
    return single, up, down, result.node_pressures_pa[down]
