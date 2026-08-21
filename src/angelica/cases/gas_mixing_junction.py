"""Two-inlet gas mixing junction case for the compositional solver.

Topology
--------

    Node 1 (P=100 bar, rich gas)  ──── Pipe A (500 m) ────┐
                                                            ├── Node 3 ─── Pipe C (1000 m) ── Node 4 (P=10 bar)
    Node 2 (P=100 bar, lean gas)  ──── Pipe B (500 m) ────┘

The two inlets carry different mole-fraction compositions; the solver
propagates them through the network and mixes them at the junction node.
"""

from __future__ import annotations

from angelica.core.case import (
    InletCompositionBC,
    NetworkCase,
    PressureBoundary,
    ThermalBoundary,
)
from angelica.core.components import Pipe
from angelica.properties.compositional_fluid import CompositionalFluid


def build_gas_mixing_junction_case(
    p_inlet_pa: float = 100e5,
    p_outlet_pa: float = 10e5,
    inlet_temperature_c: float = 60.0,
    pipe_diameter_m: float = 0.1,
) -> NetworkCase:
    """Two-inlet gas-mixing junction case for the compositional solver.

    Two high-pressure inlet nodes (nodes 1 and 2) feed a central junction
    (node 3) through 500 m pipes.  A 1000 m delivery pipe runs from the
    junction to a single low-pressure outlet (node 4).

    Inlet 1 carries a methane-rich gas (CH₄ 90%, C₂H₆ 10%).
    Inlet 2 carries a leaner gas (CH₄ 60%, C₂H₆ 40%).

    The solver will mix these compositions at node 3 and propagate the
    blended stream to the outlet.

    Args:
        p_inlet_pa: Inlet pressure at nodes 1 and 2 (Pa). Default 100 bar.
        p_outlet_pa: Outlet pressure at node 4 (Pa). Default 10 bar.
        inlet_temperature_c: Fixed temperature at both inlets (°C).
        pipe_diameter_m: Internal diameter of all three pipes (m).

    Returns:
        A :class:`~angelica.core.case.NetworkCase` ready for
        :class:`~angelica.solvers.SteadyCompositionalSolver`.
    """
    components_list = ["methane", "ethane"]

    fluid = CompositionalFluid(
        components=components_list,
        default_zs=[0.75, 0.25],  # fallback composition (blend of the two inlets)
    )

    roughness = 4.6e-5  # commercial steel, m

    pipe_a = Pipe(
        start_node=1, end_node=3,
        length_m=500.0, diameter_m=pipe_diameter_m,
        absolute_roughness_m=roughness,
    )
    pipe_b = Pipe(
        start_node=2, end_node=3,
        length_m=500.0, diameter_m=pipe_diameter_m,
        absolute_roughness_m=roughness,
    )
    pipe_c = Pipe(
        start_node=3, end_node=4,
        length_m=1000.0, diameter_m=pipe_diameter_m,
        absolute_roughness_m=roughness,
    )

    return NetworkCase(
        name="gas_mixing_junction",
        fluid_model=fluid,
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=p_inlet_pa),
            PressureBoundary(node_id=2, pressure_pa=p_inlet_pa),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=4, pressure_pa=p_outlet_pa),
        ),
        components=(pipe_a, pipe_b, pipe_c),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=inlet_temperature_c, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=2, temperature_c=inlet_temperature_c, bc_type="fixed_temperature"),
        ),
        inlet_composition_bcs=(
            InletCompositionBC(node_id=1, zs=(0.90, 0.10)),  # rich gas
            InletCompositionBC(node_id=2, zs=(0.60, 0.40)),  # leaner gas
        ),
    )
