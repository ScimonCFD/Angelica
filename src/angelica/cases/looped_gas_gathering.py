"""Two-source looped gas gathering network for the compositional solver.

Topology
--------

    Node 1 ──PipeA (5 km)──→ Node 3 ──PipeC (8 km)──→ Node 5
    (inlet A, 120 bar)         │     rich gas              │
                              PipeD (3 km, loop)          PipeE (8 km)
                               │                            │
    Node 2 ──PipeB (5 km)──→ Node 4 ─────────────────────┘
    (inlet B, 110 bar)
    lean gas                                          Node 5 = outlet (20 bar)

PipeD is the loop pipe.  Its flow direction — and thus how the two
compositions blend at Node 4 — emerges from the hydraulic solution.
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


def build_looped_gas_gathering_case(
    p_inlet_a_pa: float = 120e5,
    p_inlet_b_pa: float = 110e5,
    p_outlet_pa:  float =  20e5,
    t_inlet_a_c:  float =  70.0,
    t_inlet_b_c:  float =  60.0,
) -> NetworkCase:
    """Two-source looped gas gathering network (three-component EOS).

    Two sources with different compositions feed a looped mid-section.
    The loop pipe (PipeD, connecting the two junction nodes) carries gas
    from the higher-pressure junction to the lower-pressure one; its
    direction is determined by the hydraulic solution.

    Components:  methane, ethane, propane

    Source A (rich gas)  :  CH₄ 90 %,  C₂H₆  8 %,  C₃H₈  2 %
    Source B (leaner gas):  CH₄ 60 %,  C₂H₆ 30 %,  C₃H₈ 10 %

    Returns a :class:`~angelica.core.case.NetworkCase` ready for
    :class:`~angelica.solvers.SteadyCompositionalSolver`.
    """
    components = ["methane", "ethane", "propane"]

    fluid = CompositionalFluid(
        components  = components,
        default_zs  = [0.75, 0.19, 0.06],
    )

    roughness = 46e-6   # commercial steel
    U         = 5.0     # W/(m²·K)
    T_amb     = 15.0    # °C

    def pipe(cid, s, e, length_m, d_m):
        return Pipe(
            component_id                        = cid,
            start_node                          = s,
            end_node                            = e,
            length_m                            = length_m,
            diameter_m                          = d_m,
            absolute_roughness_m                = roughness,
            heat_transfer_coefficient_w_per_m2k = U,
            ambient_temperature_c               = T_amb,
        )

    return NetworkCase(
        name             = "Looped Gas Gathering Network",
        fluid_model      = fluid,
        pressure_inlets  = (
            PressureBoundary(node_id=1, pressure_pa=p_inlet_a_pa),
            PressureBoundary(node_id=2, pressure_pa=p_inlet_b_pa),
        ),
        pressure_outlets = (PressureBoundary(node_id=5, pressure_pa=p_outlet_pa),),
        components       = (
            pipe("pipe_A", 1, 3, 5_000.0, 0.12),
            pipe("pipe_B", 2, 4, 5_000.0, 0.12),
            pipe("pipe_C", 3, 5, 8_000.0, 0.15),
            pipe("pipe_D", 3, 4, 3_000.0, 0.10),   # loop pipe
            pipe("pipe_E", 4, 5, 8_000.0, 0.15),
        ),
        thermal_inlets   = (
            ThermalBoundary(node_id=1, temperature_c=t_inlet_a_c, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=2, temperature_c=t_inlet_b_c, bc_type="fixed_temperature"),
        ),
        inlet_composition_bcs = (
            InletCompositionBC(node_id=1, zs=(0.90, 0.08, 0.02)),   # rich gas
            InletCompositionBC(node_id=2, zs=(0.60, 0.30, 0.10)),   # leaner gas
        ),
    )
