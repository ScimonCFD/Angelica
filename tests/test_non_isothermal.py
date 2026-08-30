"""Tests for the non-isothermal incompressible solver.

Benchmark: single adiabatic-wall or heat-losing pipe against the analytical
exponential-decay solution for steady convection-diffusion with a uniform
heat-loss source term.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica import (
    FlowBoundary,
    NetworkCase,
    NonIsothermalSolverSettings,
    Pipe,
    PressureBoundary,
    SteadyNonIsothermalIncompressibleSolver,
    ThermalBoundary,
    ThermalFluid,
    build_thermal_dead_oil,
)
from angelica.cases import (
    build_crude_oil_pipeline_thermal_case,
    build_inline_heater_fixed_flow_case,
    build_symmetric_adiabatic_loop_case,
    build_symmetric_heat_loss_loop_case,
    build_thermal_mixing_junction_case,
)
from angelica.closures import HybridScheme, PowerLawScheme, UpwindScheme

# ---------------------------------------------------------------------------
# Shared fixture: single pipe with heat loss
# ---------------------------------------------------------------------------
#
# Fluid: heavy-oil-like, constant properties
#   ρ = 900 kg/m³,  μ = 0.05 Pa·s,  cp = 2000 J/(kg·K),  k = 0.15 W/(m·K)
#
# Pipe: D = 0.01 m, L = 5 m, U = 100 W/(m²K), T_amb = 20 °C
#   Hagen-Poiseuille (Re ≈ 22, well laminar) gives
#     Q ≈ 9.82 × 10⁻⁶ m³/s  for ΔP = 10 000 Pa
#     ṁ ≈ 8.84 × 10⁻³ kg/s
#   NTU = U π D L / (ṁ cp) ≈ 0.89
#   T_out_analytical ≈ 44.7 °C  (inlet 80 °C, ambient 20 °C)
# ---------------------------------------------------------------------------

_FLUID = ThermalFluid.from_constants(
    density_kg_per_m3=900.0,
    viscosity_pa_s=0.05,
    specific_heat_j_per_kg_k=2000.0,
    thermal_conductivity_w_per_m_k=0.15,
)

_PIPE_WITH_HEAT_LOSS = Pipe(
    start_node=1,
    end_node=2,
    diameter_m=0.01,
    length_m=5.0,
    absolute_roughness_m=0.0,
    heat_transfer_coefficient_w_per_m2k=100.0,
    ambient_temperature_c=20.0,
    n_thermal_segments=50,
)

_PIPE_ADIABATIC = Pipe(
    start_node=1,
    end_node=2,
    diameter_m=0.01,
    length_m=5.0,
    absolute_roughness_m=0.0,
    heat_transfer_coefficient_w_per_m2k=0.0,
    ambient_temperature_c=20.0,
    n_thermal_segments=20,
)

_BOUNDARIES = {
    "pressure_inlets": (PressureBoundary(node_id=1, pressure_pa=111_325.0),),
    "pressure_outlets": (PressureBoundary(node_id=2, pressure_pa=101_325.0),),
    "thermal_inlets": (ThermalBoundary(node_id=1, temperature_c=80.0),),
}


def _analytical_T_out(mdot: float) -> float:
    """Continuous-domain analytical outlet temperature (°C)."""
    U, D, L, cp = 100.0, 0.01, 5.0, 2000.0
    T_in, T_amb = 80.0, 20.0
    ntu = U * math.pi * D * L / (mdot * cp)
    return T_amb + (T_in - T_amb) * math.exp(-ntu)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_outlet_temperature_upwind():
    """Upwind scheme: outlet T within 1 °C of analytical exponential decay."""
    case = NetworkCase(
        name="heat_loss_upwind",
        fluid_model=_FLUID,
        components=(_PIPE_WITH_HEAT_LOSS,),
        **_BOUNDARIES,
    )
    solver = SteadyNonIsothermalIncompressibleSolver(
        non_isothermal_settings=NonIsothermalSolverSettings(
            max_temperature_iterations=30,
            temperature_tolerance_k=0.001,
        ),
        convection_scheme=UpwindScheme(),
    )
    result = solver.solve(case)

    assert result.converged

    mdot = result.component_flows[0].mass_flow_kg_per_s
    assert mdot > 0.0  # flow from node 1 → node 2

    T_out = result.node_temperatures_c[2]
    T_expected = _analytical_T_out(mdot)
    assert abs(T_out - T_expected) < 1.0, (
        f"outlet T = {T_out:.2f} °C, analytical = {T_expected:.2f} °C"
    )


@pytest.mark.parametrize("scheme", [HybridScheme(), PowerLawScheme()])
def test_outlet_temperature_alternative_schemes(scheme):
    """Hybrid and power-law schemes give similarly accurate outlet temperature."""
    case = NetworkCase(
        name="heat_loss_scheme",
        fluid_model=_FLUID,
        components=(_PIPE_WITH_HEAT_LOSS,),
        **_BOUNDARIES,
    )
    solver = SteadyNonIsothermalIncompressibleSolver(
        non_isothermal_settings=NonIsothermalSolverSettings(
            temperature_tolerance_k=0.001,
        ),
        convection_scheme=scheme,
    )
    result = solver.solve(case)

    assert result.converged
    mdot = result.component_flows[0].mass_flow_kg_per_s
    T_out = result.node_temperatures_c[2]
    T_expected = _analytical_T_out(mdot)
    assert abs(T_out - T_expected) < 1.5


def test_adiabatic_pipe_temperature_unchanged():
    """With U = 0 W/(m²K), outlet temperature equals inlet temperature."""
    case = NetworkCase(
        name="adiabatic",
        fluid_model=_FLUID,
        components=(_PIPE_ADIABATIC,),
        **_BOUNDARIES,
    )
    solver = SteadyNonIsothermalIncompressibleSolver()
    result = solver.solve(case)

    assert result.converged
    T_out = result.node_temperatures_c[2]
    assert abs(T_out - 80.0) < 0.5, f"adiabatic pipe: T_out = {T_out:.3f} °C, expected 80 °C"


def test_result_has_temperature_field():
    """SolveResult carries node_temperatures_c for all nodes."""
    case = NetworkCase(
        name="temp_field",
        fluid_model=_FLUID,
        components=(_PIPE_WITH_HEAT_LOSS,),
        **_BOUNDARIES,
    )
    result = SteadyNonIsothermalIncompressibleSolver().solve(case)
    assert set(result.node_temperatures_c.keys()) == {1, 2}
    assert result.node_temperatures_c[1] == pytest.approx(80.0, abs=0.1)


def test_no_flow_no_crash():
    """Zero ΔP (equal pressures) should converge without error."""
    case = NetworkCase(
        name="no_flow",
        fluid_model=_FLUID,
        components=(_PIPE_WITH_HEAT_LOSS,),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=101_325.0),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=101_325.0),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=80.0),),
    )
    result = SteadyNonIsothermalIncompressibleSolver().solve(case)
    # no assertion on T — just verify it doesn't raise
    assert result.node_temperatures_c is not None


def test_adiabatic_mixing_junction_matches_exact_mixing_temperature():
    """Two inlet streams should mix exactly at the junction in an adiabatic network."""
    result = SteadyNonIsothermalIncompressibleSolver().solve(
        build_thermal_mixing_junction_case()
    )

    assert result.converged
    expected_mix_c = 40.0
    assert result.node_temperatures_c[3] == pytest.approx(expected_mix_c, abs=0.1)
    assert result.node_temperatures_c[4] == pytest.approx(expected_mix_c, abs=0.1)

    flows = {cf.label: cf.mass_flow_kg_per_s for cf in result.component_flows}
    assert flows["Pipe:hot_branch"] == pytest.approx(1.0, abs=1e-6)
    assert flows["Pipe:cold_branch"] == pytest.approx(2.0, abs=1e-6)
    assert flows["Pipe:mixed_outlet"] == pytest.approx(3.0, abs=1e-6)


def test_inline_heater_fixed_flow_matches_exact_delta_t():
    """With fixed mass flow and adiabatic pipes, ΔT must match Q/(ṁcp)."""
    result = SteadyNonIsothermalIncompressibleSolver().solve(
        build_inline_heater_fixed_flow_case()
    )

    assert result.converged
    expected_delta_t = 50_000.0 / (1.0 * 4182.0)
    actual_delta_t = result.node_temperatures_c[4] - result.node_temperatures_c[1]
    assert actual_delta_t == pytest.approx(expected_delta_t, abs=0.2)

    flows = {cf.label: cf.mass_flow_kg_per_s for cf in result.component_flows}
    assert flows["Pipe:feed_pipe"] == pytest.approx(1.0, abs=1e-4)
    assert flows["HeatSource:heater"] == pytest.approx(1.0, abs=1e-4)
    assert flows["Pipe:exit_pipe"] == pytest.approx(1.0, abs=1e-4)


def test_symmetric_adiabatic_loop_stays_isothermal_and_splits_flow_evenly():
    result = SteadyNonIsothermalIncompressibleSolver().solve(
        build_symmetric_adiabatic_loop_case()
    )

    assert result.converged
    for node_id, temperature_c in result.node_temperatures_c.items():
        assert temperature_c == pytest.approx(60.0, abs=0.1), node_id

    flows = {cf.label: cf.mass_flow_kg_per_s for cf in result.component_flows}
    assert flows["Pipe:upper_branch"] == pytest.approx(1.0, abs=1e-4)
    assert flows["Pipe:lower_branch"] == pytest.approx(1.0, abs=1e-4)
    assert flows["Pipe:upper_return"] == pytest.approx(1.0, abs=1e-4)
    assert flows["Pipe:lower_return"] == pytest.approx(1.0, abs=1e-4)


def test_symmetric_heat_loss_loop_matches_exact_branch_temperature():
    result = SteadyNonIsothermalIncompressibleSolver().solve(
        build_symmetric_heat_loss_loop_case()
    )

    assert result.converged
    mdot_branch = 1.0
    cp = 4182.0
    U = 50.0
    D = 0.025
    L = 500.0
    T_in = 80.0
    T_amb = 20.0
    ntu = U * math.pi * D * L / (mdot_branch * cp)
    expected_branch_out = T_amb + (T_in - T_amb) * math.exp(-ntu)

    assert result.node_temperatures_c[3] == pytest.approx(expected_branch_out, abs=0.3)
    assert result.node_temperatures_c[4] == pytest.approx(expected_branch_out, abs=0.3)
    assert result.node_temperatures_c[5] == pytest.approx(expected_branch_out, abs=0.3)
    assert result.node_temperatures_c[6] == pytest.approx(expected_branch_out, abs=0.3)

    flows = {cf.label: cf.mass_flow_kg_per_s for cf in result.component_flows}
    assert flows["Pipe:upper_branch"] == pytest.approx(1.0, abs=1e-4)
    assert flows["Pipe:lower_branch"] == pytest.approx(1.0, abs=1e-4)


def test_crude_oil_thermal_pipeline_converges_and_cools():
    """Non-isothermal crude oil pipeline: oil cools from 80 °C toward ambient,
    and the solver converges using temperature-dependent Beggs & Robinson properties."""
    result = SteadyNonIsothermalIncompressibleSolver().solve(
        build_crude_oil_pipeline_thermal_case()
    )

    assert result.converged

    # Inlet node is a fixed-temperature thermal boundary
    assert result.node_temperatures_c[1] == pytest.approx(80.0, abs=0.5)

    # All outlet nodes must be cooler than the inlet (heat loss to 15 °C ambient)
    for outlet_nid in (3, 4):
        assert result.node_temperatures_c[outlet_nid] < result.node_temperatures_c[1]

    # Flow must be positive (from inlet to outlets)
    for cf in result.component_flows:
        assert cf.mass_flow_kg_per_s > 0.0


def test_crude_oil_thermal_colder_outlet_than_isothermal_hot():
    """Temperature-dependent viscosity: outlet temperature must lie strictly between
    the inlet temperature (80 °C) and ambient (15 °C), confirming heat loss."""
    result = SteadyNonIsothermalIncompressibleSolver().solve(
        build_crude_oil_pipeline_thermal_case()
    )

    assert result.converged
    T_in = 80.0
    T_amb = 15.0
    for nid in (3, 4):
        T_out = result.node_temperatures_c[nid]
        assert T_amb < T_out < T_in, (
            f"Outlet node {nid}: T={T_out:.2f} °C not between ambient and inlet"
        )


def test_non_convergence_still_returns_consistent_results():
    """When max_temperature_iterations is exhausted without meeting tolerance,
    the returned flow and temperature fields must come from the same final
    synchronous solve — not from mismatched iterations."""
    case = NetworkCase(
        name="forced_nonconvergence",
        fluid_model=_FLUID,
        components=(_PIPE_WITH_HEAT_LOSS,),
        **_BOUNDARIES,
    )
    solver = SteadyNonIsothermalIncompressibleSolver(
        non_isothermal_settings=NonIsothermalSolverSettings(
            max_temperature_iterations=1,   # single outer pass — never enough
            temperature_tolerance_k=1e-15,  # impossibly tight
        ),
    )
    result = solver.solve(case)

    assert not result.converged
    assert len(result.node_temperatures_c) == 2
    assert len(result.node_pressures_pa) == 2
    assert len(result.component_flows) == 1
    # Temperatures must be finite and physically plausible
    for T in result.node_temperatures_c.values():
        assert 10.0 < T < 100.0, f"Temperature {T:.2f} °C out of plausible range"
    # Inlet boundary must be honoured even without convergence
    assert abs(result.node_temperatures_c[1] - 80.0) < 0.5


def test_neumann_zero_gradient_outlet():
    """Zero-gradient Neumann BC at the outlet allows the solver to find the
    natural exit temperature without imposing a Dirichlet constraint.
    For a heat-loss pipe the outlet must lie strictly between ambient and inlet."""
    case = NetworkCase(
        name="neumann_outlet",
        fluid_model=_FLUID,
        components=(_PIPE_WITH_HEAT_LOSS,),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=111_325.0),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=101_325.0),),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=80.0, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=2, temperature_c=0.0, bc_type="zero_gradient"),
        ),
    )
    result = SteadyNonIsothermalIncompressibleSolver(
        non_isothermal_settings=NonIsothermalSolverSettings(temperature_tolerance_k=0.001),
    ).solve(case)

    assert result.converged
    T_out = result.node_temperatures_c[2]
    # Must be cooled by heat loss but not below ambient
    assert 20.0 < T_out < 80.0, f"Neumann outlet T = {T_out:.2f} °C outside plausible range"
    assert T_out < result.node_temperatures_c[1]


def test_invalid_bc_type_raises():
    """ThermalBoundary with an unknown bc_type must raise ValueError immediately."""
    with pytest.raises(ValueError, match="bc_type"):
        ThermalBoundary(node_id=1, temperature_c=80.0, bc_type="typo")


def test_initial_temperature_uses_fixed_bc_not_zero_gradient():
    """When a zero_gradient BC appears before a fixed_temperature BC in
    thermal_inlets, the solver must not seed the temperature field at 0 °C."""
    case = NetworkCase(
        name="bc_order",
        fluid_model=_FLUID,
        components=(_PIPE_WITH_HEAT_LOSS,),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=111_325.0),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=101_325.0),),
        thermal_inlets=(
            ThermalBoundary(node_id=2, temperature_c=0.0, bc_type="zero_gradient"),
            ThermalBoundary(node_id=1, temperature_c=80.0, bc_type="fixed_temperature"),
        ),
    )
    result = SteadyNonIsothermalIncompressibleSolver().solve(case)

    assert result.converged
    assert result.node_temperatures_c[1] == pytest.approx(80.0, abs=0.1)
    assert result.node_temperatures_c[2] > 20.0, (
        "Outlet temperature seeded at 0 °C — initial_temperature bug not fixed"
    )


def test_mesh_convergence_approaches_ntu_analytical():
    """FV discretisation error decreases as n_thermal_segments increases.

    At n=5 the error may be several degrees; at n=100 it must be < 0.5 °C.
    The error sequence must also be monotonically decreasing.
    """
    solver = SteadyNonIsothermalIncompressibleSolver(
        non_isothermal_settings=NonIsothermalSolverSettings(temperature_tolerance_k=0.001),
    )

    errors = []
    for n_segs in (5, 20, 100):
        pipe = Pipe(
            start_node=1,
            end_node=2,
            diameter_m=0.01,
            length_m=5.0,
            absolute_roughness_m=0.0,
            heat_transfer_coefficient_w_per_m2k=100.0,
            ambient_temperature_c=20.0,
            n_thermal_segments=n_segs,
        )
        case = NetworkCase(
            name=f"mesh_conv_n{n_segs}",
            fluid_model=_FLUID,
            components=(pipe,),
            **_BOUNDARIES,
        )
        result = solver.solve(case)
        assert result.converged
        mdot = result.component_flows[0].mass_flow_kg_per_s
        T_out = result.node_temperatures_c[2]
        errors.append(abs(T_out - _analytical_T_out(mdot)))

    # Monotonically decreasing error
    assert errors[0] > errors[1] > errors[2], (
        f"Errors not monotonically decreasing: {errors}"
    )
    # Finest mesh within 0.5 °C
    assert errors[2] < 0.5, f"Finest mesh error = {errors[2]:.3f} °C, expected < 0.5 °C"


def test_series_pipes_chained_ntu():
    """Two pipes in series: intermediate and final temperatures match chained NTU.

    Pipe 1: 1→2, U1=100 W/(m²K), T_amb1=20 °C
    Pipe 2: 2→3, U2=50  W/(m²K), T_amb2=30 °C
    Analytical: apply NTU formula sequentially using the shared mass flow.
    """
    pipe1 = Pipe(
        start_node=1,
        end_node=2,
        diameter_m=0.01,
        length_m=5.0,
        absolute_roughness_m=0.0,
        heat_transfer_coefficient_w_per_m2k=100.0,
        ambient_temperature_c=20.0,
        n_thermal_segments=50,
    )
    pipe2 = Pipe(
        start_node=2,
        end_node=3,
        diameter_m=0.01,
        length_m=5.0,
        absolute_roughness_m=0.0,
        heat_transfer_coefficient_w_per_m2k=50.0,
        ambient_temperature_c=30.0,
        n_thermal_segments=50,
    )
    case = NetworkCase(
        name="series_pipes",
        fluid_model=_FLUID,
        components=(pipe1, pipe2),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=121_325.0),),
        pressure_outlets=(PressureBoundary(node_id=3, pressure_pa=101_325.0),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=80.0),),
    )
    result = SteadyNonIsothermalIncompressibleSolver(
        non_isothermal_settings=NonIsothermalSolverSettings(temperature_tolerance_k=0.001),
    ).solve(case)

    assert result.converged
    mdot = result.component_flows[0].mass_flow_kg_per_s

    cp = 2000.0
    D, L = 0.01, 5.0
    T_in = 80.0
    T_amb1, U1 = 20.0, 100.0
    T_amb2, U2 = 30.0, 50.0

    ntu1 = U1 * math.pi * D * L / (mdot * cp)
    T_mid_exact = T_amb1 + (T_in - T_amb1) * math.exp(-ntu1)

    ntu2 = U2 * math.pi * D * L / (mdot * cp)
    T_out_exact = T_amb2 + (T_mid_exact - T_amb2) * math.exp(-ntu2)

    assert abs(result.node_temperatures_c[2] - T_mid_exact) < 1.0, (
        f"Intermediate T = {result.node_temperatures_c[2]:.2f} °C, "
        f"analytical = {T_mid_exact:.2f} °C"
    )
    assert abs(result.node_temperatures_c[3] - T_out_exact) < 1.0, (
        f"Outlet T = {result.node_temperatures_c[3]:.2f} °C, "
        f"analytical = {T_out_exact:.2f} °C"
    )


def test_thermal_coupling_warm_inlet_increases_flow():
    """Temperature-dependent viscosity: hot inlet produces more flow than cold inlet.

    Dead oil (API 30) is highly viscous at low temperature. Feeding the same
    pressure difference with T_in=80°C vs T_in=20°C should produce significantly
    more flow at 80°C because the lower viscosity reduces frictional resistance.
    """
    fluid = build_thermal_dead_oil(api_gravity=30)
    pipe = Pipe(
        start_node=1,
        end_node=2,
        diameter_m=0.05,
        length_m=100.0,
        absolute_roughness_m=0.0,
        heat_transfer_coefficient_w_per_m2k=0.0,
        ambient_temperature_c=20.0,
        n_thermal_segments=10,
    )

    def _run(t_inlet_c: float) -> float:
        case = NetworkCase(
            name=f"visc_coupling_{t_inlet_c:.0f}c",
            fluid_model=fluid,
            components=(pipe,),
            pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=200_000.0),),
            pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=101_325.0),),
            thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=t_inlet_c),),
        )
        result = SteadyNonIsothermalIncompressibleSolver().solve(case)
        assert result.converged
        return result.component_flows[0].mass_flow_kg_per_s

    mdot_hot = _run(80.0)
    mdot_cold = _run(20.0)

    assert mdot_hot > mdot_cold * 2.0, (
        f"Hot flow {mdot_hot:.4f} kg/s should be at least 2× cold flow {mdot_cold:.4f} kg/s"
    )


def test_global_energy_conservation_adiabatic_network():
    """In an adiabatic network, inlet enthalpy flux must equal outlet enthalpy flux.

    Uses the thermal mixing junction case (fixed-flow BCs, adiabatic pipes).
    Enthalpy flux H = ṁ · cp · T.  Conservation: ΣH_in = ΣH_out.
    """
    from angelica.cases import build_thermal_mixing_junction_case
    result = SteadyNonIsothermalIncompressibleSolver().solve(
        build_thermal_mixing_junction_case()
    )

    assert result.converged

    cp = 4182.0  # matches fluid in build_thermal_mixing_junction_case
    flows = {cf.label: cf.mass_flow_kg_per_s for cf in result.component_flows}

    mdot_hot = flows["Pipe:hot_branch"]   # 1→3
    mdot_cold = flows["Pipe:cold_branch"]  # 2→3
    mdot_out = flows["Pipe:mixed_outlet"]  # 3→4

    T_hot_in = result.node_temperatures_c[1]
    T_cold_in = result.node_temperatures_c[2]
    T_out = result.node_temperatures_c[4]

    H_in = cp * (mdot_hot * T_hot_in + mdot_cold * T_cold_in)
    H_out = cp * mdot_out * T_out

    assert abs(H_in - H_out) / max(abs(H_in), 1.0) < 0.01, (
        f"Energy not conserved: H_in={H_in:.1f} W, H_out={H_out:.1f} W"
    )


def test_reverse_flow_temperature_transported_correctly():
    """When flow is reversed (end_node → start_node), temperature must follow flow.

    Pipe defined as 1→2 but pressure drives flow from node 2 to node 1.
    Thermal BC is imposed at node 2 (the true inlet).  Node 1 (the true outlet)
    must reach the inlet temperature in an adiabatic pipe.
    """
    pipe = Pipe(
        start_node=1,
        end_node=2,
        diameter_m=0.01,
        length_m=5.0,
        absolute_roughness_m=0.0,
        heat_transfer_coefficient_w_per_m2k=0.0,
        ambient_temperature_c=20.0,
        n_thermal_segments=20,
    )
    case = NetworkCase(
        name="reverse_flow",
        fluid_model=_FLUID,
        components=(pipe,),
        pressure_inlets=(PressureBoundary(node_id=2, pressure_pa=111_325.0),),
        pressure_outlets=(PressureBoundary(node_id=1, pressure_pa=101_325.0),),
        thermal_inlets=(ThermalBoundary(node_id=2, temperature_c=80.0),),
    )
    result = SteadyNonIsothermalIncompressibleSolver().solve(case)

    assert result.converged

    # Flow must be negative (component direction 1→2 but actual flow 2→1)
    mdot = result.component_flows[0].mass_flow_kg_per_s
    assert mdot < 0.0, f"Expected reversed flow, got mdot = {mdot:.6f} kg/s"

    # Outlet node (node 1) should carry the inlet temperature from node 2
    T_outlet = result.node_temperatures_c[1]
    assert abs(T_outlet - 80.0) < 0.5, (
        f"Reversed adiabatic flow: T_outlet = {T_outlet:.2f} °C, expected ~80 °C"
    )


def test_cengel_example_8_3_oil_pipeline_through_lake():
    """Textbook benchmark: Cengel & Ghajar, Heat and Mass Transfer, Example 8-3.

    Oil at 20 °C flows at 2 m/s through a 30-cm-diameter, 200-m-long pipeline
    submerged in an icy lake (T_surface = 0 °C).  The book computes the
    convective coefficient from the Gnielinski-entry-length Nu correlation
    (Nu = 33.7 → h = 16.3 W/m²·K) and the exit temperature via the
    NTU/effectiveness formula:

        T_e = T_s - (T_s - T_i) exp(-h A_s / (ṁ cp)) = 19.74 °C

    We drive Angelica with the same h as a fixed overall U coefficient and a
    fixed-flow boundary condition (ṁ = 125.6 kg/s), then check that the
    solver reproduces the published exit temperature within 0.05 °C.

    Reference:
        Cengel, Y.A. & Ghajar, A.J., *Heat and Mass Transfer*, 5th ed.,
        McGraw-Hill, 2015, Example 8-3, pp. 486-488.
    """
    fluid = ThermalFluid.from_constants(
        density_kg_per_m3=888.1,
        viscosity_pa_s=888.1 * 9.429e-4,   # ρν, from Table A-13 at 20 °C
        specific_heat_j_per_kg_k=1880.0,
        thermal_conductivity_w_per_m_k=0.145,
    )
    pipe = Pipe(
        start_node=1,
        end_node=2,
        diameter_m=0.3,
        length_m=200.0,
        absolute_roughness_m=0.0,
        heat_transfer_coefficient_w_per_m2k=16.3,   # h from Nu = 33.7 (book)
        ambient_temperature_c=0.0,                   # icy lake
        n_thermal_segments=100,
    )
    case = NetworkCase(
        name="cengel_8_3_oil_lake",
        fluid_model=fluid,
        components=(pipe,),
        pressure_inlets=(),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=101_325.0),),
        flow_inlets=(FlowBoundary(node_id=1, mass_flow_kg_per_s=125.6),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=20.0),),
    )
    result = SteadyNonIsothermalIncompressibleSolver(
        non_isothermal_settings=NonIsothermalSolverSettings(temperature_tolerance_k=1e-4),
    ).solve(case)

    assert result.converged
    T_out = result.node_temperatures_c[2]
    assert abs(T_out - 19.74) < 0.05, (
        f"Cengel Ex. 8-3: T_out = {T_out:.4f} °C, expected 19.74 °C"
    )
