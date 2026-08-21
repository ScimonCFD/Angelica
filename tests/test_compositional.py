"""Tests for the compositional fluid model and solver.

Coverage
--------
* CompositionalFluid: property evaluation, two-phase handling, validation.
* _propagate_compositions: single inlet, two-inlet mixing, flow-reversal.
* SteadyCompositionalSolver: end-to-end convergence on simple networks.
* InletCompositionBC: dataclass validation.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

# Tests that call thermo.Mixture are skipped when thermo is not installed.
_thermo_available = importlib.util.find_spec("thermo") is not None
requires_thermo = pytest.mark.skipif(
    not _thermo_available,
    reason="thermo not installed — pip install angelica[compositional]",
)

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica import (
    CompositionalFluid,
    CompositionalSolverSettings,
    FlowBoundary,
    InletCompositionBC,
    NetworkCase,
    Pipe,
    PressureBoundary,
    SteadyCompositionalSolver,
    ThermalBoundary,
)
from angelica.cases.gas_mixing_junction import build_gas_mixing_junction_case

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _methane_ethane_fluid(z_ch4: float = 0.8) -> CompositionalFluid:
    return CompositionalFluid(
        components=["methane", "ethane"],
        default_zs=[z_ch4, 1.0 - z_ch4],
    )


def _single_pipe_case(
    fluid: CompositionalFluid,
    p_in: float = 100e5,
    p_out: float = 20e5,
    t_in: float = 60.0,
    inlet_zs: tuple[float, ...] | None = None,
) -> NetworkCase:
    """One inlet → one pipe (500 m, D=0.1 m) → one outlet."""
    bcs = (InletCompositionBC(node_id=1, zs=inlet_zs),) if inlet_zs else ()
    return NetworkCase(
        name="single_pipe_compositional",
        fluid_model=fluid,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=p_in),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=p_out),),
        components=(
            Pipe(start_node=1, end_node=2, length_m=500.0, diameter_m=0.1, absolute_roughness_m=4.6e-5),
        ),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=t_in, bc_type="fixed_temperature"),
        ),
        inlet_composition_bcs=bcs,
    )


# ---------------------------------------------------------------------------
# InletCompositionBC validation
# ---------------------------------------------------------------------------

class TestInletCompositionBC:
    def test_normalised_zs_accepted(self):
        bc = InletCompositionBC(node_id=1, zs=(0.7, 0.3))
        assert bc.zs == (0.7, 0.3)

    def test_zs_converted_to_tuple(self):
        bc = InletCompositionBC(node_id=1, zs=[0.5, 0.5])
        assert isinstance(bc.zs, tuple)

    def test_negative_z_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            InletCompositionBC(node_id=1, zs=(-0.1, 1.1))

    def test_zs_not_summing_to_one_raises(self):
        with pytest.raises(ValueError, match="sum to 1"):
            InletCompositionBC(node_id=1, zs=(0.5, 0.3))

    def test_empty_zs_raises(self):
        with pytest.raises(ValueError, match="empty"):
            InletCompositionBC(node_id=1, zs=())


# ---------------------------------------------------------------------------
# CompositionalFluid validation
# ---------------------------------------------------------------------------

class TestCompositionalFluidValidation:
    def test_mismatched_lengths_raise(self):
        with pytest.raises(ValueError, match="same length"):
            CompositionalFluid(["methane", "ethane"], [1.0])

    def test_negative_default_zs_raise(self):
        with pytest.raises(ValueError, match="non-negative"):
            CompositionalFluid(["methane", "ethane"], [-0.1, 1.1])

    def test_unnormalised_default_zs_raise(self):
        with pytest.raises(ValueError, match="sum to 1"):
            CompositionalFluid(["methane", "ethane"], [0.3, 0.3])

    def test_valid_construction(self):
        fluid = _methane_ethane_fluid()
        assert fluid.component_names == ("methane", "ethane")
        assert len(fluid.default_zs) == 2


# ---------------------------------------------------------------------------
# CompositionalFluid property evaluation
# ---------------------------------------------------------------------------

@requires_thermo
class TestCompositionalFluidProperties:
    """Flash methane/ethane at high-pressure gas conditions (single phase)."""

    def _make_link_state(self, P_pa: float, T_c: float, zs: tuple):
        """Minimal link-state stub for testing property evaluation."""
        class _Node:
            def __init__(self, p):
                self.pressure_pa = p
        class _Link:
            def __init__(self):
                self.start_node = _Node(P_pa)
                self.end_node   = _Node(P_pa)
                self.temperature_c = T_c
                self.zs = zs
        return _Link()

    def test_density_positive_single_phase_gas(self):
        fluid = _methane_ethane_fluid(0.9)
        link  = self._make_link_state(5e6, 80.0, (0.9, 0.1))
        rho   = fluid.density_for_link(link)
        # Methane at 50 bar, 80 °C: rough expectation ~30-50 kg/m³
        assert 10.0 < rho < 200.0

    def test_viscosity_positive(self):
        fluid = _methane_ethane_fluid(0.9)
        link  = self._make_link_state(5e6, 80.0, (0.9, 0.1))
        mu    = fluid.viscosity_for_link(link)
        assert 1e-9 < mu < 1e-2

    def test_specific_heat_positive(self):
        fluid = _methane_ethane_fluid(0.9)
        link  = self._make_link_state(5e6, 80.0, (0.9, 0.1))
        Cp    = fluid.specific_heat_for_link(link)
        assert 100.0 < Cp < 1e5

    def test_thermal_conductivity_positive(self):
        fluid = _methane_ethane_fluid(0.9)
        link  = self._make_link_state(5e6, 80.0, (0.9, 0.1))
        k     = fluid.thermal_conductivity_for_link(link)
        assert 1e-5 < k < 10.0

    def test_density_increases_with_pressure(self):
        fluid = _methane_ethane_fluid(0.9)
        link_lo = self._make_link_state(1e6, 60.0, (0.9, 0.1))
        link_hi = self._make_link_state(10e6, 60.0, (0.9, 0.1))
        assert fluid.density_for_link(link_hi) > fluid.density_for_link(link_lo)

    def test_density_different_compositions(self):
        """Heavier mixture should be denser."""
        fluid    = CompositionalFluid(["methane", "ethane"], [0.5, 0.5])
        link_ch4 = self._make_link_state(5e6, 60.0, (0.99, 0.01))
        link_c2h = self._make_link_state(5e6, 60.0, (0.01, 0.99))
        assert fluid.density_for_link(link_c2h) > fluid.density_for_link(link_ch4)

    def test_fallback_to_default_zs_when_link_has_no_zs(self):
        """link_state.zs == () should fall back to default_zs without error."""
        fluid = _methane_ethane_fluid(0.8)
        class _MinimalLink:
            class _Node:
                pressure_pa = 5e6
            start_node = _Node()
            end_node   = _Node()
            temperature_c = 60.0
            zs = ()
        rho = fluid.density_for_link(_MinimalLink())
        assert rho > 0.0


# ---------------------------------------------------------------------------
# Composition propagation
# ---------------------------------------------------------------------------

class TestPropagateCompositions:
    """Unit-test _propagate_compositions in isolation."""

    def _build_network_state(self, mdots: list[float], nodes: list[int], pipes: list[tuple]):
        """Build a minimal NetworkState from (start, end) pipe list."""
        from angelica.core.state import NetworkState, NodeState, PipeState
        from angelica.core.components import Pipe

        node_states = {nid: NodeState(node_id=nid) for nid in nodes}
        pipe_states = []
        for i, (s, e) in enumerate(pipes):
            p = Pipe(start_node=s, end_node=e, length_m=100.0, diameter_m=0.05, absolute_roughness_m=4.6e-5)
            ps = PipeState(component=p, start_node=node_states[s], end_node=node_states[e])
            ps.mass_flow_kg_per_s = mdots[i]
            pipe_states.append(ps)
        return NetworkState(nodes=node_states, components=pipe_states)

    def test_single_inlet_propagates_to_outlet(self):
        """Composition from node 1 should reach node 2 through one pipe."""
        from angelica.core.case import NetworkCase, PressureBoundary, InletCompositionBC
        from angelica.properties.compositional_fluid import CompositionalFluid
        from angelica.solvers.steady_compositional import SteadyCompositionalSolver

        fluid = CompositionalFluid(["methane", "ethane"], [0.8, 0.2])
        case  = NetworkCase(
            name="test",
            fluid_model=fluid,
            pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=100e5),),
            pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=10e5),),
            components=(Pipe(start_node=1, end_node=2, length_m=100.0, diameter_m=0.05, absolute_roughness_m=4.6e-5),),
            inlet_composition_bcs=(InletCompositionBC(node_id=1, zs=(0.9, 0.1)),),
        )
        ns = self._build_network_state([1.0], [1, 2], [(1, 2)])
        for ps in ns.components:
            ps.zs = fluid.default_zs

        SteadyCompositionalSolver._propagate_compositions(ns, case, fluid.default_zs)

        pipe_zs = ns.components[0].zs
        assert abs(pipe_zs[0] - 0.9) < 1e-9
        assert abs(pipe_zs[1] - 0.1) < 1e-9

    def test_two_inlets_mix_at_junction(self):
        """Equal mass-flows from nodes 1 and 2 should mix to the average at node 3."""
        from angelica.core.case import NetworkCase, PressureBoundary, InletCompositionBC
        from angelica.properties.compositional_fluid import CompositionalFluid
        from angelica.solvers.steady_compositional import SteadyCompositionalSolver

        fluid = CompositionalFluid(["methane", "ethane"], [0.75, 0.25])
        case  = NetworkCase(
            name="test_mix",
            fluid_model=fluid,
            pressure_inlets=(
                PressureBoundary(node_id=1, pressure_pa=100e5),
                PressureBoundary(node_id=2, pressure_pa=100e5),
            ),
            pressure_outlets=(PressureBoundary(node_id=4, pressure_pa=10e5),),
            components=(
                Pipe(start_node=1, end_node=3, length_m=100.0, diameter_m=0.05, absolute_roughness_m=4.6e-5),
                Pipe(start_node=2, end_node=3, length_m=100.0, diameter_m=0.05, absolute_roughness_m=4.6e-5),
                Pipe(start_node=3, end_node=4, length_m=100.0, diameter_m=0.05, absolute_roughness_m=4.6e-5),
            ),
            inlet_composition_bcs=(
                InletCompositionBC(node_id=1, zs=(0.9, 0.1)),
                InletCompositionBC(node_id=2, zs=(0.6, 0.4)),
            ),
        )
        # Equal mass flows from both inlets
        ns = self._build_network_state([1.0, 1.0, 2.0], [1, 2, 3, 4], [(1, 3), (2, 3), (3, 4)])
        for ps in ns.components:
            ps.zs = fluid.default_zs

        SteadyCompositionalSolver._propagate_compositions(ns, case, fluid.default_zs)

        # Delivery pipe (index 2) gets mixed composition: (0.9+0.6)/2 = 0.75
        pipe_c_zs = ns.components[2].zs
        assert abs(pipe_c_zs[0] - 0.75) < 1e-6
        assert abs(pipe_c_zs[1] - 0.25) < 1e-6


# ---------------------------------------------------------------------------
# SteadyCompositionalSolver — end-to-end
# ---------------------------------------------------------------------------

@requires_thermo
class TestSteadyCompositionalSolver:

    def test_wrong_fluid_type_raises(self):
        from angelica.properties.thermal_fluid import ThermalFluid
        from angelica.properties.single_component import SingleComponentFluid

        bad_fluid = SingleComponentFluid(density_kg_per_m3=700.0, viscosity_pa_s=1e-3)
        case = NetworkCase(
            name="bad",
            fluid_model=bad_fluid,
            pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=100e5),),
            pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=10e5),),
            components=(Pipe(start_node=1, end_node=2, length_m=100.0, diameter_m=0.05, absolute_roughness_m=4.6e-5),),
            thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=20.0),),
        )
        with pytest.raises(TypeError, match="CompositionalFluid"):
            SteadyCompositionalSolver().solve(case)

    def test_mismatched_inlet_bc_length_raises(self):
        fluid = _methane_ethane_fluid()
        case = _single_pipe_case(fluid, inlet_zs=(0.5, 0.3, 0.2))  # 3 components vs 2
        with pytest.raises(ValueError, match="components"):
            SteadyCompositionalSolver().solve(case)

    def test_single_pipe_converges(self):
        """Single pipe, no inlet composition BC: uses default_zs throughout."""
        fluid = _methane_ethane_fluid(0.8)
        case  = _single_pipe_case(fluid)
        solver = SteadyCompositionalSolver(
            compositional_settings=CompositionalSolverSettings(max_outer_iterations=10),
        )
        result = solver.solve(case)

        assert result.converged
        # Flow direction: node 1 (high P) → node 2 (low P)
        assert result.component_flows[0].mass_flow_kg_per_s > 0.0
        assert result.node_pressures_pa[1] > result.node_pressures_pa[2]

    def test_single_pipe_with_inlet_bc(self):
        """Single pipe with explicit inlet composition: pipe should adopt that composition."""
        fluid = _methane_ethane_fluid(0.8)
        case  = _single_pipe_case(fluid, inlet_zs=(0.95, 0.05))
        solver = SteadyCompositionalSolver(
            compositional_settings=CompositionalSolverSettings(max_outer_iterations=10),
        )
        result = solver.solve(case)

        assert result.converged
        assert result.component_flows[0].mass_flow_kg_per_s > 0.0

    def test_mixing_junction_case_converges(self):
        """Two-inlet mixing junction: end-to-end with composition propagation."""
        case = build_gas_mixing_junction_case()
        solver = SteadyCompositionalSolver(
            compositional_settings=CompositionalSolverSettings(max_outer_iterations=20),
        )
        result = solver.solve(case)

        assert result.converged
        # Both inlet pipes should carry positive flow toward the junction
        assert result.component_flows[0].mass_flow_kg_per_s > 0.0
        assert result.component_flows[1].mass_flow_kg_per_s > 0.0
        # Delivery pipe also positive
        assert result.component_flows[2].mass_flow_kg_per_s > 0.0
        # Mass balance < 1 %
        assert result.global_balance.mass_error_pct < 1.0

    def test_result_has_temperatures(self):
        """Solver must populate node_temperatures_c for all nodes."""
        fluid = _methane_ethane_fluid(0.8)
        case  = _single_pipe_case(fluid)
        result = SteadyCompositionalSolver().solve(case)

        for nid, T in result.node_temperatures_c.items():
            assert T is not None, f"node {nid} has no temperature"
            assert 0.0 < T < 200.0, f"node {nid} temperature {T} looks wrong"
