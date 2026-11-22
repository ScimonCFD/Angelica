"""Tests for the compositional fluid model and solver.

Coverage
--------
* CompositionalFluid: property evaluation, two-phase handling, validation.
* _propagate_compositions: single inlet, two-inlet mixing, flow-reversal.
* SteadyCompositionalSolver: end-to-end convergence on simple networks.
* InletCompositionBC: dataclass validation.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica import (
    CompositionalFluid,
    CompositionalSolverSettings,
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
# PressureBoundary validation
# ---------------------------------------------------------------------------

class TestPressureBoundary:
    def test_non_negative_pressure_accepted(self):
        pb = PressureBoundary(node_id=1, pressure_pa=101325.0)
        assert pb.pressure_pa == 101325.0

    def test_zero_pressure_accepted(self):
        pb = PressureBoundary(node_id=1, pressure_pa=0.0)
        assert pb.pressure_pa == 0.0

    def test_negative_pressure_raises(self):
        with pytest.raises(ValueError, match=">= 0"):
            PressureBoundary(node_id=1, pressure_pa=-1.0)


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
        from angelica.core.components import Pipe
        from angelica.core.state import NetworkState, NodeState, PipeState

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
        from angelica.core.case import InletCompositionBC, NetworkCase, PressureBoundary
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
        from angelica.core.case import InletCompositionBC, NetworkCase, PressureBoundary
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

class TestSteadyCompositionalSolver:

    def test_wrong_fluid_type_raises(self):
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


class TestFlashCacheResolution:
    """The flash memo rounds (P, T, z).  One rounding step must move a two-phase
    density far less than the outer-loop density tolerance (1e-4); with the
    old 0.1 K step it moved ~5e-4 and the GUI version of tutorial 06 locked
    into a two-state limit cycle and never converged."""

    def test_rounding_step_below_density_tolerance(self):
        from angelica.properties import compositional_fluid as cfm

        names = ("methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane", "water")
        zs = (0.658, 0.088, 0.078, 0.058, 0.035, 0.035, 0.048)
        P, T = 42.5e5, 26.15
        step_T = 10.0 ** (-cfm._T_DIGITS)
        rho0 = cfm._flash_properties(names, *cfm._cache_key_pt(P, T), zs, "PR")[0]
        rho1 = cfm._flash_properties(names, *cfm._cache_key_pt(P, T + step_T), zs, "PR")[0]
        rho2 = cfm._flash_properties(names, *cfm._cache_key_pt(P + cfm._P_STEP_PA, T), zs, "PR")[0]
        assert abs(rho1 - rho0) / rho0 < 1e-5
        assert abs(rho2 - rho0) / rho0 < 1e-5


class TestGasViscosityPressureDependence:
    """Gas viscosity must include the dense-gas contribution.

    Reference values: CoolProp 6.6 reference correlations (methane: Quiñones-
    Cisneros & Huber; ethane: Friend et al.).  Before v1.10.5 the compositional
    model used thermo's low-pressure mixing rule only (methane at 300 K stayed at
    11.1 μPa·s up to 150 bar, −31 %), and a dense supercritical phase labelled
    "liquid" got a 1e-5 Pa·s placeholder (ethane at 350 K, 100 bar: −59 %).
    """

    @staticmethod
    def _mu(names, zs, P_pa, T_K):
        from angelica.properties import compositional_fluid as cfm
        return cfm._flash_properties(tuple(names), P_pa, T_K - 273.15, tuple(zs), "PR")[1]

    @pytest.mark.parametrize("P_bar, ref_upas, tol", [
        (1, 11.255, 0.02), (50, 12.173, 0.02), (100, 13.815, 0.02), (150, 16.163, 0.04),
    ])
    def test_methane_300K_vs_reference(self, P_bar, ref_upas, tol):
        mu = self._mu(["methane"], [1.0], P_bar * 1e5, 300.0) * 1e6
        assert abs(mu - ref_upas) / ref_upas < tol, f"{mu:.3f} vs {ref_upas}"

    def test_methane_ethane_mixture_100bar(self):
        mu = self._mu(["methane", "ethane"], [0.9, 0.1], 100e5, 300.0) * 1e6
        assert abs(mu - 14.225) / 14.225 < 0.02

    def test_dense_supercritical_ethane_not_placeholder(self):
        mu = self._mu(["ethane"], [1.0], 100e5, 350.0) * 1e6
        assert abs(mu - 24.42) / 24.42 < 0.05, mu

    def test_viscosity_increases_with_pressure(self):
        mus = [self._mu(["methane", "ethane"], [0.8, 0.2], P, 333.15) for P in (10e5, 40e5, 70e5, 100e5)]
        assert all(b > a for a, b in zip(mus, mus[1:]))

    def test_liquid_viscosity_unchanged(self):
        """Pure HC liquids keep thermo's (DIPPR-based) value: n-hexane, 25 °C."""
        mu = self._mu(["n-hexane"], [1.0], 5e5, 298.15) * 1e3
        assert abs(mu - 0.2993) / 0.2993 < 0.02

    def test_two_phase_mcadams_uses_gas_mass_fraction(self):
        """thermo's own blend took phase 0 (here the liquid) as the gas quality."""
        from angelica.properties import compositional_fluid as cfm
        names = ("methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane")
        zs = [v / 0.93 for v in (0.55, 0.10, 0.10, 0.08, 0.05, 0.05)]
        fo = cfm._get_flash_obj(names, "PR")
        res = fo.flash(T=303.15, P=50e5, zs=zs)
        assert res.phase_count == 2
        i_gas = min(range(2), key=lambda i: res.phases[i].rho_mass())
        gas, liq = res.phases[i_gas], res.phases[1 - i_gas]
        x_g = res.betas_mass[i_gas]
        mu_g = cfm._phase_viscosity(gas, 303.15, 50e5, fo.constants, fo.correlations, True)
        mu_l = cfm._phase_viscosity(liq, 303.15, 50e5, fo.constants, fo.correlations, False)
        expected = 1.0 / (x_g / mu_g + (1.0 - x_g) / mu_l)
        mu = cfm._flash_properties(names, 50e5, 30.0, tuple(zs), "PR")[1]
        assert abs(mu - expected) / expected < 1e-9
        assert mu_g < mu < mu_l


class TestThreePhaseDensity:
    """Bulk density with free water = total mass / total volume, each HC phase
    at its own molar mass (v1.13.4)."""

    HC = ("methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane")
    Z = (0.70, 0.08, 0.07, 0.05, 0.03, 0.03)

    def _norm(self, z):
        s = sum(z)
        return tuple(x / s for x in z)

    @pytest.mark.parametrize("P, T", [(60e5, 20.0), (40e5, 30.0), (20e5, 20.0)])
    def test_trace_water_reproduces_dry_eos_density(self, P, T):
        from angelica.properties.compositional_fluid import _flash_properties
        z = self._norm(self.Z)
        dry = _flash_properties(self.HC, P, T, z, "PR")[0]
        wet = _flash_properties(self.HC + ("water",), P, T, tuple(x * (1 - 1e-9) for x in z) + (1e-9,), "PR")[0]
        assert wet == pytest.approx(dry, rel=1e-6)

    def test_mass_over_volume_with_free_water(self):
        from angelica.properties.compositional_fluid import _flash_properties, _hc_phase_split
        from angelica.properties.free_water import (
            WATER_MW,
            free_water_split,
            water_liquid_density_kg_m3,
        )
        P, T, z_w = 60e5, 20.0, 0.10
        z = self._norm(self.Z)
        zs = tuple(x * (1 - z_w) for x in z) + (z_w,)
        rho = _flash_properties(self.HC + ("water",), P, T, zs, "PR")[0]
        beta, rg, rl, mg, ml, _, _ = _hc_phase_split(self.HC, P, T, z, "PR")
        n_wl, n_wv = free_water_split(z_w, beta, 1 - z_w, T + 273.15, P)
        ng, nl = beta * (1 - z_w), (1 - beta) * (1 - z_w)
        mass = ng * mg + nl * ml + (n_wl + n_wv) * WATER_MW
        vol = (ng + n_wv) * mg / rg + nl * ml / rl + n_wl * WATER_MW / water_liquid_density_kg_m3(T + 273.15)
        assert rho == pytest.approx(mass / vol, rel=1e-9)
        assert mg < sum(a * b for a, b in zip(z, (16.04, 30.07, 44.10, 58.12, 72.15, 86.18))) < ml
