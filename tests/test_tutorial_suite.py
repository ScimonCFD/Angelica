from __future__ import annotations

import sys
import unittest
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases import (
    build_crude_oil_pipeline_case,
    build_laminar_parallel_pipes_case,
    build_steady_water_network_aggressive_elevation_case,
    build_steady_water_network_aggressive_elevation_outlet_flow_case,
    build_steady_water_network_case,
    build_steady_water_network_inlet_flow_boundary_case,
    build_steady_water_network_no_fittings_case,
    build_steady_water_network_two_flow_boundaries_case,
    build_three_reservoir_junction_case,
    build_hot_water_pipe_heat_loss_case,
    build_district_heating_branch_case,
    build_looped_network_heat_loss_case,
    build_inline_heater_case,
)
from angelica.closures import ColebrookPipeCorrelation
from angelica.core.settings import SolverSettings
from angelica.solvers import (
    SteadyIsothermalIncompressibleSolver,
    SteadyNonIsothermalIncompressibleSolver,
    NonIsothermalSolverSettings,
)


def solve(case, settings: SolverSettings):
    solver = SteadyIsothermalIncompressibleSolver(settings)
    return solver.solve(case)


class TutorialSuiteSmokeTests(unittest.TestCase):
    def test_pipe_only_case_converges(self) -> None:
        result = solve(
            build_steady_water_network_no_fittings_case(),
            SolverSettings(
                turbulent_iterations=60,
                pressure_relaxation=1.0,
                colebrook_residual_tolerance=1e-4,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-8,
            ),
        )
        self.assertTrue(result.converged)
        flows = {cf.label: cf.volumetric_flow_m3_per_h for cf in result.component_flows}
        self.assertAlmostEqual(flows["Pipe 1->2"], 6.227, delta=0.05)
        self.assertAlmostEqual(flows["Pipe 3->6"], 9.282, delta=0.05)

    def test_fittings_case_converges(self) -> None:
        result = solve(
            build_steady_water_network_case(),
            SolverSettings(
                turbulent_iterations=60,
                pressure_relaxation=1.0,
                colebrook_residual_tolerance=1e-4,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-8,
            ),
        )
        self.assertTrue(result.converged)
        flows = {cf.label: cf.volumetric_flow_m3_per_h for cf in result.component_flows}
        self.assertAlmostEqual(flows["Pipe 1->7"], 6.150, delta=0.05)
        self.assertAlmostEqual(flows["Pipe 3->6"], 9.207, delta=0.05)

    def test_elevation_case_converges(self) -> None:
        result = solve(
            build_steady_water_network_aggressive_elevation_case(),
            SolverSettings(
                turbulent_iterations=60,
                pressure_relaxation=1.0,
                colebrook_residual_tolerance=1e-4,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-8,
            ),
        )
        self.assertTrue(result.converged)
        flows = {cf.label: cf.volumetric_flow_m3_per_h for cf in result.component_flows}
        self.assertAlmostEqual(flows["Pipe 1->7"], 5.103, delta=0.05)
        self.assertAlmostEqual(flows["Pipe 3->6"], 6.388, delta=0.05)

    def test_inlet_flow_case_converges(self) -> None:
        result = solve(
            build_steady_water_network_inlet_flow_boundary_case(),
            SolverSettings(
                turbulent_iterations=60,
                pressure_relaxation=1.0,
                colebrook_residual_tolerance=1e-4,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-6,
            ),
        )
        self.assertTrue(result.converged)
        flows = {cf.label: cf.volumetric_flow_m3_per_h for cf in result.component_flows}
        self.assertAlmostEqual(flows["Pipe 1->7"], 6.150, delta=0.05)
        self.assertAlmostEqual(flows["Pipe 3->6"], 9.207, delta=0.05)

    def test_outlet_flow_case_converges(self) -> None:
        result = solve(
            build_steady_water_network_aggressive_elevation_outlet_flow_case(),
            SolverSettings(
                laminar_iterations=20,
                turbulent_iterations=120,
                pressure_relaxation=0.3,
                colebrook_residual_tolerance=1e-4,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-6,
            ),
        )
        self.assertTrue(result.converged)
        flows = {cf.label: cf.volumetric_flow_m3_per_h for cf in result.component_flows}
        self.assertAlmostEqual(flows["Pipe 1->7"], 5.103, delta=0.05)
        self.assertAlmostEqual(flows["Pipe 3->6"], 6.389, delta=0.05)

    def test_two_flow_boundaries_case_converges(self) -> None:
        result = solve(
            build_steady_water_network_two_flow_boundaries_case(),
            SolverSettings(
                turbulent_iterations=60,
                pressure_relaxation=1.0,
                colebrook_residual_tolerance=1e-4,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-7,
            ),
        )
        self.assertTrue(result.converged)
        flows = {cf.label: cf.volumetric_flow_m3_per_h for cf in result.component_flows}
        self.assertAlmostEqual(flows["Pipe 1->7"], 6.150, delta=0.05)
        self.assertAlmostEqual(flows["Pipe 3->6"], 9.207, delta=0.05)

    def test_crude_oil_pipeline_converges(self) -> None:
        result = solve(
            build_crude_oil_pipeline_case(),
            SolverSettings(
                turbulent_iterations=60,
                pressure_relaxation=1.0,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-8,
            ),
        )
        self.assertTrue(result.converged)
        flows = {cf.label: cf.volumetric_flow_m3_per_h for cf in result.component_flows}
        self.assertAlmostEqual(flows["Pipe:pipe_1"], 43.0, delta=0.1)
        self.assertAlmostEqual(flows["Pipe:pipe_2"], 26.7, delta=0.1)
        self.assertAlmostEqual(flows["Pipe:pipe_3"], 16.3, delta=0.1)

    def test_three_reservoir_junction_converges(self) -> None:
        solver = SteadyIsothermalIncompressibleSolver(
            settings=SolverSettings(
                turbulent_iterations=300,
                pressure_relaxation=0.5,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-8,
            ),
            turbulent_pipe_correlation=ColebrookPipeCorrelation(),
        )
        result = solver.solve(build_three_reservoir_junction_case())
        self.assertTrue(result.converged)
        junction_pa = result.node_pressures_pa[2]
        junction_head_m = junction_pa / (998.25 * 9.81)
        self.assertAlmostEqual(junction_head_m, 25.26, delta=0.01)
        flows = {cf.label: cf.volumetric_flow_m3_per_h for cf in result.component_flows}
        self.assertAlmostEqual(flows["Pipe 1->2"], 112.35, delta=0.1)
        self.assertAlmostEqual(flows["Pipe 2->3"], 62.36, delta=0.1)
        self.assertAlmostEqual(flows["Pipe 2->4"], 49.99, delta=0.1)

    def test_laminar_parallel_pipes_match_poiseuille(self) -> None:
        result = solve(
            build_laminar_parallel_pipes_case(),
            SolverSettings(
                laminar_iterations=50,
                turbulent_iterations=0,
                pressure_relaxation=1.0,
                pressure_correction_abs_tolerance_pa=1e-9,
                pressure_correction_rel_tolerance=1e-12,
                nodal_mass_imbalance_rel_tolerance=1e-12,
            ),
        )
        flows = {cf.label: cf.volumetric_flow_m3_per_h for cf in result.component_flows}
        self.assertAlmostEqual(flows["Pipe:pipe_1"], 0.706858347, delta=1e-6)
        self.assertAlmostEqual(flows["Pipe:pipe_2"], 0.279567712, delta=1e-6)
        self.assertAlmostEqual(flows["Pipe:pipe_3"], 1.438106675, delta=1e-6)

    def test_hot_water_pipe_heat_loss_converges(self) -> None:
        solver = SteadyNonIsothermalIncompressibleSolver()
        result = solver.solve(build_hot_water_pipe_heat_loss_case())
        self.assertTrue(result.converged)
        self.assertAlmostEqual(result.node_temperatures_c[1], 80.0, delta=0.01)
        # Outlet should be significantly cooled (well below 60 °C)
        self.assertLess(result.node_temperatures_c[2], 60.0)
        self.assertGreater(result.node_temperatures_c[2], 20.0)

    def test_district_heating_branch_converges(self) -> None:
        solver = SteadyNonIsothermalIncompressibleSolver()
        result = solver.solve(build_district_heating_branch_case())
        self.assertTrue(result.converged)
        self.assertAlmostEqual(result.node_temperatures_c[1], 85.0, delta=0.01)
        # Junction must be cooler than supply
        self.assertLess(result.node_temperatures_c[2], 85.0)
        # Long branch (node 3) must be cooler than short branch (node 4)
        self.assertLess(result.node_temperatures_c[3], result.node_temperatures_c[4])

    def test_inline_heater_converges(self) -> None:
        solver = SteadyNonIsothermalIncompressibleSolver(
            hydraulic_settings=SolverSettings(
                turbulent_iterations=200,
                pressure_relaxation=1.0,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-6,
            ),
        )
        result = solver.solve(build_inline_heater_case())
        self.assertTrue(result.converged)
        # Inlet fixed at 20 °C
        self.assertAlmostEqual(result.node_temperatures_c[1], 20.0, delta=0.01)
        # Heater adds 50 kW — outlet must be warmer than inlet
        self.assertGreater(result.node_temperatures_c[4], result.node_temperatures_c[1])
        # Energy balance: ΔT = Q / (ṁ·cₚ) — tolerance 0.5 K
        mdot = next(cf for cf in result.component_flows if "HeatSource" in cf.label).mass_flow_kg_per_s
        cp = 4182.0
        expected_delta_t = 50_000.0 / (mdot * cp)
        actual_delta_t = result.node_temperatures_c[4] - result.node_temperatures_c[1]
        self.assertAlmostEqual(actual_delta_t, expected_delta_t, delta=0.5)

    def test_looped_network_heat_loss_converges(self) -> None:
        solver = SteadyNonIsothermalIncompressibleSolver(
            hydraulic_settings=SolverSettings(
                turbulent_iterations=200,
                pressure_relaxation=0.7,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-8,
            ),
            non_isothermal_settings=NonIsothermalSolverSettings(
                temperature_relaxation=0.5,
                max_temperature_iterations=50,
            ),
        )
        result = solver.solve(build_looped_network_heat_loss_case())
        self.assertTrue(result.converged)
        self.assertAlmostEqual(result.node_temperatures_c[1], 95.0, delta=0.01)
        # Upper branch (node 3) must be cooler than the inlet
        self.assertLess(result.node_temperatures_c[3], 90.0)
        # Temperatures decrease monotonically along the main flow path
        self.assertGreater(result.node_temperatures_c[2], result.node_temperatures_c[4])
        self.assertGreater(result.node_temperatures_c[4], result.node_temperatures_c[5])
        # All nodes above ambient (10 °C)
        self.assertGreater(result.node_temperatures_c[5], 10.0)
        # Outer loop requires multiple iterations with relax=0.5
        self.assertGreater(len(result.temperature_history), 3)
