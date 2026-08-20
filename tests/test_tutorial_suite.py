from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases import (
    build_black_oil_gathering_elevation_case,
    build_crude_oil_pipeline_case,
    build_crude_oil_pipeline_thermal_case,
    build_gas_pipeline_hill_crossing_case,
    build_hilly_hot_water_network_case,
    build_laminar_parallel_pipes_case,
    build_natural_gas_pipeline_case,
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
from angelica.core.case import (
    FlowBoundary,
    InletFluidBC,
    NetworkCase,
    PressureBoundary,
    ThermalBoundary,
)
from angelica.core.components import Pipe
from angelica.core.settings import SolverSettings
from angelica.properties.black_oil import BlackOilFluid
from angelica.properties.compressible_fluid import CompressibleFluid
from angelica.properties.eos import IdealGasEOS
from angelica.solvers import (
    SteadyBlackOilSolver,
    SteadyCompressibleSolver,
    SteadyIsothermalIncompressibleSolver,
    SteadyNonIsothermalIncompressibleSolver,
    NonIsothermalSolverSettings,
)

_TUTORIALS_ROOT = (
    Path(__file__).resolve().parents[1]
    / "src" / "angelica" / "tutorials"
)


def _load_tutorial_build_case(rel_path: str):
    """Load build_case() from a tutorial run.py that has no module-level side effects."""
    spec = importlib.util.spec_from_file_location(
        "_tutorial_mod", _TUTORIALS_ROOT / rel_path
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod.build_case


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


class NonIsothermalExtendedTests(unittest.TestCase):
    """T05 and T06 of the non-isothermal incompressible tutorial series."""

    def test_t05_crude_oil_pipeline_thermal(self) -> None:
        solver = SteadyNonIsothermalIncompressibleSolver()
        result = solver.solve(build_crude_oil_pipeline_thermal_case())
        self.assertTrue(result.converged)
        self.assertAlmostEqual(result.node_temperatures_c[1], 80.0, delta=0.01)
        # Temperature decreases along the flow path
        self.assertGreater(result.node_temperatures_c[1], result.node_temperatures_c[2])
        self.assertGreater(result.node_temperatures_c[2], result.node_temperatures_c[3])
        # All nodes above ambient (15 °C)
        self.assertGreater(result.node_temperatures_c[3], 15.0)
        # Trunk throughput ≈ 9.92 kg/s
        flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
        trunk = flows.get("trunk") or flows.get("pipe_trunk") or next(iter(flows.values()))
        self.assertAlmostEqual(trunk, 9.92, delta=0.1)

    def test_t06_hilly_hot_water(self) -> None:
        solver = SteadyNonIsothermalIncompressibleSolver()
        result = solver.solve(build_hilly_hot_water_network_case())
        self.assertTrue(result.converged)
        self.assertAlmostEqual(result.node_temperatures_c[1], 85.0, delta=0.01)
        # Gravity assist: node 3 (60 m below inlet via branch_a_descending) has
        # higher pressure than the inlet (600 kPa) due to the elevation head.
        self.assertGreater(result.node_pressures_pa[3], result.node_pressures_pa[1])
        # Branch A carries exactly the prescribed outflow demand (3 kg/s)
        flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
        branch_a = flows.get("branch_a_descending", 0.0)
        self.assertAlmostEqual(branch_a, 3.0, delta=0.01)


class CompressibleTutorialTests(unittest.TestCase):
    """Steady compressible tutorials T01–T04."""

    def test_t01_natural_gas_pipeline(self) -> None:
        result = SteadyCompressibleSolver().solve(build_natural_gas_pipeline_case())
        self.assertTrue(result.converged)
        p = result.node_pressures_pa
        self.assertAlmostEqual(p[1], 800_000.0, delta=500.0)
        self.assertAlmostEqual(p[2], 667_310.0, delta=500.0)
        self.assertAlmostEqual(p[3], 500_000.0, delta=500.0)
        self.assertAlmostEqual(p[4], 500_000.0, delta=500.0)
        flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
        self.assertAlmostEqual(flows["pipe_1"], 2.8206, delta=0.01)
        self.assertAlmostEqual(flows["pipe_2"], 1.2665, delta=0.01)
        self.assertAlmostEqual(flows["pipe_3"], 1.5542, delta=0.01)
        # Mass balance
        self.assertAlmostEqual(
            flows["pipe_1"], flows["pipe_2"] + flows["pipe_3"], delta=1e-3
        )

    def test_t02_flow_bc_cross_validation(self) -> None:
        build_case = _load_tutorial_build_case(
            "steady_compressible/02_flow_bc_cross_validation/run.py"
        )
        result = SteadyCompressibleSolver().solve(build_case())
        self.assertTrue(result.converged)
        p = result.node_pressures_pa
        self.assertAlmostEqual(p[1], 800_000.0, delta=500.0)
        # Computed outlet pressures must be close to the pressure-BC values from T01
        self.assertAlmostEqual(p[3], 500_000.0, delta=1_000.0)
        self.assertAlmostEqual(p[4], 500_000.0, delta=1_000.0)
        flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
        # pipe_1 carries the sum of the two prescribed flow BCs
        self.assertAlmostEqual(flows["pipe_1"], 1.255479 + 1.540697, delta=0.01)

    def test_t03_looped_gas_heat_loss(self) -> None:
        build_case = _load_tutorial_build_case(
            "steady_compressible/03_looped_gas_pipeline_heat_loss/run.py"
        )
        result = SteadyCompressibleSolver().solve(build_case())
        self.assertTrue(result.converged)
        p = result.node_pressures_pa
        t = result.node_temperatures_c
        self.assertAlmostEqual(p[1], 700_000.0, delta=500.0)
        self.assertAlmostEqual(p[6], 500_000.0, delta=500.0)
        self.assertAlmostEqual(t[1], 40.0, delta=0.01)
        # Temperature decreases from inlet to outlet (heat loss to cold ground)
        self.assertGreater(t[1], t[6])
        # Looped network: feeder ≈ collector (mass balance)
        flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
        self.assertAlmostEqual(flows["L1-feeder"], flows["L7-collector"], delta=1e-3)

    def test_t04_hill_crossing(self) -> None:
        result = SteadyCompressibleSolver().solve(build_gas_pipeline_hill_crossing_case())
        self.assertTrue(result.converged)
        p = result.node_pressures_pa
        # Gravity + friction drop going uphill: hilltop pressure < inlet
        self.assertLess(p[2], p[1])
        # Gravity recovery going downhill past inlet elevation: node 4 > node 2
        self.assertGreater(p[4], p[2])
        flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
        # Outflow BC at branch_b_full_descent = 8 kg/s is satisfied
        self.assertAlmostEqual(flows["branch_b_full_descent"], 8.0, delta=0.01)


# ── Black-oil inline case builders ──────────────────────────────────────────

_BO_FLUID = BlackOilFluid(
    api_gravity=32.0, gas_gravity=0.65,
    gor_sc_m3_per_m3=25.0, wor_sc_m3_per_m3=0.5,
    reference_pressure_pa=5e6, reference_temperature_c=60.0,
)
_BO_ROUGHNESS = 46e-6
_BO_U = 5.0
_BO_T_AMB = 15.0


def _bo_pipe(cid, s, e, d, length):
    return Pipe(
        component_id=cid, start_node=s, end_node=e,
        diameter_m=d, length_m=length,
        absolute_roughness_m=_BO_ROUGHNESS,
        heat_transfer_coefficient_w_per_m2k=_BO_U,
        ambient_temperature_c=_BO_T_AMB,
    )


def _bo_inlet():
    return (ThermalBoundary(node_id=1, temperature_c=60.0, bc_type="fixed_temperature"),)


def _build_bo_t01() -> NetworkCase:
    return NetworkCase(
        name="T01-Black-Oil-Single-Pipe",
        fluid_model=_BO_FLUID,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=8e6),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=2e6),),
        components=(_bo_pipe("pipeline_1_2", 1, 2, 0.2, 10_000.0),),
        thermal_inlets=_bo_inlet(),
    )


def _build_bo_diamond_pipes():
    return (
        _bo_pipe("pipe_A_trunk_in", 1, 2, 0.22, 2_000.0),
        _bo_pipe("pipe_B_upper_L",  2, 3, 0.20, 5_000.0),
        _bo_pipe("pipe_C_upper_R",  3, 5, 0.18, 5_000.0),
        _bo_pipe("pipe_D_lower_L",  2, 4, 0.15, 8_000.0),
        _bo_pipe("pipe_E_lower_R",  4, 5, 0.20, 3_000.0),
        _bo_pipe("pipe_F_trunk_out", 5, 6, 0.22, 2_000.0),
    )


def _build_bo_t02() -> NetworkCase:
    return NetworkCase(
        name="T02-Black-Oil-Looped-Gathering",
        fluid_model=_BO_FLUID,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=8e6),),
        pressure_outlets=(PressureBoundary(node_id=6, pressure_pa=2e6),),
        components=_build_bo_diamond_pipes(),
        thermal_inlets=_bo_inlet(),
    )


def _build_bo_t03() -> NetworkCase:
    return NetworkCase(
        name="T03-Black-Oil-Flow-Outlet",
        fluid_model=_BO_FLUID,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=8e6),),
        pressure_outlets=(),
        flow_outlets=(FlowBoundary(node_id=6, mass_flow_kg_per_s=100.0),),
        components=_build_bo_diamond_pipes(),
        thermal_inlets=_bo_inlet(),
    )


def _build_bo_t04() -> NetworkCase:
    pipes = (
        _bo_pipe("pipe_A_trunk_in",  1, 2, 0.22, 2_000.0),
        _bo_pipe("pipe_B_upper_L",   2, 3, 0.20, 5_000.0),
        _bo_pipe("pipe_C_upper_R",   3, 5, 0.18, 5_000.0),
        _bo_pipe("pipe_D_lower_L",   2, 4, 0.18, 5_000.0),
        _bo_pipe("pipe_E_lower_R",   4, 5, 0.16, 5_000.0),
        _bo_pipe("pipe_F_trunk_out", 5, 6, 0.22, 2_000.0),
        _bo_pipe("pipe_G_satellite", 4, 7, 0.14, 2_000.0),
    )
    return NetworkCase(
        name="T04-Black-Oil-Two-Separators",
        fluid_model=_BO_FLUID,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=8e6),),
        pressure_outlets=(),
        flow_outlets=(
            FlowBoundary(node_id=6, mass_flow_kg_per_s=80.0),
            FlowBoundary(node_id=7, mass_flow_kg_per_s=20.0),
        ),
        components=pipes,
        thermal_inlets=_bo_inlet(),
    )


def _build_bo_t05() -> NetworkCase:
    default_fluid = BlackOilFluid(
        api_gravity=32.0, gas_gravity=0.65,
        gor_sc_m3_per_m3=25.0, wor_sc_m3_per_m3=0.3,
        reference_pressure_pa=5e6, reference_temperature_c=70.0,
    )
    return NetworkCase(
        name="T05-Black-Oil-Two-Reservoir-Blending",
        fluid_model=default_fluid,
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=9e6),
            PressureBoundary(node_id=2, pressure_pa=8e6),
        ),
        pressure_outlets=(PressureBoundary(node_id=4, pressure_pa=2e6),),
        components=(
            _bo_pipe("pipe_A", 1, 3, 0.18, 5_000.0),
            _bo_pipe("pipe_B", 2, 3, 0.16, 5_000.0),
            _bo_pipe("pipe_C", 3, 4, 0.22, 3_000.0),
        ),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=70.0, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=2, temperature_c=60.0, bc_type="fixed_temperature"),
        ),
        inlet_fluid_bcs=(
            InletFluidBC(node_id=1, api_gravity=32.0, gas_gravity=0.65,
                         gor_sc_m3_per_m3=25.0, wor_sc_m3_per_m3=0.3),
            InletFluidBC(node_id=2, api_gravity=22.0, gas_gravity=0.70,
                         gor_sc_m3_per_m3=10.0, wor_sc_m3_per_m3=1.5),
        ),
    )


class BlackOilTutorialTests(unittest.TestCase):
    """Steady black-oil tutorials T01–T06."""

    def test_t01_single_pipe(self) -> None:
        result = SteadyBlackOilSolver().solve(_build_bo_t01())
        self.assertTrue(result.converged)
        cf = result.component_flows[0]
        self.assertAlmostEqual(cf.mass_flow_kg_per_s, 109.79, delta=0.1)
        # Outlet at 2 MPa is below bubble point (≈5.5 MPa at 60 °C) → free gas
        pvt_out = _BO_FLUID.pvt(2e6, 60.0)
        self.assertGreater(pvt_out.holdup_gas, 0.0)

    def test_t02_looped_gathering(self) -> None:
        result = SteadyBlackOilSolver().solve(_build_bo_t02())
        self.assertTrue(result.converged)
        flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
        trunk_in = flows["pipe_A_trunk_in"]
        self.assertAlmostEqual(trunk_in, 120.42, delta=0.1)
        # Loop split: upper + lower ≈ trunk
        self.assertAlmostEqual(
            flows["pipe_B_upper_L"] + flows["pipe_D_lower_L"], trunk_in, delta=1e-3
        )
        # Both branches carry positive flow
        self.assertGreater(flows["pipe_B_upper_L"], 0.0)
        self.assertGreater(flows["pipe_D_lower_L"], 0.0)

    def test_t03_flow_outlet(self) -> None:
        result = SteadyBlackOilSolver().solve(_build_bo_t03())
        self.assertTrue(result.converged)
        flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
        # Total throughput equals the prescribed flow BC (100 kg/s)
        self.assertAlmostEqual(flows["pipe_A_trunk_in"], 100.0, delta=0.1)
        # Outlet pressure is a solver result, not prescribed
        self.assertIn(6, result.node_pressures_pa)
        self.assertGreater(result.node_pressures_pa[6], 0.0)

    def test_t04_two_separators(self) -> None:
        result = SteadyBlackOilSolver().solve(_build_bo_t04())
        self.assertTrue(result.converged)
        flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
        # Both flow BCs must be satisfied
        self.assertAlmostEqual(flows["pipe_F_trunk_out"], 80.0, delta=0.1)
        self.assertAlmostEqual(flows["pipe_G_satellite"], 20.0, delta=0.1)
        # Inlet carries the combined demand
        self.assertAlmostEqual(flows["pipe_A_trunk_in"], 100.0, delta=0.1)

    def test_t05_two_reservoir_blending(self) -> None:
        result = SteadyBlackOilSolver().solve(_build_bo_t05())
        self.assertTrue(result.converged)
        flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
        # Both reservoirs contribute positive flow
        self.assertGreater(flows["pipe_A"], 0.0)
        self.assertGreater(flows["pipe_B"], 0.0)
        # Trunk carries the blended total (mass balance)
        self.assertAlmostEqual(
            flows["pipe_C"], flows["pipe_A"] + flows["pipe_B"], delta=1e-3
        )

    def test_t06_elevation_gathering(self) -> None:
        result = SteadyBlackOilSolver().solve(build_black_oil_gathering_elevation_case())
        self.assertTrue(result.converged)
        flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
        well_a = flows["well_a_flowline"]
        well_b = flows["well_b_flowline"]
        # Well A (6 MPa, −150 m descent) out-flows Well B (9 MPa, +100 m ascent)
        # because ~1 MPa gravity gain compensates for the 3 MPa wellhead disadvantage
        self.assertGreater(well_a, well_b)
        # Manifold pressure (node 3) in expected range
        self.assertAlmostEqual(
            result.node_pressures_pa[3], 3_723_268.0, delta=50_000.0
        )
