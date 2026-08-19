from __future__ import annotations

import sys
import unittest
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

TUTORIALS_ROOT = SRC_ROOT / "angelica" / "tutorials"

from angelica.gui.io import build_network_case_from_scene, build_solver_from_scene, load_scene_from_file
from angelica.gui.model import CanvasLink, CanvasLinkComponent, CanvasNode, CanvasScene


class GuiIoTests(unittest.TestCase):
    @staticmethod
    def _pipe_only_case_path() -> Path:
        return (
            TUTORIALS_ROOT
            / "steady_isothermal_incompressible"
            / "01_pipe_only"
            / "pipe_only.gui.json"
        )

    def test_load_pipe_only_tutorial_scene(self) -> None:
        scene = load_scene_from_file(self._pipe_only_case_path())

        self.assertEqual(len(scene.nodes), 6)
        self.assertEqual(len(scene.links), 6)
        self.assertEqual(scene.material["library_key"], "water_liquid")
        self.assertEqual(scene.pressure_drop_model["library_key"], "colebrook_white")
        self.assertEqual(scene.get_node(1).properties["pressure"], "251300.0")
        self.assertEqual(scene.get_link(1).components[0].component_type, "pipe")

    def test_pipe_only_tutorial_scene_builds_and_converges(self) -> None:
        scene = load_scene_from_file(self._pipe_only_case_path())
        case = build_network_case_from_scene(scene)
        result = build_solver_from_scene(scene).solve(case)

        self.assertEqual(case.name, "Pipe-only network")
        self.assertEqual(len(case.components), 6)
        self.assertTrue(result.converged)

    def test_build_network_case_uses_scene_material(self) -> None:
        scene = load_scene_from_file(self._pipe_only_case_path())
        scene.update_material(
            {
                "library_key": "water_liquid",
                "name": "Water",
                "density_kg_per_m3": "1000.0",
                "viscosity_pa_s": "0.002",
            }
        )

        case = build_network_case_from_scene(scene)

        self.assertEqual(case.fluid_model.density_kg_per_m3, 1000.0)
        self.assertEqual(case.fluid_model.viscosity_pa_s, 0.002)

    def test_build_network_case_requires_material_definition(self) -> None:
        scene = load_scene_from_file(self._pipe_only_case_path())
        scene.material = {}

        with self.assertRaises(ValueError):
            build_network_case_from_scene(scene)

    def test_build_network_case_uses_fitting_preset_loss_coefficient(self) -> None:
        scene = load_scene_from_file(
            TUTORIALS_ROOT
            / "steady_isothermal_incompressible"
            / "02_fittings_no_elevation"
            / "fittings_no_elevation.gui.json"
        )
        first_link = scene.get_link(1)
        assert first_link is not None
        fitting_component = first_link.components[1]
        scene.update_link_component_properties(
            first_link.link_id,
            fitting_component.component_id,
            {
                "diameter_m": "0.05",
                "fitting_mode": "preset",
                "fitting_preset": "globe_valve_fully_open",
                "loss_coefficient": "10.0",
            },
        )

        case = build_network_case_from_scene(scene)
        fitting = next(component for component in case.components if type(component).__name__ == "Fitting")

        self.assertEqual(fitting.loss_coefficient, 10.0)

    def test_build_network_case_parses_pump_curve_table(self) -> None:
        scene = load_scene_from_file(self._pipe_only_case_path())
        first_link = scene.get_link(1)
        assert first_link is not None
        scene.add_link_component(first_link.link_id, "pump")
        pump_component = scene.get_link(first_link.link_id).components[-1]
        scene.update_link_component_properties(
            first_link.link_id,
            pump_component.component_id,
            {
                "diameter_m": "0.15",
                "curve_points_q_head": "0, 80\n200, 60\n400, 0",
            },
        )

        case = build_network_case_from_scene(scene)
        pump = next(component for component in case.components if type(component).__name__ == "Pump")

        self.assertEqual(pump.diameter_m, 0.15)
        self.assertEqual(pump.curve_points_q_head, ((0.0, 80.0), (200.0, 60.0), (400.0, 0.0)))

    def test_all_gui_tutorial_scenes_build_and_converge(self) -> None:
        tutorial_root = TUTORIALS_ROOT / "steady_isothermal_incompressible"
        gui_paths = sorted(tutorial_root.glob("*/*.gui.json"))

        self.assertGreaterEqual(len(gui_paths), 7)

        for gui_path in gui_paths:
            with self.subTest(gui_path=gui_path.name):
                scene = load_scene_from_file(gui_path)
                case = build_network_case_from_scene(scene)
                result = build_solver_from_scene(scene).solve(case)
                self.assertTrue(result.converged, msg=f"{gui_path} did not converge")

    def test_all_non_isothermal_gui_tutorial_scenes_build_and_converge(self) -> None:
        tutorial_root = TUTORIALS_ROOT / "steady_non_isothermal_incompressible"
        gui_paths = sorted(tutorial_root.glob("*/*.gui.json"))

        self.assertGreaterEqual(len(gui_paths), 4)

        for gui_path in gui_paths:
            with self.subTest(gui_path=gui_path.name):
                scene = load_scene_from_file(gui_path)
                case = build_network_case_from_scene(scene)
                result = build_solver_from_scene(scene).solve(case)
                self.assertTrue(result.converged, msg=f"{gui_path} did not converge")
                self.assertTrue(
                    result.node_temperatures_c,
                    msg=f"{gui_path} produced no temperature results",
                )

    def test_all_compressible_gui_tutorial_scenes_build_and_converge(self) -> None:
        tutorial_root = TUTORIALS_ROOT / "steady_compressible"
        gui_paths = sorted(tutorial_root.glob("*/*.gui.json"))

        self.assertGreaterEqual(len(gui_paths), 1)

        for gui_path in gui_paths:
            with self.subTest(gui_path=gui_path.name):
                scene = load_scene_from_file(gui_path)
                case = build_network_case_from_scene(scene)
                result = build_solver_from_scene(scene).solve(case)
                self.assertTrue(result.converged, msg=f"{gui_path} did not converge")

    def test_all_black_oil_gui_tutorial_scenes_build_and_converge(self) -> None:
        tutorial_root = TUTORIALS_ROOT / "steady_black_oil"
        gui_paths = sorted(tutorial_root.glob("*/*.gui.json"))

        self.assertGreaterEqual(len(gui_paths), 5)

        for gui_path in gui_paths:
            with self.subTest(gui_path=gui_path.name):
                scene = load_scene_from_file(gui_path)
                case = build_network_case_from_scene(scene)
                result = build_solver_from_scene(scene).solve(case)
                self.assertTrue(result.converged, msg=f"{gui_path} did not converge")

    def test_build_solver_uses_supported_default_pressure_drop_model(self) -> None:
        scene = load_scene_from_file(self._pipe_only_case_path())

        solver = build_solver_from_scene(scene)

        self.assertEqual(scene.pressure_drop_model["library_key"], "colebrook_white")
        self.assertEqual(type(solver.turbulent_pipe_correlation).__name__, "ColebrookPipeCorrelation")
        self.assertIsNone(solver.settings.laminar_iterations)
        self.assertEqual(solver.settings.turbulent_iterations, 60)
        self.assertEqual(solver.settings.colebrook_friction_strategy, "transformed")
        self.assertEqual(solver.settings.friction_factor_method, "newton")
        self.assertEqual(solver.settings.friction_factor_max_iterations, 50)
        self.assertEqual(solver.settings.velocity_loop_method, "secant")
        self.assertEqual(solver.settings.velocity_loop_max_iterations, 50)

    def test_build_solver_reads_explicit_numerics_settings(self) -> None:
        scene = load_scene_from_file(self._pipe_only_case_path())
        scene.update_solver_settings(
            {
                "laminar_iterations": 14,
                "turbulent_iterations": 140,
                "pressure_relaxation": 0.5,
                "colebrook_friction_strategy": "direct",
                "friction_factor_method": "newton",
                "friction_factor_max_iterations": 80,
                "velocity_loop_method": "secant",
                "velocity_loop_max_iterations": 120,
            }
        )

        solver = build_solver_from_scene(scene)

        self.assertEqual(solver.settings.laminar_iterations, 14)
        self.assertEqual(solver.settings.turbulent_iterations, 140)
        self.assertEqual(solver.settings.pressure_relaxation, 0.5)
        self.assertEqual(solver.settings.colebrook_friction_strategy, "direct")
        self.assertEqual(solver.settings.friction_factor_method, "newton")
        self.assertEqual(solver.settings.friction_factor_max_iterations, 80)
        self.assertEqual(solver.settings.velocity_loop_method, "secant")
        self.assertEqual(solver.settings.velocity_loop_max_iterations, 120)

    def test_build_solver_supports_hazen_williams_pressure_drop_model(self) -> None:
        scene = load_scene_from_file(self._pipe_only_case_path())
        scene.update_pressure_drop_model(
            {
                "library_key": "hazen_williams",
                "name": "Hazen-Williams",
            }
        )

        solver = build_solver_from_scene(scene)

        self.assertEqual(scene.pressure_drop_model["library_key"], "hazen_williams")
        self.assertEqual(type(solver.turbulent_pipe_correlation).__name__, "HazenWilliamsPipeCorrelation")

    def test_hazen_williams_single_pipe_matches_analytical(self) -> None:
        import math

        D = 0.05      # m
        L = 100.0     # m
        C = 130.0     # HW roughness coefficient
        P_in = 200_000.0   # Pa
        P_out = 100_000.0  # Pa
        rho = 998.25  # kg/m³

        h_f = (P_in - P_out) / (rho * 9.81)
        R = 10.67 * L / (C ** 1.852 * D ** 4.871)
        q_analytical = (h_f / R) ** (1.0 / 1.852)
        mdot_analytical = rho * q_analytical

        pipe_props = {
            "length_m": str(L),
            "diameter_m": str(D),
            "height_change_m": "0.0",
            "roughness_m": "0.000045",
            "hazen_williams_c": str(C),
            "num_segments": "1",
        }
        scene = CanvasScene()
        scene.nodes = [
            CanvasNode(1, "source", 0.0, 0.0, {"condition_type": "pressure", "pressure": str(P_in), "flow": ""}),
            CanvasNode(2, "sink", 100.0, 0.0, {"condition_type": "pressure", "pressure": str(P_out), "flow": ""}),
        ]
        scene.links = [
            CanvasLink(
                link_id=1,
                start_node_id=1,
                end_node_id=2,
                components=[CanvasLinkComponent(1, "pipe", pipe_props)],
            )
        ]
        scene.material = {"density_kg_per_m3": str(rho), "viscosity_pa_s": "0.001"}
        scene.update_pressure_drop_model({"library_key": "hazen_williams", "name": "Hazen-Williams"})

        case = build_network_case_from_scene(scene)
        result = build_solver_from_scene(scene).solve(case)

        self.assertTrue(result.converged)
        pipe_mdot = result.component_flows[0].mass_flow_kg_per_s
        relative_error = abs(pipe_mdot - mdot_analytical) / mdot_analytical
        self.assertLess(relative_error, 1e-4, msg=f"HW flow error {relative_error:.2e} exceeds 0.01%")

    # ------------------------------------------------------------------
    # Junction topology validation
    # ------------------------------------------------------------------

    @staticmethod
    def _minimal_scene_with_junction_links(link_pairs: list[tuple[int, int]]) -> CanvasScene:
        """Build a 3-node scene (source-1, junction-2, sink-3) with given link topology."""
        pipe_props = {
            "length_m": "100",
            "diameter_m": "0.05",
            "height_change_m": "0.0",
            "roughness_m": "0.000045",
            "hazen_williams_c": "130.0",
            "num_segments": "1",
        }
        scene = CanvasScene()
        scene.nodes = [
            CanvasNode(1, "source", 0.0, 0.0, {"condition_type": "pressure", "pressure": "200000", "flow": ""}),
            CanvasNode(2, "junction", 100.0, 0.0, {"label": ""}),
            CanvasNode(3, "sink", 200.0, 0.0, {"condition_type": "pressure", "pressure": "100000", "flow": ""}),
        ]
        scene.links = [
            CanvasLink(
                link_id=idx + 1,
                start_node_id=start,
                end_node_id=end,
                components=[CanvasLinkComponent(idx + 1, "pipe", pipe_props)],
            )
            for idx, (start, end) in enumerate(link_pairs)
        ]
        scene.material = {"density_kg_per_m3": "998.25", "viscosity_pa_s": "0.001"}
        return scene

    def test_heat_source_in_isothermal_mode_is_rejected(self) -> None:
        scene = CanvasScene()
        scene.nodes = [
            CanvasNode(1, "source", 0.0, 0.0, {"condition_type": "pressure", "pressure": "200000", "flow": ""}),
            CanvasNode(2, "sink", 100.0, 0.0, {"condition_type": "pressure", "pressure": "100000", "flow": ""}),
        ]
        scene.links = [
            CanvasLink(
                link_id=1,
                start_node_id=1,
                end_node_id=2,
                components=[CanvasLinkComponent(1, "heat_source", {"diameter_m": "0.05", "power_w": "1000"})],
            )
        ]
        scene.material = {"density_kg_per_m3": "998.25", "viscosity_pa_s": "0.001"}
        scene.physics_mode = "isothermal"

        with self.assertRaises(ValueError) as ctx:
            build_network_case_from_scene(scene)
        self.assertIn("Heat Source", str(ctx.exception))
        self.assertIn("Non-isothermal", str(ctx.exception))

    def test_isolated_junction_is_rejected(self) -> None:
        # junction(2) has no connections at all
        scene = self._minimal_scene_with_junction_links([(1, 3)])
        with self.assertRaises(ValueError) as ctx:
            build_network_case_from_scene(scene)
        self.assertIn("not connected to any pipe", str(ctx.exception))

    def test_dead_end_junction_is_rejected(self) -> None:
        # junction(2) has exactly one connection — dead end, mass cannot be conserved
        scene = self._minimal_scene_with_junction_links([(1, 2), (1, 3)])
        with self.assertRaises(ValueError) as ctx:
            build_network_case_from_scene(scene)
        self.assertIn("dead end", str(ctx.exception))

    def test_junction_with_two_connections_is_accepted(self) -> None:
        # junction(2) has one connection in and one out — valid
        scene = self._minimal_scene_with_junction_links([(1, 2), (2, 3)])
        case = build_network_case_from_scene(scene)
        self.assertEqual(len(case.components), 2)
