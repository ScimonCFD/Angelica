from __future__ import annotations

import sys
import unittest
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from netSim.gui.io import build_network_case_from_scene, build_solver_from_scene, load_scene_from_file


class GuiIoTests(unittest.TestCase):
    @staticmethod
    def _pipe_only_case_path() -> Path:
        return (
            Path(__file__).resolve().parents[1]
            / "tutorials"
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

        self.assertEqual(case.name, "GUI scene")
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
            Path(__file__).resolve().parents[1]
            / "tutorials"
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
        tutorial_root = (
            Path(__file__).resolve().parents[1]
            / "tutorials"
            / "steady_isothermal_incompressible"
        )
        gui_paths = sorted(tutorial_root.glob("*/*.gui.json"))

        self.assertGreaterEqual(len(gui_paths), 7)

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
        self.assertEqual(solver.settings.velocity_loop_method, "fixed_point")
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
