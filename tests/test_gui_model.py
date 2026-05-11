from __future__ import annotations

import sys
import unittest
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from netSim.gui.model import CanvasScene


class GuiModelTests(unittest.TestCase):
    def test_add_node_requires_tool(self) -> None:
        scene = CanvasScene()
        with self.assertRaises(ValueError):
            scene.add_node(10, 20)

    def test_add_node_uses_active_tool(self) -> None:
        scene = CanvasScene()
        scene.set_active_tool("source")
        node = scene.add_node(10, 20)

        self.assertEqual(node.node_type, "source")
        self.assertEqual(node.node_id, 1)
        self.assertEqual(len(scene.nodes), 1)
        self.assertEqual(node.properties["condition_type"], "pressure")
        self.assertEqual(node.properties["pressure"], "")
        self.assertEqual(node.properties["flow"], "")

    def test_clear_resets_scene(self) -> None:
        scene = CanvasScene()
        scene.set_active_tool("junction")
        scene.add_node(50, 60)
        scene.update_material(
            {
                "library_key": "water_liquid",
                "name": "Water",
                "density_kg_per_m3": "1000.0",
                "viscosity_pa_s": "0.002",
            }
        )
        scene.clear()

        self.assertEqual(scene.nodes, [])
        self.assertIsNone(scene.active_tool)
        self.assertEqual(scene.material, {})

    def test_scene_starts_with_default_material(self) -> None:
        scene = CanvasScene()

        self.assertEqual(scene.material, {})

    def test_scene_starts_with_default_pressure_drop_model(self) -> None:
        scene = CanvasScene()

        self.assertEqual(scene.pressure_drop_model["library_key"], "colebrook_white")
        self.assertEqual(scene.pressure_drop_model["name"], "Colebrook-White")

    def test_update_solver_settings_merges_values(self) -> None:
        scene = CanvasScene()

        scene.update_solver_settings({"turbulent_iterations": 80})
        scene.update_solver_settings(
            {
                "laminar_iterations": 12,
                "pressure_relaxation": 0.5,
                "colebrook_friction_strategy": "direct",
                "friction_factor_method": "newton",
                "friction_factor_max_iterations": 80,
                "velocity_loop_method": "secant",
                "velocity_loop_max_iterations": 120,
            }
        )

        self.assertEqual(scene.solver_settings["laminar_iterations"], 12)
        self.assertEqual(scene.solver_settings["turbulent_iterations"], 80)
        self.assertEqual(scene.solver_settings["pressure_relaxation"], 0.5)
        self.assertEqual(scene.solver_settings["colebrook_friction_strategy"], "direct")
        self.assertEqual(scene.solver_settings["friction_factor_method"], "newton")
        self.assertEqual(scene.solver_settings["friction_factor_max_iterations"], 80)
        self.assertEqual(scene.solver_settings["velocity_loop_method"], "secant")
        self.assertEqual(scene.solver_settings["velocity_loop_max_iterations"], 120)

    def test_add_link_connects_two_existing_nodes(self) -> None:
        scene = CanvasScene()
        scene.set_active_tool("source")
        node_1 = scene.add_node(10, 20)
        scene.set_active_tool("sink")
        node_2 = scene.add_node(100, 120)

        link = scene.add_link(node_1.node_id, node_2.node_id)

        self.assertEqual(link.link_id, 1)
        self.assertEqual(len(scene.links), 1)
        self.assertEqual(link.components, [])

    def test_duplicate_link_is_rejected(self) -> None:
        scene = CanvasScene()
        scene.set_active_tool("source")
        node_1 = scene.add_node(10, 20)
        scene.set_active_tool("sink")
        node_2 = scene.add_node(100, 120)
        scene.add_link(node_1.node_id, node_2.node_id)

        with self.assertRaises(ValueError):
            scene.add_link(node_2.node_id, node_1.node_id)

    def test_update_node_properties_replaces_properties(self) -> None:
        scene = CanvasScene()
        scene.set_active_tool("sink")
        node = scene.add_node(10, 20)

        updated = scene.update_node_properties(
            node.node_id,
            {"condition_type": "flow", "pressure": "", "flow": "2.5"},
        )

        self.assertEqual(updated.properties["condition_type"], "flow")
        self.assertEqual(updated.properties["flow"], "2.5")

    def test_add_link_component_updates_link(self) -> None:
        scene = CanvasScene()
        scene.set_active_tool("source")
        node_1 = scene.add_node(10, 20)
        scene.set_active_tool("sink")
        node_2 = scene.add_node(100, 120)
        link = scene.add_link(node_1.node_id, node_2.node_id)

        updated_link = scene.add_link_component(link.link_id, "pipe")

        self.assertEqual(updated_link.components[0].component_type, "pipe")
        self.assertEqual(updated_link.components[0].properties["height_change_m"], "0.0")
        self.assertEqual(updated_link.components[0].properties["num_segments"], "1")

    def test_update_link_component_properties(self) -> None:
        scene = CanvasScene()
        scene.set_active_tool("source")
        node_1 = scene.add_node(10, 20)
        scene.set_active_tool("sink")
        node_2 = scene.add_node(100, 120)
        link = scene.add_link(node_1.node_id, node_2.node_id)
        link = scene.add_link_component(link.link_id, "fitting")

        updated_link = scene.update_link_component_properties(
            link.link_id,
            link.components[0].component_id,
            {"loss_coefficient": "2.0"},
        )

        self.assertEqual(updated_link.components[0].properties["loss_coefficient"], "2.0")

    def test_source_cannot_have_more_than_two_connections(self) -> None:
        scene = CanvasScene()
        scene.set_active_tool("source")
        source = scene.add_node(10, 20)

        scene.set_active_tool("junction")
        j1 = scene.add_node(40, 20)
        scene.set_active_tool("junction")
        j2 = scene.add_node(70, 20)
        scene.set_active_tool("junction")
        j3 = scene.add_node(100, 20)

        scene.add_link(source.node_id, j1.node_id)
        scene.add_link(source.node_id, j2.node_id)

        with self.assertRaises(ValueError):
            scene.add_link(source.node_id, j3.node_id)

    def test_sink_cannot_have_more_than_two_connections(self) -> None:
        scene = CanvasScene()
        scene.set_active_tool("sink")
        sink = scene.add_node(10, 20)

        scene.set_active_tool("junction")
        j1 = scene.add_node(40, 20)
        scene.set_active_tool("junction")
        j2 = scene.add_node(70, 20)
        scene.set_active_tool("junction")
        j3 = scene.add_node(100, 20)

        scene.add_link(sink.node_id, j1.node_id)
        scene.add_link(sink.node_id, j2.node_id)

        with self.assertRaises(ValueError):
            scene.add_link(j3.node_id, sink.node_id)

    def test_junction_can_have_more_than_two_connections(self) -> None:
        scene = CanvasScene()
        scene.set_active_tool("junction")
        junction = scene.add_node(50, 50)

        scene.set_active_tool("source")
        s1 = scene.add_node(10, 20)
        scene.set_active_tool("source")
        s2 = scene.add_node(10, 80)
        scene.set_active_tool("sink")
        sink = scene.add_node(100, 50)

        scene.add_link(junction.node_id, s1.node_id)
        scene.add_link(junction.node_id, s2.node_id)
        scene.add_link(junction.node_id, sink.node_id)

        self.assertEqual(scene.connection_count(junction.node_id), 3)
