from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Dict
from typing import List, Optional


NODE_TYPES = ("source", "sink", "junction")
LINK_COMPONENT_TYPES = ("pipe", "fitting", "pump")
DEFAULT_LIBRARY_MATERIAL = {
    "library_key": "water_liquid",
    "definition_mode": "library",
    "name": "Water",
    "density_kg_per_m3": "998.25",
    "viscosity_pa_s": "0.001",
}
DEFAULT_PRESSURE_DROP_MODEL = {
    "library_key": "colebrook_white",
    "name": "Colebrook-White",
}


@dataclass(frozen=True)
class CanvasNode:
    node_id: int
    node_type: str
    x: float
    y: float
    properties: Dict[str, str]


@dataclass(frozen=True)
class CanvasLink:
    link_id: int
    start_node_id: int
    end_node_id: int
    components: List["CanvasLinkComponent"]


@dataclass(frozen=True)
class CanvasLinkComponent:
    component_id: int
    component_type: str
    properties: Dict[str, str]


@dataclass
class CanvasScene:
    active_tool: Optional[str] = None
    nodes: List[CanvasNode] = field(default_factory=list)
    links: List[CanvasLink] = field(default_factory=list)
    material: Dict[str, str] = field(default_factory=dict)
    pressure_drop_model: Dict[str, str] = field(default_factory=lambda: dict(DEFAULT_PRESSURE_DROP_MODEL))
    solver_settings: Dict[str, float | int] = field(default_factory=dict)
    initial_node_pressures_pa: Dict[int, float] = field(default_factory=dict)
    _next_node_id: int = 1
    _next_link_id: int = 1
    _next_component_id: int = 1

    def set_active_tool(self, tool: Optional[str]) -> None:
        if tool is not None and tool not in NODE_TYPES:
            raise ValueError(f"Unsupported tool: {tool}")
        self.active_tool = tool

    def add_node(self, x: float, y: float) -> CanvasNode:
        if self.active_tool is None:
            raise ValueError("No active tool selected.")

        node = CanvasNode(
            node_id=self._next_node_id,
            node_type=self.active_tool,
            x=float(x),
            y=float(y),
            properties=self._default_properties(self.active_tool),
        )
        self.nodes.append(node)
        self._next_node_id += 1
        return node

    def add_link(self, start_node_id: int, end_node_id: int) -> CanvasLink:
        if start_node_id == end_node_id:
            raise ValueError("A node cannot be connected to itself.")
        start_node = self.get_node(start_node_id)
        end_node = self.get_node(end_node_id)
        if start_node is None or end_node is None:
            raise ValueError("Both link endpoints must already exist.")
        if self.has_link(start_node_id, end_node_id):
            raise ValueError("This connection already exists.")
        self._validate_link_capacity(start_node)
        self._validate_link_capacity(end_node)

        link = CanvasLink(
            link_id=self._next_link_id,
            start_node_id=start_node_id,
            end_node_id=end_node_id,
            components=[],
        )
        self.links.append(link)
        self._next_link_id += 1
        return link

    def move_node(self, node_id: int, x: float, y: float) -> CanvasNode:
        for index, node in enumerate(self.nodes):
            if node.node_id == node_id:
                updated_node = replace(node, x=float(x), y=float(y))
                self.nodes[index] = updated_node
                return updated_node
        raise ValueError(f"Node {node_id} does not exist.")

    def update_node_properties(self, node_id: int, properties: Dict[str, str]) -> CanvasNode:
        for index, node in enumerate(self.nodes):
            if node.node_id == node_id:
                updated_node = replace(node, properties=dict(properties))
                self.nodes[index] = updated_node
                return updated_node
        raise ValueError(f"Node {node_id} does not exist.")

    def get_link(self, link_id: int) -> Optional[CanvasLink]:
        for link in self.links:
            if link.link_id == link_id:
                return link
        return None

    def add_link_component(self, link_id: int, component_type: str) -> CanvasLink:
        if component_type not in LINK_COMPONENT_TYPES:
            raise ValueError(f"Unsupported link component: {component_type}")

        for index, link in enumerate(self.links):
            if link.link_id == link_id:
                component = CanvasLinkComponent(
                    component_id=self._next_component_id,
                    component_type=component_type,
                    properties=self._default_link_component_properties(component_type),
                )
                updated_link = replace(
                    link,
                    components=[*link.components, component],
                )
                self.links[index] = updated_link
                self._next_component_id += 1
                return updated_link
        raise ValueError(f"Link {link_id} does not exist.")

    def update_link_component_properties(
        self,
        link_id: int,
        component_id: int,
        properties: Dict[str, str],
    ) -> CanvasLink:
        for link_index, link in enumerate(self.links):
            if link.link_id != link_id:
                continue

            updated_components: List[CanvasLinkComponent] = []
            updated_component_found = False
            for component in link.components:
                if component.component_id == component_id:
                    updated_components.append(replace(component, properties=dict(properties)))
                    updated_component_found = True
                else:
                    updated_components.append(component)

            if not updated_component_found:
                raise ValueError(f"Component {component_id} does not exist in link {link_id}.")

            updated_link = replace(link, components=updated_components)
            self.links[link_index] = updated_link
            return updated_link

        raise ValueError(f"Link {link_id} does not exist.")

    def get_node(self, node_id: int) -> Optional[CanvasNode]:
        for node in self.nodes:
            if node.node_id == node_id:
                return node
        return None

    def has_link(self, start_node_id: int, end_node_id: int) -> bool:
        for link in self.links:
            if (
                link.start_node_id == start_node_id
                and link.end_node_id == end_node_id
            ) or (
                link.start_node_id == end_node_id
                and link.end_node_id == start_node_id
            ):
                return True
        return False

    def connection_count(self, node_id: int) -> int:
        count = 0
        for link in self.links:
            if link.start_node_id == node_id or link.end_node_id == node_id:
                count += 1
        return count

    def remove_node(self, node_id: int) -> list[int]:
        """Remove a node and all links touching it. Returns the removed link IDs."""
        if self.get_node(node_id) is None:
            raise ValueError(f"Node {node_id} does not exist.")
        self.nodes = [n for n in self.nodes if n.node_id != node_id]
        removed_link_ids = [
            link.link_id
            for link in self.links
            if link.start_node_id == node_id or link.end_node_id == node_id
        ]
        self.links = [
            link for link in self.links
            if link.start_node_id != node_id and link.end_node_id != node_id
        ]
        self.initial_node_pressures_pa.pop(node_id, None)
        return removed_link_ids

    def clear(self) -> None:
        self.nodes.clear()
        self.links.clear()
        self.material = {}
        self.pressure_drop_model = dict(DEFAULT_PRESSURE_DROP_MODEL)
        self.solver_settings.clear()
        self.initial_node_pressures_pa.clear()
        self._next_node_id = 1
        self._next_link_id = 1
        self._next_component_id = 1
        self.active_tool = None

    def update_material(self, material: Dict[str, str]) -> None:
        self.material = dict(material)

    def update_pressure_drop_model(self, pressure_drop_model: Dict[str, str]) -> None:
        self.pressure_drop_model = dict(pressure_drop_model)

    def update_solver_settings(self, settings: Dict[str, float | int | str]) -> None:
        updated_settings = dict(self.solver_settings)
        updated_settings.update(settings)
        self.solver_settings = updated_settings

    @staticmethod
    def _default_properties(node_type: str) -> Dict[str, str]:
        if node_type == "source":
            return {
                "condition_type": "pressure",
                "pressure": "",
                "flow": "",
            }
        if node_type == "sink":
            return {
                "condition_type": "pressure",
                "pressure": "",
                "flow": "",
            }
        return {
            "label": "",
        }

    @staticmethod
    def _default_link_component_properties(component_type: str) -> Dict[str, str]:
        if component_type == "pipe":
            return {
                "length_m": "",
                "diameter_m": "",
                "height_change_m": "0.0",
                "roughness_m": "0.000045",
                "hazen_williams_c": "130.0",
                "num_segments": "1",
            }
        if component_type == "pump":
            return {
                "diameter_m": "",
                "curve_points_q_head": "",
            }
        return {
            "diameter_m": "",
            "fitting_mode": "manual",
            "fitting_preset": "regular_90_flanged",
            "loss_coefficient": "1.5",
        }

    def _validate_link_capacity(self, node: CanvasNode) -> None:
        if node.node_type not in {"source", "sink"}:
            return
        if self.connection_count(node.node_id) >= 2:
            raise ValueError(
                f"{node.node_type.capitalize()} #{node.node_id} cannot have more than two connections."
            )
