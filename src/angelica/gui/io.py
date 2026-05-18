from __future__ import annotations

import json
import math
from pathlib import Path

from angelica.core.case import FlowBoundary, NetworkCase, PressureBoundary
from angelica.core.components import FITTING_PRESET_LIBRARY, Fitting, Pipe, Pump
from angelica.core.settings import SolverSettings
from angelica.closures import ColebrookPipeCorrelation, HazenWilliamsPipeCorrelation
from angelica.properties.single_component import SingleComponentFluid
from angelica.solvers import SteadyIsothermalIncompressibleSolver

from .model import (
    CanvasLink,
    CanvasLinkComponent,
    CanvasNode,
    CanvasScene,
    DEFAULT_PRESSURE_DROP_MODEL,
)


def scene_from_dict(data: dict) -> CanvasScene:
    scene = CanvasScene()
    scene.nodes = [
        CanvasNode(
            node_id=int(node["node_id"]),
            node_type=str(node["node_type"]),
            x=float(node["x"]),
            y=float(node["y"]),
            properties=dict(node.get("properties", {})),
        )
        for node in data.get("nodes", [])
    ]
    scene.links = [
        CanvasLink(
            link_id=int(link["link_id"]),
            start_node_id=int(link["start_node_id"]),
            end_node_id=int(link["end_node_id"]),
            components=[
                CanvasLinkComponent(
                    component_id=int(component["component_id"]),
                    component_type=str(component["component_type"]),
                    properties=dict(component.get("properties", {})),
                )
                for component in link.get("components", [])
            ],
        )
        for link in data.get("links", [])
    ]

    scene._next_node_id = max((node.node_id for node in scene.nodes), default=0) + 1
    scene._next_link_id = max((link.link_id for link in scene.links), default=0) + 1
    scene._next_component_id = (
        max(
            (
                component.component_id
                for link in scene.links
                for component in link.components
            ),
            default=0,
        )
        + 1
    )
    scene.active_tool = None
    scene.material = {
        key: str(value) for key, value in data.get("material", {}).items()
    }
    scene.pressure_drop_model = dict(DEFAULT_PRESSURE_DROP_MODEL)
    scene.pressure_drop_model.update(
        {key: str(value) for key, value in data.get("pressure_drop_model", {}).items()}
    )
    scene.solver_settings = dict(data.get("solver_settings", {}))
    scene.initial_node_pressures_pa = {
        int(node_id): float(value)
        for node_id, value in data.get("initial_node_pressures_pa", {}).items()
    }
    return scene


def load_scene_from_file(path: str | Path) -> CanvasScene:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return scene_from_dict(data)


def build_network_case_from_scene(scene: CanvasScene) -> NetworkCase:
    if not scene.nodes:
        raise ValueError("The scene is empty. Add nodes before running the simulation.")
    if not scene.links:
        raise ValueError("The scene has no links. Add at least one connection before running.")
    if not scene.material:
        raise ValueError("No material is defined. Use Material -> Define Material before running.")
    if not scene.material.get("density_kg_per_m3", "").strip():
        raise ValueError("The material is missing density_kg_per_m3.")
    if not scene.material.get("viscosity_pa_s", "").strip():
        raise ValueError("The material is missing viscosity_pa_s.")

    pressure_inlets: list[PressureBoundary] = []
    pressure_outlets: list[PressureBoundary] = []
    flow_inlets: list[FlowBoundary] = []
    flow_outlets: list[FlowBoundary] = []

    for node in scene.nodes:
        if node.node_type == "junction":
            continue

        condition_type = node.properties.get("condition_type", "pressure")
        if condition_type == "pressure":
            value_text = node.properties.get("pressure", "").strip()
            if not value_text:
                raise ValueError(
                    f"{node.node_type.capitalize()} #{node.node_id} is missing a pressure value."
                )
            boundary = PressureBoundary(node_id=node.node_id, pressure_pa=float(value_text))
            if node.node_type == "source":
                pressure_inlets.append(boundary)
            else:
                pressure_outlets.append(boundary)
        elif condition_type == "flow":
            value_text = node.properties.get("flow", "").strip()
            if not value_text:
                raise ValueError(
                    f"{node.node_type.capitalize()} #{node.node_id} is missing a flow value."
                )
            boundary = FlowBoundary(node_id=node.node_id, mass_flow_kg_per_s=float(value_text))
            if node.node_type == "source":
                flow_inlets.append(boundary)
            else:
                flow_outlets.append(boundary)
        else:
            raise ValueError(
                f"{node.node_type.capitalize()} #{node.node_id} has unsupported boundary type '{condition_type}'."
            )

    components = []
    next_internal_node_id = max(node.node_id for node in scene.nodes) + 1
    for link in scene.links:
        if not link.components:
            raise ValueError(
                f"Connection #{link.link_id} must contain at least one component to run."
            )

        current_start = link.start_node_id
        for component_index, component in enumerate(link.components):
            is_last_component = component_index == len(link.components) - 1
            current_end = link.end_node_id if is_last_component else next_internal_node_id
            if not is_last_component:
                next_internal_node_id += 1

            if component.component_type == "pipe":
                diameter = _required_float(component, "diameter_m", link.link_id)
                length = _required_float(component, "length_m", link.link_id)
                roughness = _optional_float(component, "roughness_m", default=0.000045)
                hazen_williams_c = _optional_float(
                    component,
                    "hazen_williams_c",
                    default=130.0,
                )
                height_change = _optional_float(component, "height_change_m", default=0.0)
                components.append(
                    Pipe(
                        start_node=current_start,
                        end_node=current_end,
                        diameter_m=diameter,
                        length_m=length,
                        absolute_roughness_m=roughness,
                        hazen_williams_c=hazen_williams_c,
                        height_change_m=height_change,
                        component_id=f"link_{link.link_id}_pipe_{component.component_id}",
                    )
                )
            elif component.component_type == "fitting":
                diameter = _required_float(component, "diameter_m", link.link_id)
                fitting_mode = component.properties.get("fitting_mode", "manual").strip() or "manual"
                if fitting_mode == "preset":
                    preset_key = component.properties.get("fitting_preset", "").strip()
                    if preset_key not in FITTING_PRESET_LIBRARY:
                        raise ValueError(
                            f"Fitting in connection #{link.link_id} has unsupported preset '{preset_key}'."
                        )
                    loss_coefficient = float(FITTING_PRESET_LIBRARY[preset_key]["loss_coefficient"])
                else:
                    loss_coefficient = _required_float(component, "loss_coefficient", link.link_id)
                components.append(
                    Fitting(
                        start_node=current_start,
                        end_node=current_end,
                        diameter_m=diameter,
                        loss_coefficient=loss_coefficient,
                        component_id=f"link_{link.link_id}_fitting_{component.component_id}",
                    )
                )
            elif component.component_type == "pump":
                diameter = _required_float(component, "diameter_m", link.link_id)
                components.append(
                    Pump(
                        start_node=current_start,
                        end_node=current_end,
                        diameter_m=diameter,
                        curve_points_q_head=_parse_pump_curve_points(component, link.link_id),
                        component_id=f"link_{link.link_id}_pump_{component.component_id}",
                    )
                )
            else:
                raise ValueError(
                    f"Unsupported component type '{component.component_type}' in connection #{link.link_id}."
                )

            current_start = current_end

    visible_node_ids = {node.node_id for node in scene.nodes}
    all_node_ids = visible_node_ids.union(range(max(visible_node_ids) + 1, next_internal_node_id))

    return NetworkCase(
        name="GUI scene",
        fluid_model=SingleComponentFluid(
            density_kg_per_m3=float(scene.material["density_kg_per_m3"]),
            viscosity_pa_s=float(scene.material["viscosity_pa_s"]),
        ),
        pressure_inlets=tuple(pressure_inlets),
        pressure_outlets=tuple(pressure_outlets),
        flow_inlets=tuple(flow_inlets),
        flow_outlets=tuple(flow_outlets),
        components=tuple(components),
        node_ids=tuple(sorted(all_node_ids)),
        initial_node_pressures_pa=dict(scene.initial_node_pressures_pa),
    )


def build_solver_from_scene(scene: CanvasScene) -> SteadyIsothermalIncompressibleSolver:
    pressure_drop_model_key = scene.pressure_drop_model.get("library_key", "")
    if pressure_drop_model_key == "colebrook_white":
        turbulent_pipe_correlation = ColebrookPipeCorrelation()
    elif pressure_drop_model_key == "hazen_williams":
        turbulent_pipe_correlation = HazenWilliamsPipeCorrelation()
    else:
        raise ValueError(
            f"Unsupported pipe pressure-drop model '{pressure_drop_model_key}'."
        )

    settings = SolverSettings(**scene.solver_settings) if scene.solver_settings else SolverSettings()

    if scene.solver_settings:
        return SteadyIsothermalIncompressibleSolver(
            settings=settings,
            turbulent_pipe_correlation=turbulent_pipe_correlation,
        )
    return SteadyIsothermalIncompressibleSolver(
        settings=settings,
        turbulent_pipe_correlation=turbulent_pipe_correlation,
    )


def _required_float(
    component: CanvasLinkComponent,
    field_name: str,
    link_id: int,
) -> float:
    value_text = component.properties.get(field_name, "").strip()
    if not value_text:
        raise ValueError(
            f"{component.component_type.capitalize()} in connection #{link_id} is missing '{field_name}'."
        )
    return _parse_float(value_text)


def _optional_float(
    component: CanvasLinkComponent,
    field_name: str,
    default: float,
) -> float:
    value_text = component.properties.get(field_name, "").strip()
    if not value_text:
        return default
    return _parse_float(value_text)


def _parse_float(value_text: str) -> float:
    lowered = value_text.strip().lower()
    if lowered in {"inf", "+inf", "infinity", "+infinity"}:
        return math.inf
    if lowered in {"-inf", "-infinity"}:
        return -math.inf
    return float(value_text)


def _parse_pump_curve_points(
    component: CanvasLinkComponent,
    link_id: int,
) -> tuple[tuple[float, float], ...]:
    raw_table = component.properties.get("curve_points_q_head", "").strip()
    if not raw_table:
        raise ValueError(
            f"Pump in connection #{link_id} is missing its Q-Head curve table."
        )

    points: list[tuple[float, float]] = []
    for line_number, raw_line in enumerate(raw_table.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        tokens = [token for token in line.replace(",", " ").split() if token]
        if len(tokens) != 2:
            raise ValueError(
                f"Pump in connection #{link_id} has invalid Q-Head pair on line {line_number}: '{raw_line}'."
            )
        points.append((_parse_float(tokens[0]), _parse_float(tokens[1])))

    if not points:
        raise ValueError(
            f"Pump in connection #{link_id} must define at least one Q-Head pair."
        )

    return tuple(points)
