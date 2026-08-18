from __future__ import annotations

import json
import math
from pathlib import Path

from angelica.core.case import FlowBoundary, InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import FITTING_PRESET_LIBRARY, Fitting, HeatSource, Pipe, Pump
from angelica.core.results import ComponentFlowResult, IterationMetrics
from angelica.core.settings import SolverSettings
from angelica.closures import ColebrookPipeCorrelation, HazenWilliamsPipeCorrelation
from angelica.closures.convection_scheme import HybridScheme, PowerLawScheme, UpwindScheme
from angelica.properties.black_oil import BlackOilFluid
from angelica.properties.compressible_fluid import CompressibleFluid
from angelica.properties.eos import IdealGasEOS, PengRobinsonEOS
from angelica.properties.single_component import SingleComponentFluid
from angelica.properties.thermal_fluid import ThermalFluid
from angelica.solvers import (
    BlackOilSolverSettings,
    CompressibleSolverSettings,
    NonIsothermalSolverSettings,
    SteadyBlackOilSolver,
    SteadyCompressibleSolver,
    SteadyIsothermalIncompressibleSolver,
    SteadyNonIsothermalIncompressibleSolver,
)

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
    scene.case_name = str(data.get("case_name", ""))
    scene.physics_mode = str(data.get("physics_mode", "isothermal"))
    scene.solver_settings = dict(data.get("solver_settings", {}))
    scene.initial_node_pressures_pa = {
        int(node_id): float(value)
        for node_id, value in data.get("initial_node_pressures_pa", {}).items()
    }
    return scene


def load_scene_from_file(path: str | Path) -> CanvasScene:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return scene_from_dict(data)


def scene_to_dict(scene: CanvasScene) -> dict:
    return {
        "nodes": [
            {
                "node_id": node.node_id,
                "node_type": node.node_type,
                "x": node.x,
                "y": node.y,
                "properties": dict(node.properties),
            }
            for node in scene.nodes
        ],
        "links": [
            {
                "link_id": link.link_id,
                "start_node_id": link.start_node_id,
                "end_node_id": link.end_node_id,
                "components": [
                    {
                        "component_id": comp.component_id,
                        "component_type": comp.component_type,
                        "properties": dict(comp.properties),
                    }
                    for comp in link.components
                ],
            }
            for link in scene.links
        ],
        "case_name": scene.case_name,
        "physics_mode": scene.physics_mode,
        "material": dict(scene.material),
        "pressure_drop_model": dict(scene.pressure_drop_model),
        "solver_settings": dict(scene.solver_settings),
        "initial_node_pressures_pa": {
            str(k): v for k, v in scene.initial_node_pressures_pa.items()
        },
    }


def _metrics_to_dict(m: IterationMetrics) -> dict:
    return {
        "pressure_correction_abs_pa": m.pressure_correction_abs_pa,
        "pressure_correction_mean_abs_pa": m.pressure_correction_mean_abs_pa,
        "pressure_correction_rel": m.pressure_correction_rel,
        "max_nodal_mass_imbalance_rel": m.max_nodal_mass_imbalance_rel,
    }


def _metrics_from_dict(d: dict) -> IterationMetrics:
    return IterationMetrics(
        pressure_correction_abs_pa=float(d["pressure_correction_abs_pa"]),
        pressure_correction_mean_abs_pa=float(d["pressure_correction_mean_abs_pa"]),
        pressure_correction_rel=float(d["pressure_correction_rel"]),
        max_nodal_mass_imbalance_rel=float(d["max_nodal_mass_imbalance_rel"]),
    )


def _component_flow_to_dict(cf: ComponentFlowResult) -> dict:
    d: dict = {
        "label": cf.label,
        "mass_flow_kg_per_s": cf.mass_flow_kg_per_s,
        "volumetric_flow_m3_per_h": cf.volumetric_flow_m3_per_h,
    }
    if cf.temperature_in_c is not None:
        d["temperature_in_c"] = cf.temperature_in_c
    if cf.temperature_out_c is not None:
        d["temperature_out_c"] = cf.temperature_out_c
    return d


def _component_flow_from_dict(d: dict) -> ComponentFlowResult:
    return ComponentFlowResult(
        label=str(d["label"]),
        mass_flow_kg_per_s=float(d["mass_flow_kg_per_s"]),
        volumetric_flow_m3_per_h=float(d["volumetric_flow_m3_per_h"]),
        temperature_in_c=float(d["temperature_in_c"]) if "temperature_in_c" in d else None,
        temperature_out_c=float(d["temperature_out_c"]) if "temperature_out_c" in d else None,
    )


def save_scene_to_file(
    scene: CanvasScene,
    path: str | Path,
    boundary_results: dict[int, dict[str, float]] | None = None,
    convergence_history: dict[str, list[IterationMetrics]] | None = None,
    converged: bool = False,
    component_flows: list[ComponentFlowResult] | None = None,
) -> None:
    data = scene_to_dict(scene)
    if boundary_results:
        data["cached_results"] = {
            "converged": converged,
            "boundary_results": {str(k): v for k, v in boundary_results.items()},
            "convergence_history": {
                stage: [_metrics_to_dict(m) for m in metrics_list]
                for stage, metrics_list in (convergence_history or {}).items()
            },
            "component_flows": [_component_flow_to_dict(cf) for cf in (component_flows or [])],
        }
    Path(path).write_text(json.dumps(data, indent=2), encoding="utf-8")


def load_scene_and_results_from_file(
    path: str | Path,
) -> tuple[CanvasScene, dict | None]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    scene = scene_from_dict(data)
    raw = data.get("cached_results")
    if not raw:
        return scene, None
    cached: dict = {
        "converged": bool(raw.get("converged", False)),
        "boundary_results": {
            int(k): dict(v) for k, v in raw.get("boundary_results", {}).items()
        },
        "convergence_history": {
            stage: [_metrics_from_dict(m) for m in metrics_list]
            for stage, metrics_list in raw.get("convergence_history", {}).items()
        },
        "component_flows": [
            _component_flow_from_dict(cf) for cf in raw.get("component_flows", [])
        ],
    }
    return scene, cached


def build_network_case_from_scene(scene: CanvasScene) -> NetworkCase:
    if not scene.nodes:
        raise ValueError("The scene is empty. Add nodes before running the simulation.")
    if not scene.links:
        raise ValueError("The scene has no links. Add at least one connection before running.")

    is_compressible    = scene.physics_mode == "compressible"
    is_non_isothermal  = scene.physics_mode == "non_isothermal"
    is_black_oil       = scene.physics_mode == "black_oil"
    inlet_fluid_bcs: list[InletFluidBC] = []

    if not is_black_oil and not scene.material:
        raise ValueError("No material is defined. Use Material → Define Material before running.")

    if is_compressible:
        if not scene.material.get("molecular_weight_kg_per_mol", "").strip():
            raise ValueError(
                "Compressible mode requires Molecular Weight (M) in the material. "
                "Open Material → Define Material."
            )
        if not scene.material.get("viscosity_pa_s", "").strip():
            raise ValueError("The material is missing viscosity_pa_s.")
    elif not is_black_oil:
        if not scene.material.get("density_kg_per_m3", "").strip():
            raise ValueError("The material is missing density_kg_per_m3.")
        if not scene.material.get("viscosity_pa_s", "").strip():
            raise ValueError("The material is missing viscosity_pa_s.")

    if not (is_non_isothermal or is_compressible or is_black_oil):
        heat_source_link_ids = [
            link.link_id
            for link in scene.links
            for comp in link.components
            if comp.component_type == "heat_source"
        ]
        if heat_source_link_ids:
            ids = ", #".join(str(i) for i in heat_source_link_ids)
            raise ValueError(
                f"Connection(s) #{ids} contain a Heat Source, which requires "
                "Non-isothermal or Compressible physics mode. Change the physics mode "
                "in Settings or remove the Heat Source component."
            )

    if is_non_isothermal:
        cp_text = scene.material.get("specific_heat_j_per_kg_k", "").strip()
        k_text = scene.material.get("thermal_conductivity_w_per_m_k", "").strip()
        if not cp_text:
            raise ValueError(
                "Non-isothermal mode requires Specific Heat (cp) in the material. "
                "Open Material → Define Material."
            )
        if not k_text:
            raise ValueError(
                "Non-isothermal mode requires Thermal Conductivity (k) in the material. "
                "Open Material → Define Material."
            )

    for node in scene.nodes:
        if node.node_type != "junction":
            continue
        connection_count = sum(
            1 for link in scene.links
            if link.start_node_id == node.node_id or link.end_node_id == node.node_id
        )
        if connection_count == 0:
            raise ValueError(
                f"Junction #{node.node_id} is not connected to any pipe."
            )
        if connection_count == 1:
            raise ValueError(
                f"Junction #{node.node_id} has only one connection — "
                "it acts as a dead end. Add another connection or change this node to a source or sink."
            )

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
                heat_transfer = _optional_float(
                    component, "heat_transfer_coefficient_w_per_m2k", default=0.0
                ) if (is_non_isothermal or is_compressible or is_black_oil) else 0.0
                ambient_temp = _optional_float(
                    component, "ambient_temperature_c", default=20.0
                ) if (is_non_isothermal or is_compressible or is_black_oil) else 20.0
                num_segs = max(1, int(_optional_float(component, "num_segments", default=1.0)))
                seg_length = length / num_segs
                seg_height = height_change / num_segs
                seg_nodes = [current_start]
                for _ in range(num_segs - 1):
                    seg_nodes.append(next_internal_node_id)
                    next_internal_node_id += 1
                seg_nodes.append(current_end)
                for seg_idx in range(num_segs):
                    components.append(
                        Pipe(
                            start_node=seg_nodes[seg_idx],
                            end_node=seg_nodes[seg_idx + 1],
                            diameter_m=diameter,
                            length_m=seg_length,
                            absolute_roughness_m=roughness,
                            hazen_williams_c=hazen_williams_c,
                            height_change_m=seg_height,
                            heat_transfer_coefficient_w_per_m2k=heat_transfer,
                            ambient_temperature_c=ambient_temp,
                            n_thermal_segments=1,
                            component_id=f"link_{link.link_id}_pipe_{component.component_id}_seg{seg_idx}",
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
            elif component.component_type == "heat_source":
                diameter = _required_float(component, "diameter_m", link.link_id)
                power_w = _required_float(component, "power_w", link.link_id)
                mode = component.properties.get("pressure_drop_mode", "rated").strip() or "rated"
                dp = _optional_float(component, "pressure_drop_pa", default=0.0)
                mdot_rated = _optional_float(component, "rated_mass_flow_kg_per_s", default=1.0)
                n_segs = max(2, int(_optional_float(
                    component, "n_thermal_segments", default=10.0
                ))) if (is_non_isothermal or is_compressible or is_black_oil) else 2
                components.append(
                    HeatSource(
                        start_node=current_start,
                        end_node=current_end,
                        diameter_m=diameter,
                        power_w=power_w,
                        pressure_drop_mode=mode,
                        pressure_drop_pa=dp,
                        rated_mass_flow_kg_per_s=mdot_rated,
                        n_thermal_segments=n_segs,
                        component_id=f"link_{link.link_id}_heat_source_{component.component_id}",
                    )
                )
            else:
                raise ValueError(
                    f"Unsupported component type '{component.component_type}' in connection #{link.link_id}."
                )

            current_start = current_end

    visible_node_ids = {node.node_id for node in scene.nodes}
    all_node_ids = visible_node_ids.union(range(max(visible_node_ids) + 1, next_internal_node_id))

    if is_compressible:
        m_text = scene.material.get("molecular_weight_kg_per_mol", "0.028964").strip()
        cp_text = scene.material.get("specific_heat_j_per_kg_k", "1000.0").strip()
        k_text = scene.material.get("thermal_conductivity_w_per_m_k", "0.025").strip()
        if scene.material.get("definition_mode") == "gas_pr":
            eos = PengRobinsonEOS(
                molecular_weight_kg_per_mol=float(m_text),
                critical_temperature_k=float(scene.material["critical_temperature_k"]),
                critical_pressure_pa=float(scene.material["critical_pressure_pa"]),
                acentric_factor=float(scene.material["acentric_factor"]),
            )
        else:
            eos = IdealGasEOS(float(m_text))
        fluid_model = CompressibleFluid.from_constants(
            eos=eos,
            viscosity_pa_s=float(scene.material["viscosity_pa_s"]),
            specific_heat_j_per_kg_k=float(cp_text) if cp_text else 1000.0,
            thermal_conductivity_w_per_m_k=float(k_text) if k_text else 0.025,
        )
        thermal_inlets = tuple(
            tb
            for tb in (
                _build_thermal_boundary(node)
                for node in scene.nodes
                if node.node_type in ("source", "sink")
            )
            if tb is not None
        )
        if not any(tb.bc_type == "fixed_temperature" for tb in thermal_inlets):
            raise ValueError(
                "Compressible mode requires at least one boundary node with a fixed "
                "temperature. Open a source or sink node and set its thermal boundary "
                "condition to 'Fixed temperature'."
            )
    elif is_non_isothermal:
        fluid_model = ThermalFluid.from_constants(
            density_kg_per_m3=float(scene.material["density_kg_per_m3"]),
            viscosity_pa_s=float(scene.material["viscosity_pa_s"]),
            specific_heat_j_per_kg_k=float(scene.material["specific_heat_j_per_kg_k"]),
            thermal_conductivity_w_per_m_k=float(scene.material["thermal_conductivity_w_per_m_k"]),
        )
        thermal_inlets = tuple(
            tb
            for tb in (
                _build_thermal_boundary(node)
                for node in scene.nodes
                if node.node_type in ("source", "sink")
            )
            if tb is not None
        )
        if not any(tb.bc_type == "fixed_temperature" for tb in thermal_inlets):
            raise ValueError(
                "Non-isothermal mode requires at least one boundary node with a fixed "
                "temperature. Open a source or sink node and set its thermal boundary "
                "condition to 'Fixed temperature'."
            )
    elif is_black_oil:
        _BO_FIELDS = ("api_gravity", "gas_gravity", "gor_sc_m3_per_m3", "wor_sc_m3_per_m3")
        inlet_fluid_bcs: list[InletFluidBC] = []
        for node in scene.nodes:
            if node.node_type != "source":
                continue
            missing = [f for f in _BO_FIELDS if not node.properties.get(f, "").strip()]
            if missing:
                raise ValueError(
                    f"Source node #{node.node_id} is missing black-oil composition fields: "
                    + ", ".join(missing)
                    + ". Open the node properties and fill in all four fields."
                )
            inlet_fluid_bcs.append(InletFluidBC(
                node_id          = node.node_id,
                api_gravity      = float(node.properties["api_gravity"]),
                gas_gravity      = float(node.properties["gas_gravity"]),
                gor_sc_m3_per_m3 = float(node.properties["gor_sc_m3_per_m3"]),
                wor_sc_m3_per_m3 = float(node.properties["wor_sc_m3_per_m3"]),
            ))
        if not inlet_fluid_bcs:
            raise ValueError(
                "Black-oil mode requires at least one source node with a defined "
                "fluid composition (API, gas gravity, GOR, WOR)."
            )
        first = inlet_fluid_bcs[0]
        fluid_model = BlackOilFluid(
            api_gravity      = first.api_gravity,
            gas_gravity      = first.gas_gravity,
            gor_sc_m3_per_m3 = first.gor_sc_m3_per_m3,
            wor_sc_m3_per_m3 = first.wor_sc_m3_per_m3,
        )
        thermal_inlets = tuple(
            tb
            for tb in (
                _build_thermal_boundary(node)
                for node in scene.nodes
                if node.node_type in ("source", "sink")
            )
            if tb is not None
        )
        if not any(tb.bc_type == "fixed_temperature" for tb in thermal_inlets):
            raise ValueError(
                "Black-oil mode requires at least one boundary node with a fixed "
                "temperature. Open a source or sink node and set its thermal boundary "
                "condition to 'Fixed temperature'."
            )
    else:
        fluid_model = SingleComponentFluid(
            density_kg_per_m3=float(scene.material["density_kg_per_m3"]),
            viscosity_pa_s=float(scene.material["viscosity_pa_s"]),
        )
        thermal_inlets = ()

    return NetworkCase(
        name=scene.case_name or "GUI scene",
        fluid_model=fluid_model,
        pressure_inlets=tuple(pressure_inlets),
        pressure_outlets=tuple(pressure_outlets),
        flow_inlets=tuple(flow_inlets),
        flow_outlets=tuple(flow_outlets),
        components=tuple(components),
        node_ids=tuple(sorted(all_node_ids)),
        initial_node_pressures_pa=dict(scene.initial_node_pressures_pa),
        thermal_inlets=thermal_inlets,
        inlet_fluid_bcs=tuple(inlet_fluid_bcs) if is_black_oil else (),
    )


def build_solver_from_scene(scene: CanvasScene):
    pressure_drop_model_key = scene.pressure_drop_model.get("library_key", "")
    if pressure_drop_model_key == "colebrook_white":
        turbulent_pipe_correlation = ColebrookPipeCorrelation()
    elif pressure_drop_model_key == "hazen_williams":
        turbulent_pipe_correlation = HazenWilliamsPipeCorrelation()
    else:
        raise ValueError(
            f"Unsupported pipe pressure-drop model '{pressure_drop_model_key}'."
        )

    _NI_KEYS = {
        "max_temperature_iterations",
        "temperature_tolerance_k",
        "temperature_relaxation",
        "convection_scheme",
    }
    _COMP_KEYS = {
        "max_density_iterations", "density_rel_tolerance",
        "temperature_tolerance_k", "temperature_relaxation", "convection_scheme",
    }
    _OUTER_KEYS = _NI_KEYS | _COMP_KEYS
    ni_raw = {k: scene.solver_settings[k] for k in _NI_KEYS if k in scene.solver_settings}
    comp_raw = {k: scene.solver_settings[k] for k in _COMP_KEYS if k in scene.solver_settings}
    hyd_raw = {k: v for k, v in scene.solver_settings.items() if k not in _OUTER_KEYS}

    settings = SolverSettings(**hyd_raw) if hyd_raw else SolverSettings()

    if scene.physics_mode == "black_oil":
        return SteadyBlackOilSolver(
            hydraulic_settings=settings if hyd_raw else None,
            turbulent_pipe_correlation=turbulent_pipe_correlation,
        )

    if scene.physics_mode == "compressible":
        _CONVECTION_SCHEMES = {
            "upwind": UpwindScheme,
            "hybrid": HybridScheme,
            "power_law": PowerLawScheme,
        }
        comp_kwargs: dict = {}
        if "max_density_iterations" in comp_raw:
            comp_kwargs["max_density_iterations"] = int(comp_raw["max_density_iterations"])
        if "density_rel_tolerance" in comp_raw:
            comp_kwargs["density_rel_tolerance"] = float(comp_raw["density_rel_tolerance"])
        if "temperature_tolerance_k" in comp_raw:
            comp_kwargs["temperature_tolerance_k"] = float(comp_raw["temperature_tolerance_k"])
        if "temperature_relaxation" in comp_raw:
            comp_kwargs["temperature_relaxation"] = float(comp_raw["temperature_relaxation"])
        scheme_key = str(comp_raw.get("convection_scheme", "hybrid"))
        convection_scheme = _CONVECTION_SCHEMES.get(scheme_key, HybridScheme)()
        return SteadyCompressibleSolver(
            hydraulic_settings=settings if hyd_raw else None,
            compressible_settings=CompressibleSolverSettings(**comp_kwargs),
            turbulent_pipe_correlation=turbulent_pipe_correlation,
            convection_scheme=convection_scheme,
        )

    if scene.physics_mode == "non_isothermal":
        ni_kwargs: dict = {}
        if "max_temperature_iterations" in ni_raw:
            ni_kwargs["max_temperature_iterations"] = int(ni_raw["max_temperature_iterations"])
        if "temperature_tolerance_k" in ni_raw:
            ni_kwargs["temperature_tolerance_k"] = float(ni_raw["temperature_tolerance_k"])
        if "temperature_relaxation" in ni_raw:
            ni_kwargs["temperature_relaxation"] = float(ni_raw["temperature_relaxation"])

        _CONVECTION_SCHEMES = {
            "upwind": UpwindScheme,
            "hybrid": HybridScheme,
            "power_law": PowerLawScheme,
        }
        scheme_key = str(ni_raw.get("convection_scheme", "upwind"))
        convection_scheme = _CONVECTION_SCHEMES.get(scheme_key, UpwindScheme)()

        return SteadyNonIsothermalIncompressibleSolver(
            hydraulic_settings=settings,
            non_isothermal_settings=NonIsothermalSolverSettings(**ni_kwargs),
            turbulent_pipe_correlation=turbulent_pipe_correlation,
            convection_scheme=convection_scheme,
        )

    return SteadyIsothermalIncompressibleSolver(
        settings=settings,
        turbulent_pipe_correlation=turbulent_pipe_correlation,
    )


def _build_thermal_boundary(node) -> ThermalBoundary | None:
    """Return a ThermalBoundary for a source/sink node, or None if no thermal BC applies."""
    props = node.properties
    bc_type = props.get("thermal_bc_type", "")

    # Backward-compat: old files lack thermal_bc_type; infer from inlet_temperature_c
    if not bc_type:
        t_text = props.get("inlet_temperature_c", "").strip()
        if node.node_type == "source" and t_text:
            bc_type = "fixed_temperature"
        else:
            bc_type = "zero_gradient"

    if bc_type == "fixed_temperature":
        t_text = props.get("inlet_temperature_c", "").strip()
        if not t_text:
            return None  # no temperature set — skip (no thermal BC)
        return ThermalBoundary(
            node_id=node.node_id,
            temperature_c=float(t_text),
            bc_type="fixed_temperature",
            gradient_dc_per_m=0.0,
        )

    if bc_type == "fixed_gradient":
        g_text = props.get("thermal_gradient_dc_per_m", "0.0").strip()
        g = float(g_text) if g_text else 0.0
        return ThermalBoundary(
            node_id=node.node_id,
            temperature_c=0.0,
            bc_type="fixed_gradient",
            gradient_dc_per_m=g,
        )

    # zero_gradient (or anything unrecognised): Neumann with g=0
    return ThermalBoundary(
        node_id=node.node_id,
        temperature_c=0.0,
        bc_type="zero_gradient",
        gradient_dc_per_m=0.0,
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
