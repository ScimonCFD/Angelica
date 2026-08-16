from __future__ import annotations

from angelica.closures import (
    ColebrookPipeCorrelation,
    HeatSourceModel,
    LaminarPipeCorrelation,
    MinorLossModel,
    PumpCurveModel,
)
from angelica.closures.pressure_drop import PressureDropCorrelation
from angelica.core.case import NetworkCase
from angelica.core.components import Fitting, Pipe, Pump
from angelica.core.network import build_network_state
from angelica.core.results import ComponentFlowResult, IterationMetrics, SolveResult
from angelica.core.settings import SolverSettings
from angelica.core.state import FittingState, HeatSourceState, PipeState, PumpState
from angelica.numerics import assemble_pressure_system, max_abs_value, solve_linear_system
from .base import BaseSolver


class SteadyIsothermalIncompressibleSolver(BaseSolver):
    def __init__(
        self,
        settings: SolverSettings | None = None,
        laminar_pipe_correlation: LaminarPipeCorrelation | None = None,
        turbulent_pipe_correlation: PressureDropCorrelation | None = None,
        fitting_correlation: MinorLossModel | None = None,
        pump_correlation: PumpCurveModel | None = None,
        heat_source_correlation: HeatSourceModel | None = None,
    ):
        self.settings = settings or SolverSettings()
        self.laminar_pipe_correlation = laminar_pipe_correlation or LaminarPipeCorrelation()
        self.turbulent_pipe_correlation = turbulent_pipe_correlation or ColebrookPipeCorrelation()
        self.fitting_correlation = fitting_correlation or MinorLossModel()
        self.pump_correlation = pump_correlation or PumpCurveModel()
        self.heat_source_correlation = heat_source_correlation or HeatSourceModel()

    def solve(self, case: NetworkCase, progress_callback=None) -> SolveResult:
        network_state = build_network_state(case)
        self._initialise_pressure_field(network_state, case)

        laminar_history, laminar_metrics, laminar_converged = self._solve_laminar(
            network_state,
            case.fluid_model,
            progress_callback=progress_callback,
        )
        turbulent_history, turbulent_metrics, converged = self._solve_turbulent(
            network_state,
            case.fluid_model,
            progress_callback=progress_callback,
        )
        if self.settings.turbulent_iterations == 0:
            converged = laminar_converged

        node_pressures = {
            node_id: float(network_state.nodes[node_id].pressure_pa)
            for node_id in sorted(network_state.nodes)
        }
        component_flows = []
        for link in network_state.components:
            density = case.fluid_model.density_for_link(link)
            component_flows.append(
                ComponentFlowResult(
                    label=self._component_label(link),
                    mass_flow_kg_per_s=float(link.mass_flow_kg_per_s),
                    volumetric_flow_m3_per_h=float(3600.0 * link.mass_flow_kg_per_s / density),
                )
            )

        return SolveResult(
            case_name=case.name,
            converged=converged,
            node_pressures_pa=node_pressures,
            component_flows=component_flows,
            laminar_history=laminar_history,
            laminar_metrics=laminar_metrics,
            turbulent_history=turbulent_history,
            turbulent_metrics=turbulent_metrics,
            outer_turbulent_final_metrics=(turbulent_metrics[-1],) if turbulent_metrics else (),
        )

    def _component_label(self, link_state) -> str:
        component = link_state.component
        component_name = type(component).__name__
        if component.component_id:
            return f"{component_name}:{component.component_id}"
        return f"{component_name} {component.start_node}->{component.end_node}"

    def _initialise_pressure_field(self, network_state, case) -> None:
        boundaries = case.pressure_inlets + case.pressure_outlets
        inlet_pressures = [boundary.pressure_pa for boundary in boundaries if network_state.nodes[boundary.node_id].is_inlet]
        outlet_pressures = [boundary.pressure_pa for boundary in boundaries if not network_state.nodes[boundary.node_id].is_inlet]

        if inlet_pressures and outlet_pressures:
            pressure_max = max(inlet_pressures + outlet_pressures)
            pressure_min = min(inlet_pressures + outlet_pressures)
            delta_p = (pressure_max - pressure_min) / (2.0 * len(network_state.nodes))
            counter = 1
            for node in network_state.nodes.values():
                if node.pressure_pa is None:
                    node.pressure_pa = pressure_min + counter * delta_p
                counter += 1
            return

        if inlet_pressures:
            pressure_max = max(inlet_pressures)
            pressure_min = pressure_max / len(network_state.nodes)
            delta_p = (pressure_max - pressure_min) / (2.0 * len(network_state.nodes))
            counter = 0
            for node in network_state.nodes.values():
                if node.pressure_pa is None:
                    node.pressure_pa = pressure_max - counter * delta_p
                counter += 1
            return

        if outlet_pressures:
            pressure_min = min(outlet_pressures)
            delta_p = pressure_min / (2.0 * len(network_state.nodes))
            counter = 1
            for node in network_state.nodes.values():
                if node.pressure_pa is None:
                    node.pressure_pa = pressure_min + counter * delta_p
                counter += 1
            return

        raise ValueError("At least one pressure boundary is required.")

    def _solve_laminar(
        self,
        network_state,
        fluid_model,
        progress_callback=None,
    ) -> tuple[list[float], list[IterationMetrics], bool]:
        history = []
        metrics_history = []
        converged = False
        for iteration_index in range(self._get_laminar_iteration_count(network_state)):
            self._update_velocities(network_state, fluid_model, laminar=True)
            self._update_mass_flows(network_state, fluid_model)
            couplings = self._compute_couplings(network_state, fluid_model, laminar=True)
            matrix, vector = assemble_pressure_system(network_state, couplings)
            correction = self._solve_pressure_correction(matrix, vector)
            correction_abs, correction_mean_abs, correction_rel = self._apply_pressure_correction(
                network_state,
                correction,
            )
            history.append(correction_abs)
            self._update_velocities(network_state, fluid_model, laminar=True)
            self._update_mass_flows(network_state, fluid_model)
            max_mass_imbalance_rel = self._compute_max_nodal_mass_imbalance(network_state)
            metrics_history.append(
                self._build_iteration_metrics(
                    correction_abs,
                    correction_mean_abs,
                    correction_rel,
                    max_mass_imbalance_rel,
                )
            )
            if progress_callback is not None:
                progress_callback(
                    "laminar",
                    iteration_index + 1,
                    metrics_history[-1],
                )
            if (
                (
                    correction_abs <= self.settings.pressure_correction_abs_tolerance_pa
                    or correction_rel <= self.settings.pressure_correction_rel_tolerance
                )
                and max_mass_imbalance_rel <= self.settings.nodal_mass_imbalance_rel_tolerance
            ):
                converged = True
                break
        return history, metrics_history, converged

    def _get_laminar_iteration_count(self, network_state) -> int:
        if self.settings.laminar_iterations is not None:
            return self.settings.laminar_iterations

        has_special_components = any(
            isinstance(link_state, (FittingState, PumpState, HeatSourceState))
            for link_state in network_state.components
        )
        if has_special_components:
            return self.settings.laminar_iterations_with_fittings
        return self.settings.laminar_iterations_without_fittings

    def _solve_turbulent(
        self,
        network_state,
        fluid_model,
        progress_callback=None,
    ) -> tuple[list[float], list[IterationMetrics], bool]:
        history = []
        metrics_history = []
        converged = False
        for iteration_index in range(self.settings.turbulent_iterations):
            self._update_velocities(network_state, fluid_model, laminar=False)
            self._update_mass_flows(network_state, fluid_model)
            couplings = self._compute_couplings(network_state, fluid_model, laminar=False)
            matrix, vector = assemble_pressure_system(network_state, couplings)
            correction = self._solve_pressure_correction(matrix, vector)
            correction_abs, correction_mean_abs, correction_rel = self._apply_pressure_correction(
                network_state,
                correction,
            )
            history.append(correction_abs)
            self._update_velocities(network_state, fluid_model, laminar=False)
            self._update_mass_flows(network_state, fluid_model)
            max_mass_imbalance_rel = self._compute_max_nodal_mass_imbalance(network_state)
            metrics_history.append(
                self._build_iteration_metrics(
                    correction_abs,
                    correction_mean_abs,
                    correction_rel,
                    max_mass_imbalance_rel,
                )
            )
            if progress_callback is not None:
                progress_callback(
                    "turbulent",
                    iteration_index + 1,
                    metrics_history[-1],
                )
            if (
                (
                    correction_abs <= self.settings.pressure_correction_abs_tolerance_pa
                    or correction_rel <= self.settings.pressure_correction_rel_tolerance
                )
                and max_mass_imbalance_rel <= self.settings.nodal_mass_imbalance_rel_tolerance
            ):
                converged = True
                break
        return history, metrics_history, converged

    def _solve_pressure_correction(self, matrix, vector):
        if float(self.settings.pressure_relaxation) <= 0.0:
            raise ValueError("pressure_relaxation must be greater than zero.")
        return solve_linear_system(matrix, vector)

    def _build_iteration_metrics(
        self,
        correction_abs: float,
        correction_mean_abs: float,
        correction_rel: float,
        max_mass_imbalance_rel: float,
    ) -> IterationMetrics:
        return IterationMetrics(
            pressure_correction_abs_pa=correction_abs,
            pressure_correction_mean_abs_pa=correction_mean_abs,
            pressure_correction_rel=correction_rel,
            max_nodal_mass_imbalance_rel=max_mass_imbalance_rel,
        )

    def _update_velocities(self, network_state, fluid_model, laminar: bool) -> None:
        for link_state in network_state.components:
            density = fluid_model.density_for_link(link_state)
            viscosity = fluid_model.viscosity_for_link(link_state)
            delta_p = link_state.start_node.pressure_pa - link_state.end_node.pressure_pa
            correlation = self._get_pressure_drop_correlation(link_state, laminar)

            if isinstance(link_state, PipeState):
                self._update_reynolds(link_state, density, viscosity)
                if laminar:
                    link_state.velocity_m_per_s = correlation.calculate_velocity(
                        link_state,
                        delta_p,
                        density,
                        viscosity,
                    )
                else:
                    link_state.velocity_m_per_s = correlation.calculate_velocity(
                        link_state,
                        delta_p,
                        density,
                        viscosity,
                        self.settings.colebrook_residual_tolerance,
                        self.settings.colebrook_friction_strategy,
                        self.settings.friction_factor_method,
                        self.settings.friction_factor_max_iterations,
                        self.settings.velocity_loop_method,
                        self.settings.velocity_loop_max_iterations,
                        self.settings.velocity_loop_tolerance,
                    )
            elif isinstance(link_state, FittingState):
                link_state.velocity_m_per_s = correlation.calculate_velocity(
                    link_state,
                    delta_p,
                    density,
                    viscosity,
                )
                self._update_reynolds(link_state, density, viscosity)
            elif isinstance(link_state, PumpState):
                link_state.velocity_m_per_s = correlation.calculate_velocity(
                    link_state,
                    delta_p,
                    density,
                    viscosity,
                )
                self._update_reynolds(link_state, density, viscosity)
            elif isinstance(link_state, HeatSourceState):
                link_state.velocity_m_per_s = correlation.calculate_velocity(
                    link_state,
                    delta_p,
                    density,
                    viscosity,
                )
                self._update_reynolds(link_state, density, viscosity)
            else:
                raise TypeError(f"Unsupported state type: {type(link_state).__name__}")

    def _update_reynolds(self, link_state, density: float, viscosity: float) -> None:
        link_state.reynolds = density * abs(link_state.velocity_m_per_s) * link_state.diameter_m / viscosity

    def _update_mass_flows(self, network_state, fluid_model) -> None:
        for link_state in network_state.components:
            density = fluid_model.density_for_link(link_state)
            link_state.mass_flow_kg_per_s = density * link_state.velocity_m_per_s * link_state.area_m2

    def _compute_couplings(self, network_state, fluid_model, laminar: bool) -> list[float]:
        couplings = []
        for link_state in network_state.components:
            density = fluid_model.density_for_link(link_state)
            viscosity = fluid_model.viscosity_for_link(link_state)
            correlation = self._get_pressure_drop_correlation(link_state, laminar)
            coupling = correlation.calculate_coupling(link_state, density, viscosity)
            couplings.append(coupling)
        return couplings

    def _get_pressure_drop_correlation(self, link_state, laminar: bool):
        if isinstance(link_state, PipeState):
            if laminar:
                return self.laminar_pipe_correlation
            return self.turbulent_pipe_correlation
        if isinstance(link_state, FittingState):
            return self.fitting_correlation
        if isinstance(link_state, PumpState):
            return self.pump_correlation
        if isinstance(link_state, HeatSourceState):
            return self.heat_source_correlation
        raise TypeError(f"Unsupported state type: {type(link_state).__name__}")

    def _apply_pressure_correction(self, network_state, correction) -> tuple[float, float, float]:
        scaled_correction = []
        relative_corrections = []
        alpha = float(self.settings.pressure_relaxation)
        for idx, node_id in enumerate(sorted(network_state.nodes)):
            correction_value = float(correction[idx])
            delta_p = alpha * correction_value
            old_pressure = float(network_state.nodes[node_id].pressure_pa)
            new_pressure = old_pressure + delta_p
            network_state.nodes[node_id].pressure_pa = new_pressure
            scaled_correction.append(delta_p)
            denominator = max(abs(new_pressure), 1.0)
            relative_corrections.append(abs(delta_p) / denominator)
        abs_corrections = [abs(value) for value in scaled_correction]
        mean_abs_correction = sum(abs_corrections) / len(abs_corrections)
        return max_abs_value(scaled_correction), mean_abs_correction, max(relative_corrections)

    def _compute_max_nodal_mass_imbalance(self, network_state) -> float:
        """Return max relative mass imbalance across junction nodes: |in-out|/(in+out)."""
        imbalance_by_node = {node_id: 0.0 for node_id in network_state.nodes}
        throughput_by_node = {node_id: 0.0 for node_id in network_state.nodes}

        for link_state in network_state.components:
            mass_flow = float(link_state.mass_flow_kg_per_s)
            abs_flow = abs(mass_flow)
            imbalance_by_node[link_state.start_node.node_id] -= mass_flow
            imbalance_by_node[link_state.end_node.node_id] += mass_flow
            throughput_by_node[link_state.start_node.node_id] += abs_flow
            throughput_by_node[link_state.end_node.node_id] += abs_flow

        max_rel = 0.0
        for node_id, node_state in network_state.nodes.items():
            if node_state.is_boundary:
                continue
            imbalance = abs(imbalance_by_node[node_id])
            throughput = throughput_by_node[node_id]
            max_rel = max(max_rel, imbalance / max(throughput, 1e-10))

        return max_rel
