from __future__ import annotations

import math

from angelica.closures import (
    ColebrookPipeCorrelation,
    HeatSourceModel,
    LaminarPipeCorrelation,
    MinorLossModel,
    PumpCurveModel,
)
from angelica.closures.pressure_drop import PressureDropCorrelation
from angelica.core.case import NetworkCase
from angelica.core.network import build_network_state
from angelica.core.results import ComponentFlowResult, IterationMetrics, SolveResult
from angelica.core.settings import SolverSettings
from angelica.core.state import CompressorState, FittingState, HeatSourceState, PipeState, PumpState
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
        from angelica.closures.compressor import CompressorModel
        self.compressor_correlation = CompressorModel()
        self.pressure_floor_pa: float = 0.0     # gas / black-oil / compositional solvers set 1 kPa

    def solve(self, case: NetworkCase, progress_callback=None) -> SolveResult:
        from angelica.core.components import Compressor
        if any(isinstance(c, Compressor) for c in case.components):
            raise ValueError("A Compressor needs a gas, black-oil or compositional solver (its discharge "
                             "temperature and power come from the fluid's enthalpy); use a Pump for liquids.")
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
            node_id: float(network_state.nodes[node_id].pressure_pa)  # type: ignore[arg-type]
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

        _balance_source = turbulent_metrics if turbulent_metrics else laminar_metrics
        mass_balance_history = [
            m.global_mass_imbalance_rel * 100.0 for m in _balance_source
        ]

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
            global_balance=self._compute_global_balance(network_state),
            global_energy_balance=None,
            mass_balance_history=mass_balance_history,
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
        return self._run_iteration_loop(
            network_state,
            fluid_model,
            laminar=True,
            n_iterations=self._get_laminar_iteration_count(network_state),
            phase_label="laminar",
            progress_callback=progress_callback,
        )

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
        return self._run_iteration_loop(
            network_state,
            fluid_model,
            laminar=False,
            n_iterations=self.settings.turbulent_iterations,
            phase_label="turbulent",
            progress_callback=progress_callback,
        )

    def _run_iteration_loop(
        self,
        network_state,
        fluid_model,
        laminar: bool,
        n_iterations: int,
        phase_label: str,
        progress_callback=None,
    ) -> tuple[list[float], list[IterationMetrics], bool]:
        history: list[float] = []
        metrics_history: list[IterationMetrics] = []
        converged = False
        for iteration_index in range(n_iterations):
            self._update_velocities(network_state, fluid_model, laminar=laminar)
            self._update_mass_flows(network_state, fluid_model)
            couplings = self._compute_couplings(network_state, fluid_model, laminar=laminar)
            matrix, vector = assemble_pressure_system(network_state, couplings)
            correction = self._solve_pressure_correction(matrix, vector)
            correction_abs, correction_mean_abs, correction_rel = self._apply_pressure_correction(
                network_state, correction,
            )
            history.append(correction_abs)
            self._update_velocities(network_state, fluid_model, laminar=laminar)
            self._update_mass_flows(network_state, fluid_model)
            max_mass_imbalance_rel = self._compute_max_nodal_mass_imbalance(network_state)
            global_balance = self._compute_global_balance(network_state)
            metrics_history.append(
                self._build_iteration_metrics(
                    correction_abs,
                    correction_mean_abs,
                    correction_rel,
                    max_mass_imbalance_rel,
                    abs(global_balance.mass_inlet_kg_per_s - global_balance.mass_outlet_kg_per_s),
                    global_balance.mass_error_pct / 100.0,
                )
            )
            if progress_callback is not None:
                progress_callback(phase_label, iteration_index + 1, metrics_history[-1])
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
        global_mass_imbalance_kg_per_s: float,
        global_mass_imbalance_rel: float,
    ) -> IterationMetrics:
        return IterationMetrics(
            pressure_correction_abs_pa=correction_abs,
            pressure_correction_mean_abs_pa=correction_mean_abs,
            pressure_correction_rel=correction_rel,
            max_nodal_mass_imbalance_rel=max_mass_imbalance_rel,
            global_mass_imbalance_kg_per_s=global_mass_imbalance_kg_per_s,
            global_mass_imbalance_rel=global_mass_imbalance_rel,
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
                        self.settings.laminar_turbulent_transition_re,
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
        if isinstance(link_state, CompressorState):
            return self.compressor_correlation
        if isinstance(link_state, PumpState):
            return self.pump_correlation
        if isinstance(link_state, HeatSourceState):
            return self.heat_source_correlation
        raise TypeError(f"Unsupported state type: {type(link_state).__name__}")

    def _apply_pressure_correction(self, network_state, correction) -> tuple[float, float, float]:
        # Pressure floor: 0 Pa for liquids (gauge pressures are common, e.g. the
        # Poiseuille benchmark's outlet at 0); the gas, black-oil and
        # compositional solvers set 1 kPa absolute — at exactly 0 the equation
        # of state divides by P (a 30 MW compressor emptying its suction line
        # raised ZeroDivisionError, v1.30.0).
        # The floor never exceeds half the lowest boundary pressure, and
        # boundary nodes keep their own pressure (review of v1.30.0: a 0.5 kPa
        # outlet was lifted to the 1 kPa floor).
        node_ids = sorted(network_state.nodes)
        _P_MIN = float(getattr(self, "pressure_floor_pa", 0.0))
        if _P_MIN > 0.0:
            bc = [float(network_state.nodes[n].pressure_pa) for n in node_ids  # type: ignore[arg-type]
                  if network_state.nodes[n].is_pressure_boundary and network_state.nodes[n].pressure_pa is not None]
            if bc:
                _P_MIN = min(_P_MIN, 0.5 * max(min(bc), 0.0))
        # Pressure boundaries keep their pressure; every other node (flow
        # boundaries too: their pressure is free) stays above the floor.
        floors = [(-math.inf if network_state.nodes[n].is_pressure_boundary else _P_MIN) for n in node_ids]
        old_pressures = [float(network_state.nodes[nid].pressure_pa) for nid in node_ids]  # type: ignore[arg-type]
        raw_corrections = [float(correction[i]) for i in range(len(node_ids))]

        # Adaptive relaxation: halve alpha until all new pressures are above the floor.
        alpha = float(self.settings.pressure_relaxation)
        alpha_floor = alpha * (0.5 ** 20)  # ~20 halvings before giving up
        while alpha > alpha_floor:
            if all(p + alpha * dp >= fl for p, dp, fl in zip(old_pressures, raw_corrections, floors)):
                break
            alpha *= 0.5

        new_pressures = [max(p + alpha * dp, fl)
                         for p, dp, fl in zip(old_pressures, raw_corrections, floors)]
        for nid, new_p in zip(node_ids, new_pressures):
            network_state.nodes[nid].pressure_pa = new_p

        deltas = [new_p - old_p for new_p, old_p in zip(new_pressures, old_pressures)]
        abs_deltas = [abs(d) for d in deltas]
        relative_corrections = [abs(d) / max(abs(new_p), 1.0)
                                 for d, new_p in zip(deltas, new_pressures)]
        return max_abs_value(deltas), sum(abs_deltas) / len(abs_deltas), max(relative_corrections)

    def _compute_max_nodal_mass_imbalance(self, network_state) -> float:
        """Return max relative mass imbalance across junction nodes: |in-out|/(in+out)."""
        return max(nodal_mass_imbalances(network_state).values(), default=0.0)


def nodal_mass_imbalances(network_state) -> dict[int, float]:
    """Relative mass imbalance |in − out|/(in + out) of each interior node
    (the hydraulic convergence measure; see the comments below)."""
    imbalance_by_node = dict.fromkeys(network_state.nodes, 0.0)
    throughput_by_node = dict.fromkeys(network_state.nodes, 0.0)

    for link_state in network_state.components:
        mass_flow = float(link_state.mass_flow_kg_per_s)
        abs_flow = abs(mass_flow)
        imbalance_by_node[link_state.start_node.node_id] -= mass_flow
        imbalance_by_node[link_state.end_node.node_id] += mass_flow
        throughput_by_node[link_state.start_node.node_id] += abs_flow
        throughput_by_node[link_state.end_node.node_id] += abs_flow
    for node_id, node_state in network_state.nodes.items():
        src = getattr(node_state, "source_mass_flow_kg_per_s", 0.0)
        if src:
            imbalance_by_node[node_id] += src
            throughput_by_node[node_id] += abs(src)

    # Nodes with (almost) no flow — e.g. behind a closed check valve or on a
    # dead-end branch — carry round-off-level flows; dividing by their own
    # throughput would report ~100 % imbalance forever.  Measure them
    # against the typical throughput of the network instead.
    max_throughput = max(throughput_by_node.values(), default=0.0)
    if max_throughput < 1e-6:
        # The whole network is at rest (e.g. its only well is dead): the
        # remaining flows are velocity-tolerance noise, not imbalance.
        return {}
    throughput_floor = max(1e-10, 1e-3 * max_throughput)
    return {node_id: abs(imbalance_by_node[node_id]) / max(throughput_by_node[node_id], throughput_floor)
            for node_id, node_state in network_state.nodes.items() if not node_state.is_boundary}
