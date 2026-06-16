from __future__ import annotations

from dataclasses import dataclass

from angelica.closures.convection_scheme import ConvectionScheme, UpwindScheme
from angelica.core.network import build_network_state
from angelica.core.results import SolveResult
from angelica.core.settings import SolverSettings
from angelica.core.state import PipeState
from angelica.numerics.energy import solve_energy_system
from .base import BaseSolver
from .steady_isothermal_incompressible import SteadyIsothermalIncompressibleSolver


@dataclass
class NonIsothermalSolverSettings:
    """Settings for the outer temperature iteration loop."""

    max_temperature_iterations: int = 50
    temperature_tolerance_k: float = 0.01
    temperature_relaxation: float = 1.0


class SteadyNonIsothermalIncompressibleSolver(BaseSolver):
    """Steady, incompressible, non-isothermal pipe network solver.

    Outer loop:
        1. Solve hydraulics (SIMPLE) with current ρ(T) and μ(T).
        2. Solve energy equation for nodal temperatures.
        3. Update T on each pipe state.
        4. Repeat until max |ΔT| < tolerance.

    The inner hydraulic solver is the existing
    SteadyIsothermalIncompressibleSolver — it runs to full convergence on
    every outer iteration.
    """

    def __init__(
        self,
        hydraulic_settings: SolverSettings | None = None,
        non_isothermal_settings: NonIsothermalSolverSettings | None = None,
        convection_scheme: ConvectionScheme | None = None,
    ) -> None:
        self.hydraulic_settings = hydraulic_settings or SolverSettings()
        self.non_isothermal_settings = non_isothermal_settings or NonIsothermalSolverSettings()
        self.convection_scheme = convection_scheme or UpwindScheme()
        self._hydraulic_solver = SteadyIsothermalIncompressibleSolver(
            settings=self.hydraulic_settings
        )

    def solve(self, case, _progress_callback=None) -> SolveResult:
        network_state = build_network_state(case)

        # ── initialise thermal boundary conditions on nodes ────────────────
        for tb in case.thermal_inlets:
            if tb.node_id in network_state.nodes:
                network_state.nodes[tb.node_id].temperature_c = tb.temperature_c
                network_state.nodes[tb.node_id].is_thermal_inlet = True

        # ── initialise temperature field ───────────────────────────────────
        # Use the temperature of the inlet with the most prescribed T as seed.
        T_init = self._initial_temperature(case, network_state)
        for node in network_state.nodes.values():
            if node.temperature_c is None:
                node.temperature_c = T_init
        for ps in network_state.components:
            if isinstance(ps, PipeState):
                ps.temperature_c = T_init

        # ── outer temperature loop ─────────────────────────────────────────
        settings = self.non_isothermal_settings
        temperature_converged = False

        for _outer in range(settings.max_temperature_iterations):
            # 1. Solve hydraulics with current T
            self._hydraulic_solver._initialise_pressure_field(network_state, case)
            self._hydraulic_solver._solve_laminar(network_state, case.fluid_model)
            self._hydraulic_solver._solve_turbulent(network_state, case.fluid_model)

            # 2. Solve energy equation
            new_node_temps = solve_energy_system(
                network_state, case.fluid_model, self.convection_scheme
            )

            # 3. Check convergence and update temperatures
            max_delta_t = 0.0
            alpha = settings.temperature_relaxation
            for nid, T_new in new_node_temps.items():
                T_old = network_state.nodes[nid].temperature_c or T_init
                delta = alpha * (T_new - T_old)
                max_delta_t = max(max_delta_t, abs(delta))
                if not network_state.nodes[nid].is_thermal_inlet:
                    network_state.nodes[nid].temperature_c = T_old + delta

            # 4. Update pipe representative temperature
            for ps in network_state.components:
                if isinstance(ps, PipeState):
                    T_start = network_state.nodes[ps.start_node.node_id].temperature_c or T_init
                    T_end = network_state.nodes[ps.end_node.node_id].temperature_c or T_init
                    ps.temperature_c = 0.5 * (T_start + T_end)

            if max_delta_t < settings.temperature_tolerance_k:
                temperature_converged = True
                break

        # ── final hydraulic solve with converged temperatures ──────────────
        self._hydraulic_solver._initialise_pressure_field(network_state, case)
        lam_hist, lam_metrics, _ = self._hydraulic_solver._solve_laminar(
            network_state, case.fluid_model
        )
        turb_hist, turb_metrics, hydraulic_converged = self._hydraulic_solver._solve_turbulent(
            network_state, case.fluid_model
        )

        converged = hydraulic_converged and temperature_converged

        # ── build results ──────────────────────────────────────────────────
        node_pressures = {
            nid: float(network_state.nodes[nid].pressure_pa)
            for nid in sorted(network_state.nodes)
        }
        node_temperatures = {
            nid: float(network_state.nodes[nid].temperature_c or T_init)
            for nid in sorted(network_state.nodes)
        }
        from angelica.core.results import ComponentFlowResult
        component_flows = []
        for link in network_state.components:
            density = case.fluid_model.density_for_link(link)
            component_flows.append(
                ComponentFlowResult(
                    label=self._hydraulic_solver._component_label(link),
                    mass_flow_kg_per_s=float(link.mass_flow_kg_per_s),
                    volumetric_flow_m3_per_h=float(
                        3600.0 * link.mass_flow_kg_per_s / density
                    ),
                )
            )

        return SolveResult(
            case_name=case.name,
            converged=converged,
            node_pressures_pa=node_pressures,
            node_temperatures_c=node_temperatures,
            component_flows=component_flows,
            laminar_history=lam_hist,
            laminar_metrics=lam_metrics,
            turbulent_history=turb_hist,
            turbulent_metrics=turb_metrics,
        )

    @staticmethod
    def _initial_temperature(case, network_state) -> float:
        if case.thermal_inlets:
            return case.thermal_inlets[0].temperature_c
        # fallback: ambient of first pipe, or 20 °C
        for comp in case.components:
            if hasattr(comp, "ambient_temperature_c"):
                return comp.ambient_temperature_c
        return 20.0
