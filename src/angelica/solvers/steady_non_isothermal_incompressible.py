from __future__ import annotations

from dataclasses import dataclass

from angelica.closures.convection_scheme import ConvectionScheme, UpwindScheme
from angelica.closures.pressure_drop import PressureDropCorrelation
from angelica.core.case import NetworkCase
from angelica.core.network import build_network_state
from angelica.core.results import SolveResult
from angelica.core.settings import SolverSettings
from angelica.core.state import HeatSourceState, PipeState
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

    Applicability and assumptions
    ------------------------------
    - **Incompressible fluid**: density is not a function of pressure.
      Density and all other fluid properties (μ, cp, k) may vary with
      temperature but are evaluated at the local temperature only, not at the
      local pressure.  This is appropriate for liquids (water, thermal oils,
      crude oil) under moderate pressures.
    - **Not valid for compressible gases** where ρ = ρ(P, T), or for
      liquids at pressures high enough that pressure-induced property changes
      are significant.
    - **Steady-state only**: transient effects (thermal capacitance, hydraulic
      inertia) are not modelled.
    """

    def __init__(
        self,
        hydraulic_settings: SolverSettings | None = None,
        non_isothermal_settings: NonIsothermalSolverSettings | None = None,
        convection_scheme: ConvectionScheme | None = None,
        turbulent_pipe_correlation: PressureDropCorrelation | None = None,
    ) -> None:
        self.hydraulic_settings = hydraulic_settings or SolverSettings()
        self.non_isothermal_settings = non_isothermal_settings or NonIsothermalSolverSettings()
        self.convection_scheme = convection_scheme or UpwindScheme()
        self._hydraulic_solver = SteadyIsothermalIncompressibleSolver(
            settings=self.hydraulic_settings,
            turbulent_pipe_correlation=turbulent_pipe_correlation,
        )

    def solve(self, case: NetworkCase, progress_callback=None) -> SolveResult:
        from dataclasses import replace as _replace

        from angelica.numerics.energy import solve_energy_system
        from angelica.properties.thermal_fluid import ThermalFluid

        self._require_fixed_temperature(case)
        network_state = build_network_state(case)

        # ── initialise thermal boundary conditions on nodes ────────────────
        for tb in case.thermal_inlets:
            if tb.node_id not in network_state.nodes:
                continue
            node_st = network_state.nodes[tb.node_id]
            if tb.bc_type == "fixed_temperature":
                node_st.temperature_c = tb.temperature_c
                node_st.is_thermal_inlet = True
            else:
                # "zero_gradient" and "fixed_gradient" both set the Neumann
                # gradient; zero_gradient is the special case gradient = 0.
                node_st.thermal_gradient_dc_per_m = tb.gradient_dc_per_m

        # ── initialise temperature field ───────────────────────────────────
        # Use the temperature of the inlet with the most prescribed T as seed.
        T_init = self._initial_temperature(case)
        fluid_model = case.fluid_model
        if isinstance(fluid_model, ThermalFluid) and fluid_model.reference_temperature_c != T_init:
            fluid_model = _replace(fluid_model, reference_temperature_c=T_init)

        for node in network_state.nodes.values():
            if node.temperature_c is None:
                node.temperature_c = T_init
        for ps in network_state.components:
            if isinstance(ps, PipeState):
                ps.temperature_c = T_init

        # ── outer temperature loop ─────────────────────────────────────────
        settings = self.non_isothermal_settings
        temperature_converged = False
        temperature_history: list[float] = []
        mass_balance_history: list[float] = []
        energy_balance_history: list[float] = []
        all_lam_hist: list[float] = []
        all_lam_metrics = []
        all_turb_hist: list[float] = []
        all_turb_metrics = []
        outer_turb_final = []
        outer_boundaries: list[int] = []
        lam_hist: list[float] = []
        lam_metrics = []
        turb_hist: list[float] = []
        turb_metrics = []
        hydraulic_converged = False

        for _outer in range(settings.max_temperature_iterations):
            # 1. Solve hydraulics with current T
            self._hydraulic_solver._initialise_pressure_field(network_state, case)
            lam_hist, lam_metrics, _ = self._hydraulic_solver._solve_laminar(
                network_state,
                fluid_model,
            )
            turb_hist, turb_metrics, hydraulic_converged = self._hydraulic_solver._solve_turbulent(
                network_state,
                fluid_model,
            )
            all_lam_hist.extend(lam_hist)
            all_lam_metrics.extend(lam_metrics)
            all_turb_hist.extend(turb_hist)
            all_turb_metrics.extend(turb_metrics)
            outer_boundaries.append(len(all_turb_metrics))
            if turb_metrics:
                outer_turb_final.append(turb_metrics[-1])

            # 2. Solve energy equation
            new_node_temps, pipe_mean_temps = solve_energy_system(
                network_state, fluid_model, self.convection_scheme, T_ref=T_init
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

            temperature_history.append(max_delta_t)

            # 4. Update representative temperature for pipes and heat sources
            # Use internal FV node mean temperatures for accurate property evaluation.
            self._update_component_temperatures(network_state, T_init, pipe_mean_temps)

            _gb = self._compute_global_balance(network_state)
            mass_balance_history.append(_gb.mass_error_pct)
            _geb = self._compute_global_energy_balance(network_state, fluid_model)
            energy_balance_history.append(_geb.energy_error_pct)

            if max_delta_t < settings.temperature_tolerance_k:
                temperature_converged = True
                break

        # Final synchronous pass — always runs (converged or not) so that the
        # reported flow field and temperature field come from the same solve.
        self._hydraulic_solver._initialise_pressure_field(network_state, case)
        lam_hist, lam_metrics, _ = self._hydraulic_solver._solve_laminar(
            network_state,
            fluid_model,
        )
        turb_hist, turb_metrics, hydraulic_converged = self._hydraulic_solver._solve_turbulent(
            network_state,
            fluid_model,
        )
        all_lam_hist.extend(lam_hist)
        all_lam_metrics.extend(lam_metrics)
        all_turb_hist.extend(turb_hist)
        all_turb_metrics.extend(turb_metrics)
        outer_boundaries.append(len(all_turb_metrics))
        if turb_metrics:
            outer_turb_final.append(turb_metrics[-1])
        final_node_temps, final_pipe_mean_temps = solve_energy_system(
            network_state,
            fluid_model,
            self.convection_scheme,
            T_ref=T_init,
        )
        self._set_node_temperatures(network_state, final_node_temps)
        self._update_component_temperatures(network_state, T_init, final_pipe_mean_temps)

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
            density = fluid_model.density_for_link(link)
            component_flows.append(
                ComponentFlowResult(
                    label=self._hydraulic_solver._component_label(link),
                    mass_flow_kg_per_s=float(link.mass_flow_kg_per_s),
                    volumetric_flow_m3_per_h=float(
                        3600.0 * link.mass_flow_kg_per_s / density
                    ),
                    temperature_in_c=node_temperatures.get(link.start_node.node_id),
                    temperature_out_c=node_temperatures.get(link.end_node.node_id),
                )
            )

        return SolveResult(
            case_name=case.name,
            converged=converged,
            node_pressures_pa=node_pressures,
            node_temperatures_c=node_temperatures,
            component_flows=component_flows,
            laminar_history=all_lam_hist,
            laminar_metrics=all_lam_metrics,
            turbulent_history=all_turb_hist,
            turbulent_metrics=all_turb_metrics,
            temperature_history=temperature_history,
            outer_turbulent_final_metrics=tuple(outer_turb_final),
            outer_iteration_boundaries=tuple(outer_boundaries),
            global_balance=self._compute_global_balance(network_state),
            global_energy_balance=self._compute_global_energy_balance(
                network_state, fluid_model
            ),
            mass_balance_history=mass_balance_history,
            energy_balance_history=energy_balance_history,
        )

    @staticmethod
    def _set_node_temperatures(network_state, node_temperatures: dict[int, float]) -> None:
        for node_id, temperature_c in node_temperatures.items():
            node_state = network_state.nodes[node_id]
            if not node_state.is_thermal_inlet:
                node_state.temperature_c = temperature_c
