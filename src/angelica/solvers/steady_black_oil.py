from __future__ import annotations

from dataclasses import dataclass

from angelica.closures.convection_scheme import ConvectionScheme, HybridScheme
from angelica.closures.pressure_drop import PressureDropCorrelation
from angelica.core.case import NetworkCase
from angelica.core.network import build_network_state
from angelica.core.state import PipeState, HeatSourceState
from angelica.core.results import ComponentFlowResult, SolveResult
from angelica.core.settings import SolverSettings
from .base import BaseSolver
from .steady_isothermal_incompressible import SteadyIsothermalIncompressibleSolver


@dataclass
class BlackOilSolverSettings:
    """Settings for the black-oil outer iteration loop.

    Each outer iteration: update mixture PVT → SIMPLE hydraulics →
    energy equation → temperature update.  Convergence is declared when
    both the maximum relative mixture-density change and the maximum
    nodal temperature change fall below their tolerances.

    Args:
        max_outer_iterations: Maximum outer loop iterations.
        density_rel_tolerance: Max |Δρ_m/ρ_m| across all links for
            convergence.
        temperature_tolerance_k: Max |ΔT| across all junction nodes (K)
            for convergence.
        temperature_relaxation: Under-relaxation factor for temperature
            updates (1.0 = no relaxation).
    """

    max_outer_iterations: int = 50
    density_rel_tolerance: float = 1e-4
    temperature_tolerance_k: float = 0.01
    temperature_relaxation: float = 1.0


class SteadyBlackOilSolver(BaseSolver):
    """Steady-state black-oil three-phase pipe network solver.

    Models gas, oil, and water as a homogeneous (no-slip) mixture whose
    density and viscosity are computed from black-oil PVT correlations
    (Standing, Hall-Yarborough, Beggs-Robinson, Lee-Gonzalez-Eakin) at
    the local pressure and temperature.

    Algorithm
    ---------
    Outer loop (PVT + temperature coupling):
        1. Snapshot mixture density ρ_m at the current (P, T) field.
        2. Run the SIMPLE hydraulic solver (inner loop).
           BlackOilFluid.density_for_link re-evaluates ρ_m(P, T) as the
           pressure field evolves.
        3. Solve the steady energy equation on the converged flow field.
        4. Update nodal temperatures with under-relaxation.
        5. Compute max |Δρ_m/ρ_m| and max |ΔT|.
        6. Converge when both fall below their tolerances.

    Applicability
    -------------
    - Three-phase (oil + dissolved/free gas + water) steady-state flow.
    - No-slip homogeneous mixture — valid when phase velocities are
      similar (high flow rates, turbulent flow).
    - Valid for both undersaturated (P > Pb, no free gas) and two-phase
      (P < Pb, free gas present) conditions.
    - Uses a single GOR and WOR for the entire network (uniform inlet
      composition); junction mixing with multiple different compositions
      is deferred to a future version.
    """

    def __init__(
        self,
        hydraulic_settings: SolverSettings | None = None,
        black_oil_settings: BlackOilSolverSettings | None = None,
        turbulent_pipe_correlation: PressureDropCorrelation | None = None,
        convection_scheme: ConvectionScheme | None = None,
    ) -> None:
        self.hydraulic_settings = hydraulic_settings or SolverSettings(
            pressure_correction_abs_tolerance_pa=10.0
        )
        self.black_oil_settings = black_oil_settings or BlackOilSolverSettings()
        self.convection_scheme = convection_scheme or HybridScheme()
        self._hydraulic_solver = SteadyIsothermalIncompressibleSolver(
            settings=self.hydraulic_settings,
            turbulent_pipe_correlation=turbulent_pipe_correlation,
        )

    def solve(self, case: NetworkCase, progress_callback=None) -> SolveResult:
        from angelica.numerics.energy import solve_energy_system

        network_state = build_network_state(case)
        settings      = self.black_oil_settings
        fluid_model   = case.fluid_model

        # ── initialise thermal BCs ────────────────────────────────────────────
        for tb in case.thermal_inlets:
            if tb.node_id not in network_state.nodes:
                continue
            node_st = network_state.nodes[tb.node_id]
            if tb.bc_type == "fixed_temperature":
                node_st.temperature_c    = tb.temperature_c
                node_st.is_thermal_inlet = True
            else:
                node_st.thermal_gradient_dc_per_m = tb.gradient_dc_per_m

        # ── initialise temperature field ──────────────────────────────────────
        T_init = self._initial_temperature(case)
        for node in network_state.nodes.values():
            if node.temperature_c is None:
                node.temperature_c = T_init
        for ps in network_state.components:
            if isinstance(ps, (PipeState, HeatSourceState)):
                ps.temperature_c = T_init

        # ── outer loop ────────────────────────────────────────────────────────
        density_history:     list[float] = []
        temperature_history: list[float] = []
        all_lam_hist:   list[float] = []
        all_lam_metrics        = []
        all_turb_hist:  list[float] = []
        all_turb_metrics       = []
        outer_turb_final       = []
        outer_boundaries: list[int] = []
        hydraulic_converged    = False
        density_converged      = False
        temperature_converged  = False

        for _outer in range(settings.max_outer_iterations):
            old_densities = [fluid_model.density_for_link(link) for link in network_state.components]

            self._hydraulic_solver._initialise_pressure_field(network_state, case)
            lam_hist, lam_metrics, _ = self._hydraulic_solver._solve_laminar(
                network_state, fluid_model, progress_callback=progress_callback
            )
            turb_hist, turb_metrics, hydraulic_converged = self._hydraulic_solver._solve_turbulent(
                network_state, fluid_model, progress_callback=progress_callback
            )
            all_lam_hist.extend(lam_hist)
            all_lam_metrics.extend(lam_metrics)
            all_turb_hist.extend(turb_hist)
            all_turb_metrics.extend(turb_metrics)
            outer_boundaries.append(len(all_turb_metrics))
            if turb_metrics:
                outer_turb_final.append(turb_metrics[-1])

            new_node_temps, pipe_mean_temps = solve_energy_system(
                network_state, fluid_model, self.convection_scheme, T_ref=T_init
            )

            max_delta_t = 0.0
            alpha = settings.temperature_relaxation
            for nid, T_new in new_node_temps.items():
                T_old = network_state.nodes[nid].temperature_c or T_init
                delta = alpha * (T_new - T_old)
                max_delta_t = max(max_delta_t, abs(delta))
                if not network_state.nodes[nid].is_thermal_inlet:
                    network_state.nodes[nid].temperature_c = T_old + delta
            temperature_history.append(max_delta_t)
            self._update_component_temperatures(network_state, T_init, pipe_mean_temps)

            new_densities = [fluid_model.density_for_link(link) for link in network_state.components]
            max_rel_delta = max(
                abs(n - o) / max(abs(o), 1e-10)
                for n, o in zip(new_densities, old_densities)
            )
            density_history.append(max_rel_delta)

            temperature_converged = max_delta_t   < settings.temperature_tolerance_k
            density_converged     = max_rel_delta  < settings.density_rel_tolerance
            if temperature_converged and density_converged:
                break

        # ── build result ──────────────────────────────────────────────────────
        node_pressures = {
            nid: float(network_state.nodes[nid].pressure_pa)
            for nid in sorted(network_state.nodes)
        }
        node_temperatures = {
            nid: float(network_state.nodes[nid].temperature_c or T_init)
            for nid in sorted(network_state.nodes)
        }
        component_flows = []
        for link in network_state.components:
            density = fluid_model.density_for_link(link)
            t_in  = node_temperatures.get(link.start_node.node_id)
            t_out = node_temperatures.get(link.end_node.node_id)
            component_flows.append(
                ComponentFlowResult(
                    label=self._hydraulic_solver._component_label(link),
                    mass_flow_kg_per_s=float(link.mass_flow_kg_per_s),
                    volumetric_flow_m3_per_h=float(
                        3600.0 * link.mass_flow_kg_per_s / density
                    ),
                    temperature_in_c=t_in,
                    temperature_out_c=t_out,
                )
            )

        return SolveResult(
            case_name=case.name,
            converged=hydraulic_converged and density_converged and temperature_converged,
            node_pressures_pa=node_pressures,
            node_temperatures_c=node_temperatures,
            component_flows=component_flows,
            laminar_history=all_lam_hist,
            laminar_metrics=all_lam_metrics,
            turbulent_history=all_turb_hist,
            turbulent_metrics=all_turb_metrics,
            temperature_history=temperature_history,
            density_history=density_history,
            outer_turbulent_final_metrics=tuple(outer_turb_final),
            outer_iteration_boundaries=tuple(outer_boundaries),
        )
