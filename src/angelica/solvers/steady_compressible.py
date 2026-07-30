from __future__ import annotations

from dataclasses import dataclass

from angelica.closures.pressure_drop import PressureDropCorrelation
from angelica.core.network import build_network_state
from angelica.core.results import ComponentFlowResult, SolveResult
from angelica.core.settings import SolverSettings
from .base import BaseSolver
from .steady_isothermal_incompressible import SteadyIsothermalIncompressibleSolver


@dataclass
class CompressibleSolverSettings:
    """Settings for the outer density iteration loop.

    The compressible solver wraps SIMPLE in an outer loop that re-evaluates
    density from the equation of state at the updated pressure field after
    each full hydraulic solve.  The loop repeats until the maximum relative
    density change across all links falls below density_rel_tolerance.

    Args:
        max_density_iterations: Maximum number of outer density iterations.
        density_rel_tolerance: Convergence criterion — max |Δρ/ρ| across all
            links between successive outer iterations.
    """

    max_density_iterations: int = 50
    density_rel_tolerance: float = 1e-4


class SteadyCompressibleSolver(BaseSolver):
    """Steady, single-phase compressible pipe network solver.

    The fluid density is computed from an EquationOfState as ρ(P, T) via
    CompressibleFluid.  Viscosity, specific heat, and thermal conductivity
    are treated as functions of temperature only.

    Algorithm
    ---------
    The solver uses a segregated outer-inner structure:

    Outer loop (density):
        1. Snapshot ρ at the current pressure field.
        2. Run the SIMPLE hydraulic solver to convergence (inner loop).
           Within each SIMPLE iteration, CompressibleFluid.density_for_link
           re-evaluates ρ(P, T) at the current node pressures automatically.
        3. Compute max relative density change: max |Δρ/ρ|.
        4. If max |Δρ/ρ| < density_rel_tolerance → converged.  Otherwise repeat.

    Applicability
    -------------
    - Single-phase gas or mildly compressible liquid where ρ = ρ(P, T).
    - Steady-state, single-component flow.
    - Not valid for two-phase flow or near the critical point where the EOS
      may return non-physical roots.
    """

    def __init__(
        self,
        hydraulic_settings: SolverSettings | None = None,
        compressible_settings: CompressibleSolverSettings | None = None,
        turbulent_pipe_correlation: PressureDropCorrelation | None = None,
    ) -> None:
        # Gas systems operate at 100 kPa – 10 MPa; a 10 Pa absolute correction
        # tolerance is < 0.01 % of the lowest typical pressure — tight enough
        # in practice while avoiding oscillation near convergence.
        self.hydraulic_settings = hydraulic_settings or SolverSettings(
            pressure_correction_abs_tolerance_pa=10.0
        )
        self.compressible_settings = compressible_settings or CompressibleSolverSettings()
        self._hydraulic_solver = SteadyIsothermalIncompressibleSolver(
            settings=self.hydraulic_settings,
            turbulent_pipe_correlation=turbulent_pipe_correlation,
        )

    def solve(self, case, progress_callback=None) -> SolveResult:
        network_state = build_network_state(case)
        settings = self.compressible_settings
        fluid_model = case.fluid_model

        density_history: list[float] = []
        lam_hist: list[float] = []
        lam_metrics = []
        turb_hist: list[float] = []
        turb_metrics = []
        hydraulic_converged = False
        density_converged = False

        for _outer in range(settings.max_density_iterations):
            old_densities = [fluid_model.density_for_link(link) for link in network_state.components]

            self._hydraulic_solver._initialise_pressure_field(network_state, case)
            lam_hist, lam_metrics, _ = self._hydraulic_solver._solve_laminar(
                network_state, fluid_model
            )
            turb_hist, turb_metrics, hydraulic_converged = self._hydraulic_solver._solve_turbulent(
                network_state, fluid_model
            )

            new_densities = [fluid_model.density_for_link(link) for link in network_state.components]

            max_rel_delta = max(
                abs(n - o) / max(abs(o), 1e-10)
                for n, o in zip(new_densities, old_densities)
            )
            density_history.append(max_rel_delta)

            if max_rel_delta < settings.density_rel_tolerance:
                density_converged = True
                break

        node_pressures = {
            nid: float(network_state.nodes[nid].pressure_pa)
            for nid in sorted(network_state.nodes)
        }
        component_flows = []
        for link in network_state.components:
            density = fluid_model.density_for_link(link)
            component_flows.append(
                ComponentFlowResult(
                    label=self._hydraulic_solver._component_label(link),
                    mass_flow_kg_per_s=float(link.mass_flow_kg_per_s),
                    volumetric_flow_m3_per_h=float(3600.0 * link.mass_flow_kg_per_s / density),
                )
            )

        return SolveResult(
            case_name=case.name,
            converged=hydraulic_converged and density_converged,
            node_pressures_pa=node_pressures,
            component_flows=component_flows,
            laminar_history=lam_hist,
            laminar_metrics=lam_metrics,
            turbulent_history=turb_hist,
            turbulent_metrics=turb_metrics,
            density_history=density_history,
        )
