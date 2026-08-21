from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Optional

from angelica.closures.convection_scheme import ConvectionScheme, HybridScheme
from angelica.closures.pressure_drop import PressureDropCorrelation
from angelica.core.case import NetworkCase
from angelica.core.network import build_network_state
from angelica.core.state import HeatSourceState, NetworkState, PipeState
from angelica.core.results import ComponentFlowResult, SolveResult
from angelica.core.settings import SolverSettings
from angelica.properties.compositional_fluid import CompositionalFluid
from .base import BaseSolver
from .steady_isothermal_incompressible import SteadyIsothermalIncompressibleSolver


@dataclass
class CompositionalSolverSettings:
    """Settings for the compositional outer iteration loop.

    Each outer iteration: propagate compositions → SIMPLE hydraulics →
    energy equation → temperature update.  Convergence is declared when
    both the maximum relative mixture-density change and the maximum
    nodal temperature change fall below their tolerances.

    Args:
        max_outer_iterations: Maximum outer loop iterations.
        density_rel_tolerance: Max |Δρ_m/ρ_m| across all links for convergence.
        temperature_tolerance_k: Max |ΔT| across all junction nodes (K).
        temperature_relaxation: Under-relaxation factor for temperature updates
            (1.0 = no relaxation).
    """

    max_outer_iterations: int = 50
    density_rel_tolerance: float = 1e-4
    temperature_tolerance_k: float = 0.01
    temperature_relaxation: float = 1.0


class SteadyCompositionalSolver(BaseSolver):
    """Steady-state compositional pipe network solver.

    Models fluid mixtures of arbitrary components using the ``thermo`` library
    for equation-of-state PT flashes at each pipe's local pressure and
    temperature.  Two-phase regions are treated with the homogeneous no-slip
    model (volumetric-fraction-weighted mixture properties).

    ``case.fluid_model`` must be a :class:`~angelica.properties.CompositionalFluid`
    instance, which carries the component list and the default mole-fraction
    vector used before compositions are propagated.

    Multi-inlet composition
    -----------------------
    When ``NetworkCase.inlet_composition_bcs`` is provided, each inlet node
    carries its own mole-fraction boundary condition.  In each outer iteration
    the solver propagates compositions from the inlet nodes through the
    network following the current mass-flow field, mixing streams at
    junctions by molar-flow-weighted average, then writes the resulting
    per-pipe composition to ``PipeState.zs``.

    If ``inlet_composition_bcs`` is empty, every pipe uses
    ``CompositionalFluid.default_zs`` throughout the simulation.

    Algorithm
    ---------
    Outer loop (EOS + temperature + composition coupling):

        1. Propagate compositions from inlets → ``PipeState.zs``.
        2. Snapshot mixture density ρ_m at the current (P, T, z) field.
        3. Run the SIMPLE hydraulic solver (inner loop).
        4. Solve the steady energy equation on the converged flow field.
        5. Update nodal temperatures with under-relaxation.
        6. Compute max |Δρ_m/ρ_m| and max |ΔT|.
        7. Converge when both fall below their tolerances.
    """

    def __init__(
        self,
        hydraulic_settings: Optional[SolverSettings] = None,
        compositional_settings: Optional[CompositionalSolverSettings] = None,
        turbulent_pipe_correlation: Optional[PressureDropCorrelation] = None,
        convection_scheme: Optional[ConvectionScheme] = None,
    ) -> None:
        self.hydraulic_settings = hydraulic_settings or SolverSettings(
            pressure_correction_abs_tolerance_pa=10.0
        )
        self.compositional_settings = compositional_settings or CompositionalSolverSettings()
        self.convection_scheme = convection_scheme or HybridScheme()
        self._hydraulic_solver = SteadyIsothermalIncompressibleSolver(
            settings=self.hydraulic_settings,
            turbulent_pipe_correlation=turbulent_pipe_correlation,
        )

    # ── Composition propagation ───────────────────────────────────────────────

    @staticmethod
    def _propagate_compositions(
        network_state: NetworkState,
        case: NetworkCase,
        default_zs: tuple[float, ...],
    ) -> None:
        """Propagate inlet mole fractions through the network and write to PipeState.zs.

        Seeded at inlet nodes from ``case.inlet_composition_bcs``, compositions
        are propagated downstream following the current mass-flow field.  At
        each junction the incoming streams are mixed by molar-flow-weighted
        average.  Pipes unreachable from an inlet node retain ``default_zs``.

        The algorithm is a fixed-point iteration that converges in at most
        O(network diameter) passes — identical to the black-oil propagator.
        """
        # ── seed inlet node compositions ──────────────────────────────────────
        inlet_ids: set[int] = set()
        node_zs: Dict[int, tuple[float, ...]] = {}
        for ibc in case.inlet_composition_bcs:
            node_zs[ibc.node_id] = tuple(ibc.zs)
            inlet_ids.add(ibc.node_id)

        n_comp = len(default_zs)
        pipe_states = [
            (idx, ps)
            for idx, ps in enumerate(network_state.components)
            if isinstance(ps, PipeState)
        ]

        # ── fixed-point iteration over junction mixing ────────────────────────
        for _ in range(len(network_state.nodes) + 1):
            changed = False
            incoming: Dict[int, list] = defaultdict(list)

            for _, ps in pipe_states:
                mdot = ps.mass_flow_kg_per_s
                if mdot >= 0.0:
                    up_id, down_id = ps.start_node.node_id, ps.end_node.node_id
                else:
                    up_id, down_id = ps.end_node.node_id, ps.start_node.node_id
                    mdot = -mdot
                z_up = node_zs.get(up_id, default_zs)
                incoming[down_id].append((mdot, z_up))

            for nid, contributions in incoming.items():
                if nid in inlet_ids:
                    continue  # Dirichlet BC — never overwrite inlet composition
                if not contributions:
                    continue
                total_m = sum(m for m, _ in contributions)
                if total_m <= 1e-12:
                    continue
                z_mix = tuple(
                    sum(mdot * z[i] for mdot, z in contributions) / total_m
                    for i in range(n_comp)
                )
                prev = node_zs.get(nid)
                node_zs[nid] = z_mix
                if prev is None or max(abs(z_mix[i] - prev[i]) for i in range(n_comp)) > 1e-10:
                    changed = True

            if not changed:
                break

        # ── write per-pipe composition from upstream node ─────────────────────
        for _, ps in pipe_states:
            mdot = ps.mass_flow_kg_per_s
            up_id = ps.start_node.node_id if mdot >= 0.0 else ps.end_node.node_id
            ps.zs = node_zs.get(up_id, default_zs)

    # ── Main solve ────────────────────────────────────────────────────────────

    def solve(self, case: NetworkCase, progress_callback=None) -> SolveResult:
        from angelica.numerics.energy import solve_energy_system

        if not isinstance(case.fluid_model, CompositionalFluid):
            raise TypeError(
                "SteadyCompositionalSolver requires a CompositionalFluid as "
                f"case.fluid_model, got {type(case.fluid_model).__name__}"
            )

        self._require_fixed_temperature(case)

        fluid: CompositionalFluid = case.fluid_model
        default_zs = fluid.default_zs

        # Validate that all inlet BCs have consistent component count
        for ibc in case.inlet_composition_bcs:
            if len(ibc.zs) != len(fluid.component_names):
                raise ValueError(
                    f"InletCompositionBC for node {ibc.node_id} has {len(ibc.zs)} "
                    f"components but CompositionalFluid has {len(fluid.component_names)}"
                )

        network_state = build_network_state(case)
        settings = self.compositional_settings

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

        T_init = self._initial_temperature(case)
        for node in network_state.nodes.values():
            if node.temperature_c is None:
                node.temperature_c = T_init
        for ps in network_state.components:
            if isinstance(ps, (PipeState, HeatSourceState)):
                ps.temperature_c = T_init

        # ── initialise per-pipe compositions ─────────────────────────────────
        for ps in network_state.components:
            if isinstance(ps, PipeState):
                ps.zs = default_zs

        # ── outer loop ────────────────────────────────────────────────────────
        density_history:      list[float] = []
        temperature_history:  list[float] = []
        mass_balance_history: list[float] = []
        energy_balance_history: list[float] = []
        all_lam_hist:   list[float] = []
        all_lam_metrics        = []
        all_turb_hist:  list[float] = []
        all_turb_metrics       = []
        outer_turb_final       = []
        outer_boundaries: list[int] = []
        hydraulic_converged   = False
        density_converged     = False
        temperature_converged = False

        for _outer in range(settings.max_outer_iterations):

            # Propagate compositions from inlets → PipeState.zs
            if case.inlet_composition_bcs:
                self._propagate_compositions(network_state, case, default_zs)

            old_densities = [fluid.density_for_link(link) for link in network_state.components]

            self._hydraulic_solver._initialise_pressure_field(network_state, case)
            lam_hist, lam_metrics, _ = self._hydraulic_solver._solve_laminar(
                network_state, fluid, progress_callback=progress_callback
            )
            turb_hist, turb_metrics, hydraulic_converged = self._hydraulic_solver._solve_turbulent(
                network_state, fluid, progress_callback=progress_callback
            )
            all_lam_hist.extend(lam_hist)
            all_lam_metrics.extend(lam_metrics)
            all_turb_hist.extend(turb_hist)
            all_turb_metrics.extend(turb_metrics)
            outer_boundaries.append(len(all_turb_metrics))
            if turb_metrics:
                outer_turb_final.append(turb_metrics[-1])

            new_node_temps, pipe_mean_temps = solve_energy_system(
                network_state, fluid, self.convection_scheme, T_ref=T_init
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

            _gb  = self._compute_global_balance(network_state)
            _geb = self._compute_global_energy_balance(network_state, fluid)
            mass_balance_history.append(_gb.mass_error_pct)
            energy_balance_history.append(_geb.energy_error_pct)

            new_densities = [fluid.density_for_link(link) for link in network_state.components]
            max_rel_delta = max(
                abs(n - o) / max(abs(o), 1e-10)
                for n, o in zip(new_densities, old_densities)
            )
            density_history.append(max_rel_delta)

            temperature_converged = max_delta_t  < settings.temperature_tolerance_k
            density_converged     = max_rel_delta < settings.density_rel_tolerance
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
            density = fluid.density_for_link(link)
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
                    zs=tuple(link.zs) if getattr(link, "zs", ()) else (),
                )
            )

        # Derive per-node compositions (mass-flow-weighted mix of incoming pipes).
        _incoming: Dict[int, list] = {}
        for link in network_state.components:
            _zs = getattr(link, "zs", ())
            if not _zs:
                continue
            _mdot = link.mass_flow_kg_per_s
            _nid = link.end_node.node_id if _mdot >= 0.0 else link.start_node.node_id
            _incoming.setdefault(_nid, []).append((abs(_mdot), tuple(_zs)))
        node_compositions: Dict[int, tuple] = {}
        for _nid, _arr in _incoming.items():
            _tot = sum(m for m, _ in _arr)
            if _tot <= 0.0:
                continue
            _nc = len(_arr[0][1])
            node_compositions[_nid] = tuple(
                sum(m * zs[i] for m, zs in _arr) / _tot for i in range(_nc)
            )

        return SolveResult(
            case_name=case.name,
            converged=(hydraulic_converged and density_converged and temperature_converged),
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
            global_balance=self._compute_global_balance(network_state),
            global_energy_balance=self._compute_global_energy_balance(network_state, fluid),
            mass_balance_history=mass_balance_history,
            energy_balance_history=energy_balance_history,
            component_names=fluid.component_names,
            component_mws=fluid.component_mws,
            node_compositions=node_compositions,
        )
