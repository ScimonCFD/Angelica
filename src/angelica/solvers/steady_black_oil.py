from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Dict

from angelica.closures.convection_scheme import ConvectionScheme, HybridScheme
from angelica.closures.pressure_drop import PressureDropCorrelation
from angelica.core.case import NetworkCase
from angelica.core.network import build_network_state
from angelica.core.state import PipeState, HeatSourceState, NetworkState
from angelica.core.results import ComponentFlowResult, SolveResult
from angelica.core.settings import SolverSettings
from angelica.properties.black_oil import BlackOilComposition, compute_pvt
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

    Multi-inlet composition
    -----------------------
    When ``NetworkCase.inlet_fluid_bcs`` is provided, each inlet node
    carries its own four-parameter composition (API, gas gravity, GOR, WOR).
    In each outer iteration the solver propagates compositions from the inlet
    nodes through the network following the flow direction, mixing streams at
    junctions by mass-weighted average.  The result is a per-pipe composition
    used to evaluate density and viscosity.

    If ``inlet_fluid_bcs`` is empty the solver falls back to the single global
    ``fluid_model`` composition (backward-compatible behaviour).

    Algorithm
    ---------
    Outer loop (PVT + temperature + composition coupling):
        1. Propagate compositions from inlets → per-pipe ``BlackOilComposition``.
        2. Snapshot mixture density ρ_m at the current (P, T, composition) field.
        3. Run the SIMPLE hydraulic solver (inner loop).
        4. Solve the steady energy equation on the converged flow field.
        5. Update nodal temperatures with under-relaxation.
        6. Compute max |Δρ_m/ρ_m| and max |ΔT|.
        7. Converge when both fall below their tolerances.
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

    # ── Composition helpers ───────────────────────────────────────────────────

    @staticmethod
    def _propagate_compositions(
        network_state: NetworkState,
        case: NetworkCase,
        default_comp: BlackOilComposition,
    ) -> Dict[int, BlackOilComposition]:
        """Return a per-pipe-index composition map after one propagation pass.

        Compositions are seeded at inlet nodes and propagated downstream
        following the current flow field.  At each junction the incoming
        streams are mixed by mass-weighted average.  Pipes in which flow
        is ambiguous or zero inherit the default composition.

        Returns
        -------
        dict mapping index-in-network_state.components → BlackOilComposition
        """
        # ── seed inlet node compositions ──────────────────────────────────────
        node_comp: Dict[int, BlackOilComposition] = {
            ibc.node_id: BlackOilComposition(
                api_gravity      = ibc.api_gravity,
                gas_gravity      = ibc.gas_gravity,
                gor_sc_m3_per_m3 = ibc.gor_sc_m3_per_m3,
                wor_sc_m3_per_m3 = ibc.wor_sc_m3_per_m3,
            )
            for ibc in case.inlet_fluid_bcs
        }

        # ── determine flow directions from current mass-flow field ────────────
        # For each node accumulate (mass_flow, composition) of incoming pipes.
        # Incoming means: flow is from the pipe's upstream end TO this node.
        incoming: Dict[int, list] = defaultdict(list)
        pipe_states = [
            (idx, ps)
            for idx, ps in enumerate(network_state.components)
            if isinstance(ps, PipeState)
        ]

        # ── first pass: propagate from seeded nodes into their downstream pipes
        # We iterate several times to handle networks where compositions must
        # travel through multiple junctions (convergence in O(diameter) passes).
        for _ in range(len(network_state.nodes) + 1):
            changed = False
            incoming = defaultdict(list)

            for _, ps in pipe_states:
                ṁ = ps.mass_flow_kg_per_s
                if ṁ >= 0.0:
                    up_id, down_id = ps.start_node.node_id, ps.end_node.node_id
                else:
                    up_id, down_id = ps.end_node.node_id, ps.start_node.node_id
                    ṁ = -ṁ

                comp = node_comp.get(up_id, default_comp)
                incoming[down_id].append((ṁ, comp))

            # Mix at each non-inlet junction
            for nid, contributions in incoming.items():
                if nid in {ibc.node_id for ibc in case.inlet_fluid_bcs}:
                    continue  # inlet nodes keep their prescribed composition
                if not contributions:
                    continue
                total_m = sum(m for m, _ in contributions)
                if total_m <= 1e-12:
                    continue
                mixed = contributions[0][1]
                w_acc = contributions[0][0]
                for m_i, c_i in contributions[1:]:
                    mixed = mixed.mix(c_i, w_acc, m_i)
                    w_acc += m_i
                prev = node_comp.get(nid)
                node_comp[nid] = mixed
                if prev is None or (
                    abs(mixed.api_gravity - prev.api_gravity) > 1e-9
                    or abs(mixed.gas_gravity - prev.gas_gravity) > 1e-9
                    or abs(mixed.gor_sc_m3_per_m3 - prev.gor_sc_m3_per_m3) > 1e-9
                    or abs(mixed.wor_sc_m3_per_m3 - prev.wor_sc_m3_per_m3) > 1e-9
                ):
                    changed = True

            if not changed:
                break

        # ── assign per-pipe composition from upstream node ────────────────────
        pipe_comp: Dict[int, BlackOilComposition] = {}
        for idx, ps in pipe_states:
            ṁ = ps.mass_flow_kg_per_s
            if ṁ >= 0.0:
                up_id = ps.start_node.node_id
            else:
                up_id = ps.end_node.node_id
            pipe_comp[idx] = node_comp.get(up_id, default_comp)

        # Non-pipe components inherit the default
        for idx in range(len(network_state.components)):
            if idx not in pipe_comp:
                pipe_comp[idx] = default_comp

        return pipe_comp

    @staticmethod
    def _pipe_density(ps: PipeState, comp: BlackOilComposition) -> float:
        P = 0.5 * ((ps.start_node.pressure_pa or 0.0) + (ps.end_node.pressure_pa or 0.0))
        T = ps.temperature_c if ps.temperature_c is not None else 20.0
        return compute_pvt(P, T, comp.api_gravity, comp.gas_gravity,
                           comp.gor_sc_m3_per_m3, comp.wor_sc_m3_per_m3
                           ).mixture_density_kg_per_m3

    @staticmethod
    def _pipe_viscosity(ps: PipeState, comp: BlackOilComposition) -> float:
        P = 0.5 * ((ps.start_node.pressure_pa or 0.0) + (ps.end_node.pressure_pa or 0.0))
        T = ps.temperature_c if ps.temperature_c is not None else 20.0
        return compute_pvt(P, T, comp.api_gravity, comp.gas_gravity,
                           comp.gor_sc_m3_per_m3, comp.wor_sc_m3_per_m3
                           ).mixture_viscosity_pa_s

    # ── FluidModel proxy ──────────────────────────────────────────────────────

    @staticmethod
    def _make_pipe_fluid(comp: BlackOilComposition, ref_p: float, ref_t: float):
        """Wrap a per-pipe composition into a minimal FluidModel for the hydraulic solver."""
        from angelica.properties.black_oil import BlackOilFluid
        return BlackOilFluid(
            api_gravity             = comp.api_gravity,
            gas_gravity             = comp.gas_gravity,
            gor_sc_m3_per_m3        = comp.gor_sc_m3_per_m3,
            wor_sc_m3_per_m3        = comp.wor_sc_m3_per_m3,
            reference_pressure_pa   = ref_p,
            reference_temperature_c = ref_t,
        )

    # ── Mixed-composition FluidModel ──────────────────────────────────────────

    def _build_mixed_fluid(
        self,
        network_state: NetworkState,
        pipe_comps: Dict[int, BlackOilComposition],
        default_comp: BlackOilComposition,
        ref_t: float,
    ):
        """Return a FluidModel-compatible object that dispatches per link."""
        from angelica.properties.black_oil import BlackOilFluid

        # Build a lookup: component_object_id → per-pipe fluid
        id_to_comp: Dict[int, BlackOilComposition] = {}
        for idx, ps in enumerate(network_state.components):
            id_to_comp[id(ps)] = pipe_comps.get(idx, default_comp)

        ref_p = 5e6  # typical midpoint pressure
        default_fluid = BlackOilFluid(
            api_gravity             = default_comp.api_gravity,
            gas_gravity             = default_comp.gas_gravity,
            gor_sc_m3_per_m3        = default_comp.gor_sc_m3_per_m3,
            wor_sc_m3_per_m3        = default_comp.wor_sc_m3_per_m3,
            reference_pressure_pa   = ref_p,
            reference_temperature_c = ref_t,
        )

        class _PerPipeFluid:
            def _fluid_for(self, link_state) -> BlackOilFluid:
                comp = id_to_comp.get(id(link_state), default_comp)
                return BlackOilFluid(
                    api_gravity             = comp.api_gravity,
                    gas_gravity             = comp.gas_gravity,
                    gor_sc_m3_per_m3        = comp.gor_sc_m3_per_m3,
                    wor_sc_m3_per_m3        = comp.wor_sc_m3_per_m3,
                    reference_pressure_pa   = ref_p,
                    reference_temperature_c = ref_t,
                )

            def density_for_link(self, link_state) -> float:
                return self._fluid_for(link_state).density_for_link(link_state)

            def viscosity_for_link(self, link_state) -> float:
                return self._fluid_for(link_state).viscosity_for_link(link_state)

            def specific_heat_for_link(self, link_state) -> float:
                return self._fluid_for(link_state).specific_heat_for_link(link_state)

            def thermal_conductivity_for_link(self, link_state) -> float:
                return self._fluid_for(link_state).thermal_conductivity_for_link(link_state)

        return _PerPipeFluid() if id_to_comp else default_fluid

    # ── Main solve ────────────────────────────────────────────────────────────

    def solve(self, case: NetworkCase, progress_callback=None) -> SolveResult:
        from angelica.numerics.energy import solve_energy_system

        self._require_fixed_temperature(case)
        network_state = build_network_state(case)
        settings      = self.black_oil_settings
        fluid_model   = case.fluid_model

        # Determine whether per-inlet composition is active
        use_per_inlet = bool(case.inlet_fluid_bcs)
        if use_per_inlet:
            # Build default composition from the global fluid_model as fallback
            from angelica.properties.black_oil import BlackOilFluid
            if isinstance(fluid_model, BlackOilFluid):
                default_comp = BlackOilComposition(
                    api_gravity      = fluid_model.api_gravity,
                    gas_gravity      = fluid_model.gas_gravity,
                    gor_sc_m3_per_m3 = fluid_model.gor_sc_m3_per_m3,
                    wor_sc_m3_per_m3 = fluid_model.wor_sc_m3_per_m3,
                )
            else:
                # Use first inlet's composition as default
                ibc0 = case.inlet_fluid_bcs[0]
                default_comp = BlackOilComposition(
                    api_gravity      = ibc0.api_gravity,
                    gas_gravity      = ibc0.gas_gravity,
                    gor_sc_m3_per_m3 = ibc0.gor_sc_m3_per_m3,
                    wor_sc_m3_per_m3 = ibc0.wor_sc_m3_per_m3,
                )
        else:
            default_comp = None

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

        # ── outer loop ────────────────────────────────────────────────────────
        density_history:     list[float] = []
        temperature_history: list[float] = []
        mass_balance_history: list[float] = []
        energy_balance_history: list[float] = []
        all_lam_hist:   list[float] = []
        all_lam_metrics        = []
        all_turb_hist:  list[float] = []
        all_turb_metrics       = []
        outer_turb_final       = []
        outer_boundaries: list[int] = []
        hydraulic_converged    = False
        density_converged      = False
        temperature_converged  = False

        # pipe_comps: per-pipe composition (None → use global fluid_model)
        pipe_comps: Dict[int, BlackOilComposition] | None = None

        for _outer in range(settings.max_outer_iterations):

            # ── build effective fluid model for this iteration ────────────────
            if use_per_inlet:
                pipe_comps = self._propagate_compositions(
                    network_state, case, default_comp  # type: ignore[arg-type]
                )
                effective_fluid = self._build_mixed_fluid(
                    network_state, pipe_comps, default_comp, T_init  # type: ignore[arg-type]
                )
            else:
                effective_fluid = fluid_model

            old_densities = [
                effective_fluid.density_for_link(link)
                for link in network_state.components
            ]

            self._hydraulic_solver._initialise_pressure_field(network_state, case)
            lam_hist, lam_metrics, _ = self._hydraulic_solver._solve_laminar(
                network_state, effective_fluid, progress_callback=progress_callback
            )
            turb_hist, turb_metrics, hydraulic_converged = self._hydraulic_solver._solve_turbulent(
                network_state, effective_fluid, progress_callback=progress_callback
            )
            all_lam_hist.extend(lam_hist)
            all_lam_metrics.extend(lam_metrics)
            all_turb_hist.extend(turb_hist)
            all_turb_metrics.extend(turb_metrics)
            outer_boundaries.append(len(all_turb_metrics))
            if turb_metrics:
                outer_turb_final.append(turb_metrics[-1])

            new_node_temps, pipe_mean_temps = solve_energy_system(
                network_state, effective_fluid, self.convection_scheme, T_ref=T_init
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

            _gb = self._compute_global_balance(network_state)
            mass_balance_history.append(_gb.mass_error_pct)
            _geb = self._compute_global_energy_balance(network_state, effective_fluid)
            energy_balance_history.append(abs(_geb.energy_error_kw))

            new_densities = [
                effective_fluid.density_for_link(link)
                for link in network_state.components
            ]
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
        for idx, link in enumerate(network_state.components):
            if use_per_inlet and pipe_comps is not None:
                comp = pipe_comps.get(idx, default_comp)  # type: ignore[arg-type]
                density = effective_fluid.density_for_link(link)
            else:
                comp    = None
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
            global_balance=self._compute_global_balance(network_state),
            global_energy_balance=self._compute_global_energy_balance(
                network_state, effective_fluid
            ),
            mass_balance_history=mass_balance_history,
            energy_balance_history=energy_balance_history,
        )
