from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from angelica.closures.choke import ChokeFlowModel
from angelica.closures.convection_scheme import ConvectionScheme, HybridScheme
from angelica.closures.pressure_drop import PressureDropCorrelation
from angelica.core.case import NetworkCase
from angelica.core.network import build_network_state
from angelica.core.results import ComponentFlowResult, SolveResult
from angelica.core.settings import SolverSettings
from angelica.core.state import HeatSourceState, NetworkState, PipeState
from angelica.properties.compositional_fluid import CompositionalFluid

from .base import BaseSolver
from .separators import SEPARATOR_FRACTION_TOLERANCE, SeparatorNetworkMixin
from .steady_isothermal_incompressible import SteadyIsothermalIncompressibleSolver
from .wells import WellNetworkMixin


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
        composition_relaxation: Under-relaxation factor for junction composition
            updates (1.0 = no relaxation).  Values below 1.0 damp composition
            oscillations in looped networks with mixing junctions; 0.5 is a
            safe starting point when a 2-cycle is observed.
        pipe_enthalpy_terms: Joule-Thomson cooling, latent heat and potential
            energy along pipes, with enthalpy-conserving junctions (default
            True).  False drops them, for comparison with the pre-1.18 model
            (valves stay isenthalpic); the energy balance is still computed
            with real enthalpy, so it then shows what that model leaves out.
    """

    max_outer_iterations: int = 50
    density_rel_tolerance: float = 1e-4
    temperature_tolerance_k: float = 0.01
    temperature_relaxation: float = 1.0
    composition_relaxation: float = 1.0
    pipe_enthalpy_terms: bool = True


class SteadyCompositionalSolver(WellNetworkMixin, SeparatorNetworkMixin, BaseSolver):
    """Steady-state compositional pipe network solver.

    Models fluid mixtures of arbitrary components using the ``thermo`` library
    for equation-of-state PT flashes at each pipe's local pressure and
    temperature.  Two-phase regions are treated with the homogeneous no-slip
    model (volumetric-fraction-weighted mixture properties) unless a slip
    correlation is passed as ``turbulent_pipe_correlation``: with
    ``BeggsBrillCorrelation()`` the pipes use Beggs-Brill holdup and friction
    with the gas and liquid (condensate + free water) properties of each
    pipe's flash.

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
        hydraulic_settings: SolverSettings | None = None,
        compositional_settings: CompositionalSolverSettings | None = None,
        turbulent_pipe_correlation: PressureDropCorrelation | None = None,
        convection_scheme: ConvectionScheme | None = None,
    ) -> None:
        self.hydraulic_settings = hydraulic_settings or SolverSettings(
            pressure_correction_abs_tolerance_pa=10.0
        )
        self.compositional_settings = compositional_settings or CompositionalSolverSettings()
        self.convection_scheme = convection_scheme or HybridScheme()
        self._hydraulic_solver = SteadyIsothermalIncompressibleSolver(
            settings=self.hydraulic_settings,
            turbulent_pipe_correlation=turbulent_pipe_correlation,
            # K-factor fittings as before; also IPR (WellInflow) and IPR +
            # tubing (Well) links with the compositional fluid.
            fitting_correlation=ChokeFlowModel(),
        )

    # ── Composition propagation ───────────────────────────────────────────────

    @staticmethod
    def _propagate_compositions(
        network_state: NetworkState,
        case: NetworkCase,
        default_zs: tuple[float, ...],
        node_zs_prev: dict[int, tuple[float, ...]] | None = None,
        relaxation: float = 1.0,
        component_mws: tuple[float, ...] | None = None,
        outgoing=None,
    ) -> dict[int, tuple[float, ...]]:
        """Propagate inlet mole fractions through the network and write to PipeState.zs.

        Seeded at inlet nodes from ``case.inlet_composition_bcs``, compositions
        are propagated downstream following the current mass-flow field.  At
        each junction the incoming streams are mixed by molar flow, ṁ/MW(z)
        (mixing mole fractions by mass flow does not conserve the species when
        the streams' molar masses differ).  Pipes unreachable from an inlet
        node retain ``default_zs``.

        When ``relaxation < 1.0`` and ``node_zs_prev`` is supplied, the newly
        computed junction composition is blended with the value from the
        previous outer iteration:
            z_eff = relaxation * z_new + (1 - relaxation) * z_prev
        This damps composition oscillations in looped networks with mixing
        junctions (classic 2-cycle instability).

        ``outgoing`` (separators, see
        ``SeparatorNetworkMixin._separator_outgoing``) gives the composition
        of a link leaving a separator node.

        Returns the per-node composition dict so the caller can store it for
        the next outer iteration.
        """
        # ── seed inlet node compositions ──────────────────────────────────────
        inlet_ids: set[int] = set()
        node_zs: dict[int, tuple[float, ...]] = {}
        for ibc in case.inlet_composition_bcs:
            node_zs[ibc.node_id] = tuple(ibc.zs)
            inlet_ids.add(ibc.node_id)

        n_comp = len(default_zs)
        # Every link type carries the fluid: pipes, heat sources, valves,
        # chokes and pumps (a fitting must not cut the propagation).
        pipe_states = list(enumerate(network_state.components))

        # ── fixed-point iteration over junction mixing ────────────────────────
        for _ in range(len(network_state.nodes) + 1):
            changed = False
            incoming: dict[int, list] = defaultdict(list)

            for _, ps in pipe_states:
                mdot = ps.mass_flow_kg_per_s
                if mdot >= 0.0:
                    up_id, down_id = ps.start_node.node_id, ps.end_node.node_id
                else:
                    up_id, down_id = ps.end_node.node_id, ps.start_node.node_id
                    mdot = -mdot
                z_up = outgoing(ps, up_id, node_zs, default_zs) if outgoing is not None else None
                if z_up is None:
                    z_up = node_zs.get(up_id, default_zs)
                if component_mws:
                    mw = sum(z * m for z, m in zip(z_up, component_mws))
                    mdot = mdot / mw if mw > 0.0 else mdot     # molar flow (kmol/s)
                incoming[down_id].append((mdot, z_up))
            for nid_s, node_s in network_state.nodes.items():
                for m_s, z_s in node_s.source_streams:        # lift gas (mole fractions)
                    if component_mws:
                        mw_s = sum(z * m for z, m in zip(z_s, component_mws))
                        m_s = m_s / mw_s if mw_s > 0.0 else m_s
                    incoming[nid_s].append((m_s, tuple(z_s)))

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
                # Inter-outer-iteration relaxation to damp mixing oscillations
                if relaxation < 1.0 and node_zs_prev is not None:
                    z_old = node_zs_prev.get(nid)
                    if z_old is not None and len(z_old) == n_comp:
                        z_mix = tuple(
                            relaxation * z_mix[i] + (1.0 - relaxation) * z_old[i]
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
            z_out = outgoing(ps, up_id, node_zs, default_zs) if outgoing is not None else None
            ps.zs = z_out if z_out is not None else node_zs.get(up_id, default_zs)

        return node_zs

    # ── Main solve ────────────────────────────────────────────────────────────

    @staticmethod
    def _hydrate_report(link, fluid, node_temperatures) -> dict:
        """Hydrate formation temperature and subcooling at the pipe end where
        the fluid is coldest relative to the formation curve.  Empty (fields
        stay None) for links that are not pipes, carry no hydrate former, or
        where no hydrate is stable above 0 °C at either end's pressure."""
        from angelica.properties.hydrate import hydrate_formation_temperature
        if not isinstance(link, PipeState) or abs(link.mass_flow_kg_per_s) < 1e-12:
            return {}
        zs = getattr(link, "zs", ()) or fluid.default_zs
        worst: dict = {}
        try:
            for node in (link.start_node, link.end_node):
                T = node_temperatures.get(node.node_id)
                if T is None or node.pressure_pa is None:
                    continue
                pt = hydrate_formation_temperature(fluid.component_names, zs, float(node.pressure_pa),
                                                   fluid.eos_name)
                if pt is None:
                    continue                  # no hydrate stable above 0 °C at this end
                sub = pt.temperature_c - T
                if not worst or sub > worst["hydrate_subcooling_k"]:
                    worst = {"hydrate_subcooling_k": sub, "hydrate_temperature_c": pt.temperature_c,
                             "hydrate_structure": pt.structure}
        except Exception:      # a check must never break a converged solve
            return {}
        return worst

    def _multiphase_report(self, link, fluid):
        """Holdup, flow pattern and velocities of a solved pipe (None otherwise)."""
        from angelica.closures.beggs_brill import multiphase_flow_report
        try:
            return multiphase_flow_report(link, fluid, self._hydraulic_solver.turbulent_pipe_correlation)
        except Exception:  # a report must never break a converged solve
            return None

    def _solve_once(self, case: NetworkCase, progress_callback=None,
                    latch_disabled: frozenset[str] = frozenset()) -> SolveResult:
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
        self._disable_well_latch(network_state, latch_disabled)
        self._setup_separators(network_state)
        settings = self.compositional_settings

        # Fluid-aware correlations (Beggs-Brill, wells) read phase properties from the fluid.
        for corr in (self._hydraulic_solver.turbulent_pipe_correlation,
                     self._hydraulic_solver.fitting_correlation):
            if hasattr(corr, "set_fluid"):
                corr.set_fluid(fluid)

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

        # ── initialise per-link compositions ─────────────────────────────────
        for ps in network_state.components:
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
        node_zs_outer: dict[int, tuple[float, ...]] = {}

        from angelica.closures.separator import compositional_split

        def split_fn(node, zs):
            T = node.temperature_c if node.temperature_c is not None else T_init
            return compositional_split(fluid, zs, float(node.pressure_pa), T)

        outgoing = self._separator_outgoing(network_state, split_fn)
        separator_change = 0.0

        self._hydraulic_solver._initialise_pressure_field(network_state, case)
        for _outer in range(settings.max_outer_iterations):

            # Lift gas of gas-lifted wells: internal sources at their wellheads.
            self._set_gas_lift_sources(network_state, fluid)

            # Propagate compositions from inlets → PipeState.zs
            # (separators also need it: their outlets carry different phases)
            if case.inlet_composition_bcs or network_state.separators:
                node_zs_outer = self._propagate_compositions(
                    network_state, case, default_zs,
                    node_zs_prev=node_zs_outer,
                    relaxation=settings.composition_relaxation,
                    component_mws=fluid.component_mws,
                    outgoing=outgoing,
                )
                separator_change = self._update_separator_outlets(
                    network_state, node_zs_outer, default_zs, split_fn)

            old_densities = [fluid.density_for_link(link) for link in network_state.components]

            lam_hist, lam_metrics, turb_hist, turb_metrics, hydraulic_converged = self._hydraulic_pass(
                network_state, fluid, progress_callback)
            wells_started = self._start_wells(network_state, fluid)
            if network_state.separators:
                # Compositions consistent with the new flow field for the
                # energy solve: links downstream of a separator outlet that
                # carried no flow before this pass still had the default.
                node_zs_outer = self._propagate_compositions(
                    network_state, case, default_zs, node_zs_prev=node_zs_outer,
                    relaxation=settings.composition_relaxation,
                    component_mws=fluid.component_mws, outgoing=outgoing)
            all_lam_hist.extend(lam_hist)
            all_lam_metrics.extend(lam_metrics)
            all_turb_hist.extend(turb_hist)
            all_turb_metrics.extend(turb_metrics)
            outer_boundaries.append(len(all_turb_metrics))
            if turb_metrics:
                outer_turb_final.append(turb_metrics[-1])

            new_node_temps, pipe_mean_temps = solve_energy_system(
                network_state, fluid, self.convection_scheme, T_ref=T_init,
                pipe_enthalpy_terms=settings.pipe_enthalpy_terms,
            )

            max_delta_t = 0.0
            alpha = settings.temperature_relaxation
            for nid, T_new in new_node_temps.items():
                _tc = network_state.nodes[nid].temperature_c
                T_old = _tc if _tc is not None else T_init
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
            separators_settled = separator_change < SEPARATOR_FRACTION_TOLERANCE
            if temperature_converged and density_converged and separators_settled and not wells_started:
                break

        # Final synchronous pass — always runs so that the reported flow field
        # and temperature field come from the same solve with the final compositions.
        self._hydraulic_solver._initialise_pressure_field(network_state, case)
        lam_hist, lam_metrics, turb_hist, turb_metrics, hydraulic_converged = self._hydraulic_pass(
            network_state, fluid, progress_callback)
        all_lam_hist.extend(lam_hist)
        all_lam_metrics.extend(lam_metrics)
        all_turb_hist.extend(turb_hist)
        all_turb_metrics.extend(turb_metrics)
        outer_boundaries.append(len(all_turb_metrics))
        if turb_metrics:
            outer_turb_final.append(turb_metrics[-1])
        final_node_temps, final_pipe_mean_temps = solve_energy_system(
            network_state, fluid, self.convection_scheme, T_ref=T_init,
            pipe_enthalpy_terms=settings.pipe_enthalpy_terms,
        )
        for nid, T_new in final_node_temps.items():
            if not network_state.nodes[nid].is_thermal_inlet:
                network_state.nodes[nid].temperature_c = T_new
        self._update_component_temperatures(network_state, T_init, final_pipe_mean_temps)

        # ── build result ──────────────────────────────────────────────────────
        node_pressures = {
            nid: float(network_state.nodes[nid].pressure_pa)  # type: ignore[arg-type]
            for nid in sorted(network_state.nodes)
        }
        node_temperatures = {
            nid: float(network_state.nodes[nid].temperature_c if network_state.nodes[nid].temperature_c is not None else T_init)  # type: ignore[arg-type]
            for nid in sorted(network_state.nodes)
        }
        component_flows = []
        for link in network_state.components:
            density = fluid.density_for_link(link)
            t_in  = node_temperatures.get(link.start_node.node_id)
            t_out = node_temperatures.get(link.end_node.node_id)
            vf = fluid.vapor_fraction_for_link(link)
            fw = fluid.free_water_fraction_for_link(link)
            liq_frac, gas_zs, liquid_zs = fluid.flash_phase_detail_for_link(link)
            component_flows.append(
                ComponentFlowResult(
                    label=self._hydraulic_solver._component_label(link),
                    mass_flow_kg_per_s=float(link.mass_flow_kg_per_s),
                    volumetric_flow_m3_per_h=float(
                        3600.0 * link.mass_flow_kg_per_s / density
                    ),
                    temperature_in_c=t_in,
                    temperature_out_c=t_out,
                    zs=tuple(getattr(link, "zs", ())),
                    vapor_fraction=vf,
                    free_water_fraction=fw,
                    liquid_fraction=liq_frac,
                    gas_phase_zs=gas_zs,
                    liquid_phase_zs=liquid_zs,
                    multiphase=self._multiphase_report(link, fluid),
                    **self._hydrate_report(link, fluid, node_temperatures),
                )
            )

        # Derive per-node compositions (molar-flow-weighted mix of incoming links).
        _incoming: dict[int, list] = {}
        for link in network_state.components:
            _zs = getattr(link, "zs", ())
            if not _zs:
                continue
            _mdot = link.mass_flow_kg_per_s
            _nid = link.end_node.node_id if _mdot >= 0.0 else link.start_node.node_id
            _mw = sum(z * m for z, m in zip(_zs, fluid.component_mws))
            _incoming.setdefault(_nid, []).append((abs(_mdot) / _mw if _mw > 0 else abs(_mdot), tuple(_zs)))
        node_compositions: dict[int, tuple] = {}
        for _nid, _arr in _incoming.items():
            _tot = sum(m for m, _ in _arr)
            if _tot <= 0.0:
                continue
            _nc = len(_arr[0][1])
            node_compositions[_nid] = tuple(
                sum(m * zs[i] for m, zs in _arr) / _tot for i in range(_nc)
            )

        well_states = self._well_states(network_state, fluid)
        wells_settled = not any(ws.state == "unstable" for ws in well_states)
        separator_results = self._separator_results(
            network_state, frozenset(bc.node_id for bc in case.inlet_composition_bcs))
        return SolveResult(
            case_name=case.name,
            converged=(hydraulic_converged and density_converged and temperature_converged
                       and wells_settled and separator_change < SEPARATOR_FRACTION_TOLERANCE
                       and self._separators_settled(separator_results)),
            node_pressures_pa=node_pressures,
            node_temperatures_c=node_temperatures,
            component_flows=component_flows,
            laminar_history=all_lam_hist,
            laminar_metrics=all_lam_metrics,
            turbulent_history=all_turb_hist,
            turbulent_metrics=all_turb_metrics,
            temperature_history=temperature_history,
            density_rel_change_history=density_history,
            outer_turbulent_final_metrics=tuple(outer_turb_final),
            outer_iteration_boundaries=tuple(outer_boundaries),
            global_balance=self._compute_global_balance(network_state),
            global_energy_balance=self._compute_global_energy_balance(network_state, fluid),
            mass_balance_history=mass_balance_history,
            energy_balance_history=energy_balance_history,
            component_names=fluid.component_names,
            component_mws=fluid.component_mws,
            node_compositions=node_compositions,
            well_states=well_states,
            separators=separator_results,
        )
