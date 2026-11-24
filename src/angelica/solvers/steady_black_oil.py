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
from angelica.properties.base import FluidModel
from angelica.properties.black_oil import BlackOilComposition, mix_black_oil_streams

from .base import BaseSolver
from .steady_isothermal_incompressible import SteadyIsothermalIncompressibleSolver
from .wells import WellNetworkMixin


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
        pipe_enthalpy_terms: Joule-Thomson cooling of the gas, friction
            heating of the liquids and potential energy along pipes, with
            enthalpy-conserving junctions and isenthalpic chokes/valves
            (default True, since v1.20.0).  False drops them for comparison;
            the energy balance is then still computed with real enthalpy.
    """

    max_outer_iterations: int = 50
    density_rel_tolerance: float = 1e-4
    temperature_tolerance_k: float = 0.01
    temperature_relaxation: float = 1.0
    composition_relaxation: float = 1.0
    pipe_enthalpy_terms: bool = True


class SteadyBlackOilSolver(WellNetworkMixin, BaseSolver):
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
            fitting_correlation=ChokeFlowModel(),
        )

    def _set_correlation_fluid(self, fluid) -> None:
        """Pass the current (possibly per-pipe) fluid to fluid-aware correlations."""
        for corr in (
            self._hydraulic_solver.turbulent_pipe_correlation,
            self._hydraulic_solver.fitting_correlation,
        ):
            if hasattr(corr, "set_fluid"):
                corr.set_fluid(fluid)

    # ── Composition helpers ───────────────────────────────────────────────────

    @staticmethod
    def _propagate_compositions(
        network_state: NetworkState,
        case: NetworkCase,
        default_comp: BlackOilComposition,
        node_comp_prev: dict[int, BlackOilComposition] | None = None,
        relaxation: float = 1.0,
    ) -> tuple[dict[int, BlackOilComposition], dict[int, BlackOilComposition]]:
        """Return a per-pipe-index composition map after one propagation pass.

        Compositions are seeded at inlet nodes and propagated downstream
        following the current flow field.  At each junction the incoming
        streams are mixed conserving standard volumes
        (:func:`~angelica.properties.black_oil.mix_black_oil_streams`).  Pipes in which flow
        is ambiguous or zero inherit the default composition.

        Returns
        -------
        dict mapping index-in-network_state.components → BlackOilComposition
        """
        # ── seed inlet node compositions ──────────────────────────────────────
        node_comp: dict[int, BlackOilComposition] = {
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
        incoming: dict[int, list] = defaultdict(list)
        # Every link carries composition downstream — pipes and heat sources,
        # but also chokes, check valves, well inflows and pumps (a device with
        # no volume still delivers its upstream fluid to the next node).
        pipe_states = list(enumerate(network_state.components))

        inlet_node_ids = frozenset(ibc.node_id for ibc in case.inlet_fluid_bcs)

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
            for nid_s, node_s in network_state.nodes.items():
                for m_s, stream in node_s.source_streams:     # lift gas (GasStream)
                    incoming[nid_s].append((m_s, stream))

            # Mix at each non-inlet junction
            for nid, contributions in incoming.items():
                if nid in inlet_node_ids:
                    continue  # inlet nodes keep their prescribed composition
                if not contributions:
                    continue
                total_m = sum(m for m, _ in contributions)
                if total_m <= 1e-12:
                    continue
                mixed = mix_black_oil_streams(contributions, default_comp)
                # Apply inter-iteration relaxation toward the previous outer result
                if node_comp_prev is not None and relaxation < 1.0 and nid in node_comp_prev:
                    old = node_comp_prev[nid]
                    a = relaxation
                    mixed = BlackOilComposition(
                        api_gravity      = a * mixed.api_gravity      + (1 - a) * old.api_gravity,
                        gas_gravity      = a * mixed.gas_gravity      + (1 - a) * old.gas_gravity,
                        gor_sc_m3_per_m3 = a * mixed.gor_sc_m3_per_m3 + (1 - a) * old.gor_sc_m3_per_m3,
                        wor_sc_m3_per_m3 = a * mixed.wor_sc_m3_per_m3 + (1 - a) * old.wor_sc_m3_per_m3,
                    )

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
        pipe_comp: dict[int, BlackOilComposition] = {}
        for idx, ps in pipe_states:
            ṁ = ps.mass_flow_kg_per_s
            if ṁ >= 0.0:
                up_id = ps.start_node.node_id
            else:
                up_id = ps.end_node.node_id
            pipe_comp[idx] = node_comp.get(up_id, default_comp)


        return pipe_comp, node_comp

    # ── FluidModel proxy ──────────────────────────────────────────────────────

    # ── Mixed-composition FluidModel ──────────────────────────────────────────

    def _build_mixed_fluid(
        self,
        network_state: NetworkState,
        pipe_comps: dict[int, BlackOilComposition],
        default_comp: BlackOilComposition,
        ref_t: float,
        base_fluid=None,
    ):
        """Return a FluidModel-compatible object that dispatches per link.

        The per-link fluids keep the case fluid's PVT correlation and gas
        thermal properties (they used to fall back to Standing and the
        defaults)."""
        from angelica.properties.black_oil import BlackOilFluid

        base_correlation = getattr(base_fluid, "correlation", "STANDING")
        base_gas_cp = getattr(base_fluid, "gas_specific_heat_j_per_kg_k", 2_200.0)
        base_gas_k = getattr(base_fluid, "gas_thermal_conductivity_w_per_mk", 0.035)

        # Build a lookup: component_object_id → per-pipe fluid
        id_to_comp: dict[int, BlackOilComposition] = {}
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
            correlation             = base_correlation,
            gas_specific_heat_j_per_kg_k     = base_gas_cp,
            gas_thermal_conductivity_w_per_mk = base_gas_k,
        )

        class _PerPipeFluid:
            def __init__(self):
                self._cache: dict = {}

            def _fluid_for(self, link_state) -> BlackOilFluid:
                key = id(link_state)
                if key not in self._cache:
                    comp = id_to_comp.get(key, default_comp)
                    self._cache[key] = BlackOilFluid(
                        api_gravity             = comp.api_gravity,
                        gas_gravity             = comp.gas_gravity,
                        gor_sc_m3_per_m3        = comp.gor_sc_m3_per_m3,
                        wor_sc_m3_per_m3        = comp.wor_sc_m3_per_m3,
                        reference_pressure_pa   = ref_p,
                        reference_temperature_c = ref_t,
                        correlation             = base_correlation,
                        gas_specific_heat_j_per_kg_k     = base_gas_cp,
                        gas_thermal_conductivity_w_per_mk = base_gas_k,
                    )
                return self._cache[key]

            def fluid_for_link(self, link_state) -> BlackOilFluid:
                """Per-link BlackOilFluid (used by fluid-aware correlations)."""
                return self._fluid_for(link_state)

            def density_for_link(self, link_state) -> float:
                return self._fluid_for(link_state).density_for_link(link_state)

            def viscosity_for_link(self, link_state) -> float:
                return self._fluid_for(link_state).viscosity_for_link(link_state)

            def specific_heat_for_link(self, link_state) -> float:
                return self._fluid_for(link_state).specific_heat_for_link(link_state)

            def thermal_conductivity_for_link(self, link_state) -> float:
                return self._fluid_for(link_state).thermal_conductivity_for_link(link_state)

            def enthalpy_for_link(self, link_state, pressure_pa: float, temperature_c: float) -> float:
                return self._fluid_for(link_state).enthalpy_for_link(link_state, pressure_pa, temperature_c)

            def reference_enthalpy_for_link(self, link_state) -> float:
                return 0.0

            def specific_heat_at(self, link_state, pressure_pa: float, temperature_c: float) -> float:
                return self._fluid_for(link_state).specific_heat_at(link_state, pressure_pa, temperature_c)

            def temperature_from_enthalpy_for_link(self, link_state, h, pressure_pa, T_guess_c=20.0):
                return self._fluid_for(link_state).temperature_from_enthalpy_for_link(
                    link_state, h, pressure_pa, T_guess_c)

        return _PerPipeFluid() if id_to_comp else default_fluid

    # ── Main solve ────────────────────────────────────────────────────────────

    @staticmethod
    def _hydrate_report(link, fluid, node_temperatures) -> dict:
        """Hydrate check of a black-oil pipe with the gas-gravity (Katz) chart,
        at the end where the fluid is coldest relative to the formation curve.
        Also reports the in-situ water volume fraction as ``free_water_fraction``
        (black-oil: volume, not mole, fraction), which the GUI uses to flag
        links where hydrate can actually form."""
        from angelica.core.state import PipeState
        from angelica.properties.hydrate import katz_hydrate_temperature
        if not isinstance(link, PipeState) or abs(link.mass_flow_kg_per_s) < 1e-12:
            return {}
        f = fluid.fluid_for_link(link) if hasattr(fluid, "fluid_for_link") else fluid
        if getattr(f, "gor_sc_m3_per_m3", 0.0) <= 0.0:
            return {}
        out: dict = {}
        try:
            for node in (link.start_node, link.end_node):
                T = node_temperatures.get(node.node_id)
                if T is None or node.pressure_pa is None:
                    continue
                pt = katz_hydrate_temperature(f.gas_gravity, float(node.pressure_pa))
                if pt is None:
                    continue
                sub = pt.temperature_c - T
                if "hydrate_subcooling_k" not in out or sub > out["hydrate_subcooling_k"]:
                    out.update(hydrate_subcooling_k=sub, hydrate_temperature_c=pt.temperature_c,
                               hydrate_structure=pt.structure)
            p_s, p_e = link.start_node.pressure_pa, link.end_node.pressure_pa
            if f.wor_sc_m3_per_m3 > 0.0 and p_s is not None and p_e is not None:
                P = 0.5 * (float(p_s) + float(p_e))
                T_m = link.temperature_c if link.temperature_c is not None else T
                out["free_water_fraction"] = f.pvt(P, T_m).holdup_water
        except Exception:      # a check must never break a converged solve
            return {}
        return out

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

        self._require_fixed_temperature(case)
        network_state = build_network_state(case)
        self._disable_well_latch(network_state, latch_disabled)
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
        pipe_comps: dict[int, BlackOilComposition] | None = None
        node_comp_prev: dict[int, BlackOilComposition] | None = None
        effective_fluid: FluidModel = fluid_model  # default; updated inside loop when use_per_inlet

        self._hydraulic_solver._initialise_pressure_field(network_state, case)
        for _outer in range(settings.max_outer_iterations):

            # Lift gas of gas-lifted wells: internal sources at their wellheads.
            self._set_gas_lift_sources(network_state, effective_fluid)

            # ── build effective fluid model for this iteration ────────────────
            if use_per_inlet:
                pipe_comps, node_comp_prev = self._propagate_compositions(
                    network_state, case, default_comp,  # type: ignore[arg-type]
                    node_comp_prev=node_comp_prev,
                    relaxation=settings.composition_relaxation,
                )
                effective_fluid = self._build_mixed_fluid(
                    network_state, pipe_comps, default_comp, T_init,  # type: ignore[arg-type]
                    base_fluid=fluid_model,
                )
            else:
                effective_fluid = fluid_model

            old_densities = [
                effective_fluid.density_for_link(link)
                for link in network_state.components
            ]

            # Give fluid-aware correlations (Beggs-Brill, choke) the current fluid.
            self._set_correlation_fluid(effective_fluid)

            lam_hist, lam_metrics, turb_hist, turb_metrics, hydraulic_converged = self._hydraulic_pass(
                network_state, effective_fluid, progress_callback)
            wells_started = self._start_wells(network_state, effective_fluid)
            all_lam_hist.extend(lam_hist)
            all_lam_metrics.extend(lam_metrics)
            all_turb_hist.extend(turb_hist)
            all_turb_metrics.extend(turb_metrics)
            outer_boundaries.append(len(all_turb_metrics))
            if turb_metrics:
                outer_turb_final.append(turb_metrics[-1])

            new_node_temps, pipe_mean_temps = solve_energy_system(
                network_state, effective_fluid, self.convection_scheme, T_ref=T_init,
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

            _gb = self._compute_global_balance(network_state)
            mass_balance_history.append(_gb.mass_error_pct)
            _geb = self._compute_global_energy_balance(network_state, effective_fluid)
            energy_balance_history.append(_geb.energy_error_pct)

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
            if temperature_converged and density_converged and not wells_started:
                break

        # ── final synchronous pass ────────────────────────────────────────────
        # Always runs so that the reported flow field and temperature field
        # come from the same solve with the final compositions and PVT.
        # Pressure field is kept as-is (converged state) — no re-initialisation.
        self._set_correlation_fluid(effective_fluid)
        lam_hist, lam_metrics, turb_hist, turb_metrics, hydraulic_converged = self._hydraulic_pass(
            network_state, effective_fluid, progress_callback)
        all_lam_hist.extend(lam_hist)
        all_lam_metrics.extend(lam_metrics)
        all_turb_hist.extend(turb_hist)
        all_turb_metrics.extend(turb_metrics)
        outer_boundaries.append(len(all_turb_metrics))
        if turb_metrics:
            outer_turb_final.append(turb_metrics[-1])
        final_node_temps, final_pipe_mean_temps = solve_energy_system(
            network_state, effective_fluid, self.convection_scheme, T_ref=T_init,
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
        for idx, link in enumerate(network_state.components):
            if use_per_inlet and pipe_comps is not None:
                density = effective_fluid.density_for_link(link)
            else:
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
                    multiphase=self._multiphase_report(link, effective_fluid),
                    **self._hydrate_report(link, effective_fluid, node_temperatures),
                )
            )

        well_states = self._well_states(network_state, effective_fluid)
        wells_settled = not any(ws.state == "unstable" for ws in well_states)
        return SolveResult(
            case_name=case.name,
            converged=hydraulic_converged and density_converged and temperature_converged and wells_settled,
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
            global_energy_balance=self._compute_global_energy_balance(
                network_state, effective_fluid
            ),
            mass_balance_history=mass_balance_history,
            energy_balance_history=energy_balance_history,
            well_states=well_states,
        )
