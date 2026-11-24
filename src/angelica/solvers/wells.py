"""
Well handling shared by the black-oil and compositional network solvers.

A ``Well`` (IPR + tubing) is bistable between its static wellhead pressure
and the maximum of its deliverability curve: dead, or producing on the
stable branch.  The solvers treat this as an active set:

- a producing well dies when its wellhead reaches the curve's maximum —
  only during the turbulent hydraulic solve, since the laminar warm-start
  pressures are not physical;
- a dead well starts again, between outer iterations, when its wellhead
  falls below the static pressure (at most ``MAX_WELL_STARTS`` times per
  solve);
- a well left dead inside the band is retried forced to produce; producing
  wins when that solve converges, otherwise the well is dead;
- a dead well whose wellhead is below its static pressure has no steady
  state ("unstable") and the solve is not converged.

A solver using the mixin implements ``_solve_once(case, progress_callback,
latch_disabled)`` and builds its result with ``well_states``.
"""
from __future__ import annotations

from typing import Any

from angelica.core.results import SolveResult, WellOperatingState

MAX_WELL_STARTS = 3   # restarts of one well per solve before it is reported unstable


class WellNetworkMixin:
    """Bistable-well logic (see the module docstring)."""

    _hydraulic_solver: Any

    def solve(self, case, progress_callback=None) -> SolveResult:
        """Solve the network.

        Wells (IPR + tubing) are bistable between their static wellhead
        pressure and the maximum of their deliverability curve.  The solve
        follows each well's hysteresis; a well that ends dead inside that band
        is then given the chance to produce (re-solve with it forced to
        produce when it can).  Producing wins when it is consistent; otherwise
        the well is dead — it could not start against the network's
        pressure — and the solution stands.
        """
        result = self._solve_once(case, progress_callback)   # type: ignore[attr-defined]
        bistable_dead = frozenset(ws.component_id for ws in result.well_states
                                  if ws.state == "dead_minimum_rate")
        if not bistable_dead:
            return result
        retry = self._solve_once(case, progress_callback, latch_disabled=bistable_dead)  # type: ignore[attr-defined]
        producing = {ws.component_id for ws in retry.well_states if ws.state == "producing"}
        if retry.converged and bistable_dead <= producing:
            return retry
        return result

    @staticmethod
    def _disable_well_latch(network_state, latch_disabled) -> None:
        for link in network_state.components:
            if getattr(link.component, "is_well_model", False) and link.component.component_id in latch_disabled:
                link.well_latch_enabled = False

    def _hydraulic_pass(self, network_state, fluid, progress_callback):
        """Laminar warm start, then the turbulent solve, during which
        producing wells may die at the maximum of their curve."""
        hs = self._hydraulic_solver
        model = hs.fitting_correlation
        lam_hist, lam_metrics, _ = hs._solve_laminar(
            network_state, fluid, progress_callback=progress_callback
        )
        if hasattr(model, "well_latch_active"):
            model.well_latch_active = True
        try:
            turb_hist, turb_metrics, converged = hs._solve_turbulent(
                network_state, fluid, progress_callback=progress_callback
            )
        finally:
            if hasattr(model, "well_latch_active"):
                model.well_latch_active = False
        return lam_hist, lam_metrics, turb_hist, turb_metrics, converged

    @staticmethod
    def _set_gas_lift_sources(network_state, fluid) -> None:
        """Lift gas of every gas-lifted ``Well`` enters the network at its
        wellhead node as an internal source (mass and composition)."""
        from angelica.closures.choke import resolve_link_fluid
        from angelica.closures.well_fluid import BlackOilWellFluid, lift_gas_fluid, well_fluid_for
        from angelica.closures.well_model import gas_lift_mass_flow
        from angelica.properties.black_oil import GasStream

        for node in network_state.nodes.values():
            node.source_mass_flow_kg_per_s = 0.0
            node.source_streams = ()
        for link in network_state.components:
            comp = link.component
            if not getattr(comp, "has_gas_lift", False):
                continue
            link_fluid = resolve_link_fluid(fluid, link)
            p_res = float(link.start_node.pressure_pa)
            T = link.start_node.temperature_c
            if T is None:
                T = getattr(link_fluid, "reference_temperature_c", 60.0)
            wf = well_fluid_for(link_fluid, link, comp, round(p_res, 0), round(float(T), 3))
            if wf is None:
                continue
            m_gl = gas_lift_mass_flow(comp, link_fluid, p_res, T, link)
            gas = lift_gas_fluid(wf, comp)
            stream = GasStream(gas.key[2]) if isinstance(gas, BlackOilWellFluid) else gas.zs
            node = link.end_node
            node.source_mass_flow_kg_per_s += m_gl
            node.source_streams = (*node.source_streams, (m_gl, stream))

    @staticmethod
    def _start_wells(network_state, fluid) -> bool:
        """Between outer iterations: a dead well whose wellhead pressure has
        fallen below its static (shut-in) wellhead pressure starts producing.
        Returns True if any well started.  A well restarted
        ``MAX_WELL_STARTS`` times (it dies again each time: no steady state)
        is left dead and reported as unstable."""
        from angelica.closures.choke import resolve_link_fluid
        from angelica.closures.well_fluid import well_fluid_for
        from angelica.closures.well_model import well_link_curve

        started = False
        for link in network_state.components:
            if not (getattr(link, "well_dead", False) and getattr(link, "well_latch_enabled", False)):
                continue
            link_fluid = resolve_link_fluid(fluid, link)
            if well_fluid_for(link_fluid, link) is None:
                continue
            if (link.well_starts < MAX_WELL_STARTS
                    and well_link_curve(link, link_fluid).starts_at(float(link.end_node.pressure_pa))):
                link.well_dead = False
                link.well_starts += 1
                started = True
        return started

    @staticmethod
    def _well_states(network_state, fluid) -> tuple[WellOperatingState, ...]:
        """Operating state of every ``Well`` link at the solution."""
        from angelica.closures.choke import resolve_link_fluid
        from angelica.closures.well_fluid import well_fluid_for
        from angelica.closures.well_model import well_link_curve

        states = []
        for link in network_state.components:
            if not getattr(link.component, "is_well_model", False):
                continue
            link_fluid = resolve_link_fluid(fluid, link)
            if well_fluid_for(link_fluid, link) is None:
                continue
            curve = well_link_curve(link, link_fluid)
            p_wh = float(link.end_node.pressure_pa)
            m = float(link.mass_flow_kg_per_s)
            p_max, p_static = curve.max_wellhead_pressure_pa, curve.static_wellhead_pressure_pa
            if m > 0.0:
                state = "producing"
            elif p_max <= 0.0:
                state = "dead_column"
            elif p_wh >= p_max:
                state = "dead_back_pressure"
            elif curve.starts_at(p_wh):
                state = "unstable"
            else:
                state = "dead_minimum_rate"
            states.append(WellOperatingState(
                component_id=link.component.component_id, state=state, wellhead_pressure_pa=p_wh,
                max_wellhead_pressure_pa=p_max, static_wellhead_pressure_pa=p_static,
                minimum_stable_rate_kg_per_s=curve.minimum_stable_rate_kg_per_s, mass_flow_kg_per_s=m,
                gas_lift_kg_per_s=sum(ms for ms, _s in link.end_node.source_streams)
                if getattr(link.component, "has_gas_lift", False) else 0.0))
        return tuple(states)
