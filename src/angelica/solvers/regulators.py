"""
Pressure regulators in the black-oil and compositional network solvers.

Each ``PressureRegulator`` is in one of three states (see
``angelica.closures.regulator``).  It starts controlling; after each outer
iteration its state is checked against the solution:

downstream (reducing) regulator, controlled node = end node:
  controlling → closed       the flow it would carry is zero or negative:
                             the downstream side is held at or above the set
                             point by something else
  controlling → open         the upstream pressure minus the wide-open
                             pressure drop at that flow is below the set point
  open        → controlling  wide open, the end node is above the set point
  closed      → controlling  closed, the end node has fallen below it

upstream (back-pressure) regulator, controlled node = start node:
  controlling → closed       the flow it would carry is zero or negative:
                             nothing arrives to hold the start node up
  controlling → open         the downstream pressure plus the wide-open
                             pressure drop is above the set point
  open        → controlling  wide open, the start node is below the set point
  closed      → controlling  closed, the start node has risen above it

The solve is converged only when no regulator changed state in the last
outer iteration.  A regulator that keeps switching (more than
``MAX_SWITCHES`` times) is left in its last state and the solve reports
converged=False, naming it.
"""
from __future__ import annotations

from angelica.closures.regulator import controlled_node, open_pressure_drop, regulator_mass_flow
from angelica.core.results import RegulatorResult

MAX_SWITCHES = 12
ZERO_FLOW_KG_PER_S = 1e-5


class RegulatorNetworkMixin:
    """Regulator logic (see the module docstring)."""

    @staticmethod
    def _setup_regulators(network_state) -> None:
        regs = [lk for lk in network_state.components
                if getattr(lk.component, "is_pressure_regulator", False)]
        held: dict[int, str] = {}
        for r in regs:
            comp = r.component
            node = controlled_node(r)
            if node.is_boundary:
                raise ValueError(f"PressureRegulator {comp.component_id!r}: the node it controls "
                                 f"({node.node_id}) is a boundary; its pressure is already fixed.")
            if node.node_id in held:
                raise ValueError(f"PressureRegulator {comp.component_id!r}: node {node.node_id} is "
                                 f"already controlled by {held[node.node_id]!r}.")
            held[node.node_id] = comp.component_id
            r.regulator_links = tuple(lk for lk in network_state.components
                                      if lk is not r and (lk.start_node is node or lk.end_node is node))
            if not r.regulator_links:
                raise ValueError(f"PressureRegulator {comp.component_id!r}: nothing else is connected "
                                 f"to the node it controls ({node.node_id}).")
            r.regulator_state = "controlling"
            r.regulator_switches = 0
            r.regulator_can_control = False
            node.regulated = True
            node.pressure_pa = comp.set_pressure_pa
        network_state.regulators = tuple(regs)

    @staticmethod
    def _set_regulator_state(r, state: str) -> None:
        node = controlled_node(r)
        r.regulator_state = state
        r.regulator_capped = False
        r.regulator_switches += 1
        node.regulated = state == "controlling"
        if node.regulated:
            node.pressure_pa = r.component.set_pressure_pa

    def _update_regulator_states(self, network_state, fluid, hydraulic_converged: bool = True) -> bool:
        """Check every regulator against the solution; True if one changed state.

        A controlling regulator whose hydraulic pass did not converge and that
        was limited by its wide-open capacity during the pass opens: the
        supply cannot meet the demand at the set point, and the pass
        oscillates between drawing the demand and being limited."""
        changed = False
        flows = [abs(float(lk.mass_flow_kg_per_s)) for lk in network_state.components]
        # zero flow: round-off of the largest flow, or below 1e-5 kg/s (under
        # 1 kg/day) — a pipe with exactly zero pressure difference still
        # reports ~1e-6 kg/s, and a network at rest has nothing larger.
        tol = max(1e-9 * max(flows, default=0.0), ZERO_FLOW_KG_PER_S)
        for r in network_state.regulators:
            if r.regulator_switches >= MAX_SWITCHES:
                continue
            comp = r.component
            p_set = comp.set_pressure_pa
            p_up, p_down = float(r.start_node.pressure_pa), float(r.end_node.pressure_pa)
            state = r.regulator_state
            new = state
            cap_hits, r.regulator_cap_hits = r.regulator_cap_hits, 0
            if state == "controlling" and not hydraulic_converged and cap_hits and not r.regulator_can_control:
                new = "open"
            elif state == "controlling":
                m = regulator_mass_flow(r)
                if m <= tol:
                    # nothing to pass (or a negative flow): closed.  A zero
                    # flow with the controlled node at the set point is the
                    # boundary between the two states; closed is the
                    # well-posed one (all flows zero left the hydraulics
                    # unsettled with the regulator controlling).
                    new = "closed"
                else:
                    dp_open = open_pressure_drop(r, m, fluid.density_for_link(r))
                    short = (p_up - dp_open < p_set) if comp.mode == "downstream" else (p_down + dp_open > p_set)
                    if short:
                        new = "open"
            elif state == "open":
                if (p_down > p_set) if comp.mode == "downstream" else (p_up < p_set):
                    # Wide open it passes the set point, so some partial
                    # opening holds it: control is feasible, and a pass that
                    # then does not settle is numerical, not a reason to open.
                    new = "controlling"
                    r.regulator_can_control = True
            elif state == "closed":
                if (p_down < p_set) if comp.mode == "downstream" else (p_up > p_set):
                    new = "controlling"
            if new != state:
                self._set_regulator_state(r, new)
                changed = True
        return changed

    @staticmethod
    def _held_temperatures(network_state, T_default: float):
        """Present node and pipe temperatures, in the form the energy solve
        returns them (pipe index counts pipes and heat sources)."""
        from angelica.core.state import HeatSourceState, PipeState
        nodes = {nid: (n.temperature_c if n.temperature_c is not None else T_default)
                 for nid, n in network_state.nodes.items()}
        pipes = {}
        idx = 0
        for link in network_state.components:
            if isinstance(link, (PipeState, HeatSourceState)):
                if link.temperature_c is not None:
                    pipes[idx] = link.temperature_c
                idx += 1
        return nodes, pipes

    @staticmethod
    def _regulators_stuck(network_state) -> list[str]:
        return [r.component.component_id for r in network_state.regulators
                if r.regulator_switches >= MAX_SWITCHES]

    @staticmethod
    def _regulator_results(network_state) -> tuple[RegulatorResult, ...]:
        out = []
        for r in network_state.regulators:
            comp = r.component
            out.append(RegulatorResult(
                component_id=comp.component_id, mode=comp.mode, set_pressure_pa=comp.set_pressure_pa,
                state=r.regulator_state,
                inlet_pressure_pa=float(r.start_node.pressure_pa),
                outlet_pressure_pa=float(r.end_node.pressure_pa),
                mass_flow_kg_per_s=float(r.mass_flow_kg_per_s),
            ))
        return tuple(out)
