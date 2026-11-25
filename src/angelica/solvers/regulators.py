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

A part of the network whose only connections to the rest are controlling
regulators — the line between a back-pressure regulator and a reducing
regulator downstream of it — has no pressure of its own: both regulators
impose its flows (the supply at the first one's set point in, the demand at
the second one's out), and the pressure correction for it is singular
(v1.27.0 drove it to 10¹⁰ bar).  Both cannot control: if more is imposed in
than out the island's pressure rises until the regulators feeding it are
wide open; if less, it falls until the ones drawing from it are.  After
every state check (and at the start, when nothing is known yet: all of
them) the regulators on that side of such an island are set open.

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
        # No flows yet to tell which side of a floating island gives way:
        # open every regulator bounding one (the state checks then let back
        # into control the ones that can).
        for _island, incoming, outgoing, _other in RegulatorNetworkMixin._floating_islands(network_state):
            for r in incoming + outgoing:
                r.regulator_state = "open"
                controlled_node(r).regulated = False

    @staticmethod
    def _unanchored_groups(network_state):
        """Groups of nodes with no fixed pressure: connected to a pressure
        boundary or a regulated node by no link whose flow follows the
        pressures (controlling and closed regulators and separator outlets
        impose their flows).  Returns (groups, imposed links)."""
        parent = {nid: nid for nid in network_state.nodes}

        def find(n):
            while parent[n] != n:
                parent[n] = parent[parent[n]]
                n = parent[n]
            return n

        imposed = []
        for lk in network_state.components:
            comp = lk.component
            if getattr(comp, "is_pressure_regulator", False):
                if lk.regulator_state != "open":
                    imposed.append(lk)
                    continue
            elif getattr(comp, "is_separator_outlet", False):
                imposed.append(lk)
                continue
            a, b = find(lk.start_node.node_id), find(lk.end_node.node_id)
            if a != b:
                parent[a] = b
        fixed = {find(nid) for nid, n in network_state.nodes.items()
                 if n.is_pressure_boundary or n.regulated}
        groups: dict[int, set] = {}
        for nid in network_state.nodes:
            root = find(nid)
            if root not in fixed:
                groups.setdefault(root, set()).add(nid)
        return list(groups.values()), imposed

    @staticmethod
    def _floating_islands(network_state):
        """Unanchored groups bounded by at least one controlling regulator.

        Yields (nodes, incoming, outgoing, imposed): the controlling
        regulators feeding the island (back-pressure regulators: their end
        node is in it), those drawing from it (reducing regulators: their
        start node is in it), and the other imposed-flow links crossing its
        boundary (separator outlets)."""
        groups, imposed = RegulatorNetworkMixin._unanchored_groups(network_state)
        out = []
        for nodes in groups:
            incoming, outgoing, other = [], [], []
            for lk in imposed:
                s_in, e_in = lk.start_node.node_id in nodes, lk.end_node.node_id in nodes
                if s_in == e_in:
                    continue
                if getattr(lk.component, "is_pressure_regulator", False):
                    if lk.regulator_state == "controlling":
                        (incoming if e_in else outgoing).append(lk)
                else:
                    other.append(lk)
            if incoming or outgoing:
                out.append((nodes, incoming, outgoing, other))
        return out

    @staticmethod
    def _flow_at_set_point(r, pressures: dict) -> float:
        """The flow a regulator would carry with its controlled node at the
        set point: the node's balance at the pass's flows, moved from the
        pass's pressure (``pressures``: before any state change reset it to
        the set point) to the set point along the secant couplings of its
        other links."""
        node = controlled_node(r)
        m = regulator_mass_flow(r)
        S = sum(abs(float(lk.coupling)) for lk in r.regulator_links)
        dp = r.component.set_pressure_pa - pressures.get(node.node_id, float(node.pressure_pa))
        # net inflow into the node falls by S per Pa its pressure rises
        m += -S * dp if r.component.mode == "upstream" else S * dp
        return max(m, 0.0)

    def _resolve_floating_islands(self, network_state, pressures: dict) -> bool:
        """Open the regulators on the giving side of each floating island
        (see the module docstring); True if a state changed."""
        changed = False
        for nodes, incoming, outgoing, other in self._floating_islands(network_state):
            m_in = sum(self._flow_at_set_point(r, pressures) for r in incoming)
            m_out = sum(self._flow_at_set_point(r, pressures) for r in outgoing)
            for lk in other:
                q = float(lk.mass_flow_kg_per_s)
                if lk.end_node.node_id in nodes:
                    m_in += q
                else:
                    m_out += q
            m_in += sum(float(network_state.nodes[n].source_mass_flow_kg_per_s) for n in nodes)
            # More imposed in than out: the island's pressure rises and the
            # regulators feeding it open; less: it falls and those drawing
            # from it open.  A side with no regulator cannot give way: then
            # the other one does (the imposed flows themselves are
            # inconsistent, e.g. a level-control valve into a reducer).
            side = incoming if m_in >= m_out else outgoing
            for r in side or incoming or outgoing:
                self._set_regulator_state(r, "open")
                changed = True
        return changed

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
        pressures = {nid: float(n.pressure_pa) for nid, n in network_state.nodes.items()}
        floating = set().union(*self._unanchored_groups(network_state)[0])
        to_control = []
        capped_open = not hydraulic_converged and any(
            r.regulator_state == "controlling" and r.regulator_cap_hits and not r.regulator_can_control
            and r.regulator_switches < MAX_SWITCHES for r in network_state.regulators)
        for r in network_state.regulators:
            if r.regulator_switches >= MAX_SWITCHES:
                continue
            comp = r.component
            p_set = comp.set_pressure_pa
            p_up, p_down = float(r.start_node.pressure_pa), float(r.end_node.pressure_pa)
            state = r.regulator_state
            new = state
            cap_hits, r.regulator_cap_hits = r.regulator_cap_hits, 0
            if state == "controlling" and capped_open:
                # Some regulator in this unsettled pass was held where it
                # cannot be (limited by its wide-open capacity): it opens, and
                # the flows of the others are its consequence — a reducer
                # behind one held at 68 bar with 10 bar upstream saw zero flow
                # and closed, cutting the line off its supply.  They are
                # checked again after the next pass.
                if cap_hits and not r.regulator_can_control:
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
                    # opening holds it — with the other regulators as they
                    # are.  Applied one at a time below.
                    to_control.append(r)
                    continue
            elif state == "closed":
                # The node that decides is the controlled one.  If nothing
                # fixes its pressure (a line between two closed regulators)
                # that pressure is arbitrary and says nothing: stay closed.
                if controlled_node(r).node_id in floating:
                    pass
                elif (p_down < p_set) if comp.mode == "downstream" else (p_up > p_set):
                    new = "controlling"
            if new != state:
                self._set_regulator_state(r, new)
                changed = True
        if to_control:
            # One open regulator back into control per check, the one with
            # the highest inlet pressure (upstream in a series).  Each one's
            # feasibility was judged with the others as they are: two in
            # series switching together each assumed the other wide open (a
            # reducer to 10 bar upstream of one to 20 bar both went back to
            # control, and the second was then held where it cannot be).
            r = max(to_control, key=lambda lk: float(lk.start_node.pressure_pa))
            self._set_regulator_state(r, "controlling")
            # control is feasible: a pass that then does not settle is
            # numerical, not a reason to open (until another state changes)
            r.regulator_can_control = True
            changed = True
        if self._resolve_floating_islands(network_state, pressures):
            changed = True
        if changed:
            # "can control" was judged with the other regulators' states of
            # then; once one of them changes it no longer holds.
            for lk in network_state.regulators:
                if lk.regulator_state != "controlling" or not (to_control and lk is r):
                    lk.regulator_can_control = False
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
    def _regulator_warnings(network_state) -> tuple[str, ...]:
        """Nodes cut off by closed regulators: no flow and no pressure of
        their own, so the pressure reported for them is arbitrary."""
        if not network_state.regulators:
            return ()
        groups, _imposed = RegulatorNetworkMixin._unanchored_groups(network_state)
        out = []
        for nodes in groups:
            names = ", ".join(str(n) for n in sorted(nodes))
            out.append(f"node(s) {names}: cut off by closed regulators, no flow and no pressure of "
                       f"their own (the pressure reported there is arbitrary)")
        return tuple(out)

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
