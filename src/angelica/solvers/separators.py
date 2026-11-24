"""
Separator handling shared by the black-oil and compositional network solvers.

A separator is a node with one or more ``SeparatorOutlet`` links leaving it
(see ``angelica.core.components.SeparatorOutlet`` and
``angelica.closures.separator``).  Each outer iteration:

1. composition propagation mixes the streams arriving at the node (the
   feed) and gives each link leaving it the composition of its phases — an
   outlet its liquid, the other links the gas;
2. ``_update_separator_outlets`` flashes the feed at the node's pressure
   and temperature and stores on every outlet the mass fraction of the feed
   it carries;
3. in the hydraulic solve an outlet's flow is that fraction of the current
   feed (``ChokeFlowModel``), and the gas flow follows from the node's mass
   balance.

The solve is converged only when the fractions have stopped changing.
"""
from __future__ import annotations

from angelica.closures.separator import (
    PHASE_SETS,
    SeparatorSplit,
    gas_outlet_phases,
    separator_feed_mass_flow,
    validate_outlet_phases,
)
from angelica.core.results import SeparatorOutletResult, SeparatorResult

SEPARATOR_FRACTION_TOLERANCE = 1e-5   # outer convergence of the outlet mass fractions


class SeparatorNetworkMixin:
    """Separator logic (see the module docstring)."""

    @staticmethod
    def _setup_separators(network_state) -> None:
        """Group ``SeparatorOutlet`` links by separator node and check the layout."""
        groups: dict[int, list] = {}
        for link in network_state.components:
            if getattr(link.component, "is_separator_outlet", False):
                groups.setdefault(link.start_node.node_id, []).append(link)
        network_state.separators = {}
        for nid, outlets in groups.items():
            node = network_state.nodes[nid]
            validate_outlet_phases(nid, [o.component.phase for o in outlets])
            if node.is_inlet:
                raise ValueError(f"Separator at node {nid}: the separator node cannot be an inlet boundary.")
            others = tuple(lk for lk in network_state.components
                           if (lk.start_node is node or lk.end_node is node)
                           and not any(lk is o for o in outlets))
            if not others:
                raise ValueError(f"Separator at node {nid}: no link brings fluid to it.")
            if not node.is_pressure_boundary and len(others) < 2:
                raise ValueError(f"Separator at node {nid}: the gas has no outlet — connect a gas line "
                                 "to the node or make it a pressure boundary.")
            for o in outlets:
                if o.end_node.node_id in groups:
                    raise ValueError(f"SeparatorOutlet {o.component.component_id!r}: its outlet node "
                                     "cannot itself be a separator; put a link between them.")
                o.separator_feed_links = others
            network_state.separators[nid] = tuple(outlets)

    @staticmethod
    def _separator_outgoing(network_state, split_fn):
        """Composition carried by a link leaving a separator node, for the
        composition propagation: ``fn(link, up_id, node_comp, default)``
        returns the outlet stream's composition, or None for a link that does
        not leave a separator."""
        seps = network_state.separators
        if not seps:
            return None

        def outgoing(link, up_id, node_comp, default):
            outlets = seps.get(up_id)
            if outlets is None:
                return None
            split = split_fn(network_state.nodes[up_id], node_comp.get(up_id, default))
            if any(link is o for o in outlets):
                return split.stream(PHASE_SETS[link.component.phase])
            return split.stream(gas_outlet_phases([o.component.phase for o in outlets]))

        return outgoing

    @staticmethod
    def _update_separator_outlets(network_state, node_comp, default, split_fn) -> float:
        """Store each outlet's mass fraction of the feed; returns the largest change."""
        change = 0.0
        splits: dict[int, SeparatorSplit] = {}
        for nid, outlets in network_state.separators.items():
            split = split_fn(network_state.nodes[nid], node_comp.get(nid, default))
            splits[nid] = split
            for o in outlets:
                f = split.fraction(PHASE_SETS[o.component.phase])
                change = max(change, abs(f - o.separator_fraction))
                o.separator_fraction = f
        network_state.separator_splits = splits
        return change

    @staticmethod
    def _separator_results(network_state) -> tuple[SeparatorResult, ...]:
        out = []
        splits = getattr(network_state, "separator_splits", {})
        for nid, outlets in network_state.separators.items():
            split = splits.get(nid)
            if split is None:
                continue
            node = network_state.nodes[nid]
            feed = separator_feed_mass_flow(outlets[0])
            day = 86_400.0
            out.append(SeparatorResult(
                node_id=nid, pressure_pa=float(node.pressure_pa),
                temperature_c=float(node.temperature_c) if node.temperature_c is not None else float("nan"),
                feed_mass_flow_kg_per_s=feed,
                gas_mass_flow_kg_per_s=feed * split.fractions["gas"],
                oil_mass_flow_kg_per_s=feed * split.fractions["oil"],
                water_mass_flow_kg_per_s=feed * split.fractions["water"],
                gas_rate_sm3_per_day=feed * split.std_per_kg["gas"] * day,
                oil_rate_sm3_per_day=feed * split.std_per_kg["oil"] * day,
                water_rate_m3_per_day=feed * split.std_per_kg["water"] * day,
                outlets=tuple(SeparatorOutletResult(
                    component_id=o.component.component_id, phase=o.component.phase,
                    mass_flow_kg_per_s=float(o.mass_flow_kg_per_s),
                    valve_pressure_drop_pa=float(node.pressure_pa) - float(o.end_node.pressure_pa))
                    for o in outlets),
            ))
        return tuple(out)
