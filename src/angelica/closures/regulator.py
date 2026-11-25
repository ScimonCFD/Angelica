"""
Pressure regulator (``PressureRegulator``, roadmap item 5.4).

A regulator holds the pressure of one of its nodes at a set point:

- ``mode="downstream"`` (pressure-reducing regulator): the end node;
- ``mode="upstream"`` (back-pressure regulator): the start node.

While it controls, the controlled node is a fixed-pressure node of the
pressure correction and the regulator's flow is whatever balances that
node's other links: the demand downstream of a reducing regulator, the
supply upstream of a back-pressure regulator.  The flow passes from start to
end only.  When the set point cannot be held the regulator is either wide
open (a K-factor fitting, ``loss_coefficient``) or closed (no flow); the
solver switches between the three states between outer iterations (see
``angelica.solvers.regulators``).
"""
from __future__ import annotations


def controlled_node(link_state):
    """The node whose pressure the regulator holds."""
    if link_state.component.mode == "downstream":
        return link_state.end_node
    return link_state.start_node


def regulator_mass_flow(link_state) -> float:
    """Flow (start → end) that balances the controlled node's other links
    and its internal sources."""
    node = controlled_node(link_state)
    net_in = float(node.source_mass_flow_kg_per_s)
    for link in link_state.regulator_links:
        q = float(link.mass_flow_kg_per_s)
        if link.end_node is node:
            net_in += q
        elif link.start_node is node:
            net_in -= q
    return -net_in if link_state.component.mode == "downstream" else net_in


def open_pressure_drop(link_state, mass_flow: float, density: float) -> float:
    """Pressure drop of the wide-open regulator at ``mass_flow``: K·ρ·v²/2."""
    if mass_flow <= 0.0:
        return 0.0
    v = mass_flow / (density * link_state.area_m2)
    return link_state.component.loss_coefficient * density * v * v / 2.0
