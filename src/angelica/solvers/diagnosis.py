"""
Diagnosis of a solve that did not converge (``SolveResult.diagnosis``).

The black-oil and compositional solvers record each pipe's flow pattern
after every outer iteration.  When the solve ends without converging, the
diagnosis names the criteria not met, the interior node with the largest
mass imbalance, and the pipes whose flow pattern changed over the last
outer iterations.  Beggs-Brill's pattern boundaries are discontinuous: a
pipe can have no solution at a flow imposed on it, which shows up as a
pattern alternating between iterations (black-oil tutorial 14 at 0.4–0.5
MPa export).
"""
from __future__ import annotations

from angelica.core.results import ConvergenceDiagnosis, PatternSwitch
from angelica.core.state import PipeState

from .steady_isothermal_incompressible import nodal_mass_imbalances

PATTERN_WINDOW = 6     # outer iterations looked at for pattern changes


class DiagnosisMixin:
    """Flow-pattern history and the end-of-solve diagnosis."""

    def _reset_flow_pattern_trace(self) -> None:
        corr = getattr(self._hydraulic_solver, "turbulent_pipe_correlation", None)   # type: ignore[attr-defined]
        if corr is not None and hasattr(corr, "pattern_trace"):
            corr.pattern_trace = {}

    def _pattern_trace(self) -> dict:
        """The correlation's last flow patterns per pipe, within the final
        hydraulic pass (a pattern can alternate between hydraulic iterations
        while each outer iteration ends in the same one)."""
        corr = getattr(self._hydraulic_solver, "turbulent_pipe_correlation", None)   # type: ignore[attr-defined]
        return getattr(corr, "pattern_trace", {}) or {}

    def _record_flow_patterns(self, network_state, fluid, history: dict) -> None:
        for link in network_state.components:
            if not isinstance(link, PipeState):
                continue
            report = self._multiphase_report(link, fluid)   # type: ignore[attr-defined]
            if report is not None:
                history.setdefault(link.component.component_id, []).append(report.flow_pattern)

    def _diagnose(self, failed: list[str], network_state, history: dict) -> ConvergenceDiagnosis | None:
        if not failed:
            return None
        imbalance = nodal_mass_imbalances(network_state)
        worst = max(imbalance, key=lambda nid: imbalance[nid]) if imbalance else None
        trace = self._pattern_trace()
        switches = []
        for link in network_state.components:
            cid = link.component.component_id
            outer = tuple(dict.fromkeys(history.get(cid, [])[-PATTERN_WINDOW:]))
            inner = tuple(dict.fromkeys(trace.get(cid, ())))
            seen = tuple(dict.fromkeys(outer + inner)) if len(outer) > 1 or len(inner) > 1 else ()
            if len(seen) > 1:
                switches.append(PatternSwitch(component_id=cid, patterns=seen))
        return ConvergenceDiagnosis(
            failed=tuple(failed),
            worst_node_id=worst,
            worst_node_imbalance_rel=imbalance.get(worst, 0.0) if worst is not None else 0.0,
            pattern_switches=tuple(switches),
        )

    @staticmethod
    def _failed_criteria(hydraulic: bool, temperature: bool, density: bool, separator_split: bool,
                         wells: bool, separator_results=()) -> list[str]:
        failed = [name for name, ok in (("hydraulics", hydraulic), ("temperature", temperature),
                                        ("density", density), ("separator split", separator_split),
                                        ("wells", wells)) if not ok]
        failed += [f"separator {s.node_id}: back-flow into a vessel with no gas phase"
                   for s in separator_results if s.backflow_without_gas]
        return failed
