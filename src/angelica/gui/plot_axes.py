"""Round axis ticks (1-2-5 steps) for the GUI's canvas plots."""
from __future__ import annotations

import math


def nice_ticks(lo: float, hi: float, target: int = 6) -> tuple[float, float, float]:
    """Expand ``[lo, hi]`` to multiples of a 1-2-5 step giving about ``target`` intervals.

    Returns ``(axis_min, axis_max, step)``.
    """
    if not (math.isfinite(lo) and math.isfinite(hi)):
        raise ValueError("axis limits must be finite")
    if hi < lo:
        lo, hi = hi, lo
    span = hi - lo
    if span <= 0.0:
        span = abs(hi) * 0.1 or 1.0
        lo, hi = lo - span / 2.0, hi + span / 2.0
    raw = span / max(target, 1)
    mag = 10.0 ** math.floor(math.log10(raw))

    def _bounds(step: float) -> tuple[float, float]:
        return math.floor(lo / step + 1e-9) * step, math.ceil(hi / step - 1e-9) * step

    def _score(step: float) -> tuple[float, float]:
        a, b = _bounds(step)
        n = round((b - a) / step)
        return abs(n - target), -n          # closest to target; ties -> more ticks

    step = min((m * mag for m in (1.0, 2.0, 5.0, 10.0)), key=_score)
    axis_min, axis_max = _bounds(step)
    return axis_min, axis_max, step


def tick_values(axis_min: float, axis_max: float, step: float) -> list[float]:
    """Tick positions from ``axis_min`` to ``axis_max`` inclusive."""
    n = int(round((axis_max - axis_min) / step))
    decimals = max(0, -math.floor(math.log10(step)) + 1)
    return [round(axis_min + i * step, decimals) + 0.0 for i in range(n + 1)]


def tick_label(value: float, step: float) -> str:
    """Format a tick with just enough decimals for the step."""
    decimals = max(0, -math.floor(math.log10(step) + 1e-9))
    return f"{value:.{decimals}f}"
