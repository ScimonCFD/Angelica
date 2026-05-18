from __future__ import annotations


def max_abs_value(values) -> float:
    return max(abs(float(value)) for value in values)
