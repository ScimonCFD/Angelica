"""Round 1-2-5 axis ticks used by the GUI plots."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.gui.plot_axes import nice_ticks, tick_label, tick_values  # noqa: E402


@pytest.mark.parametrize("lo, hi, expected", [
    (-202.0, 69.0, (-300.0, 100.0, 100.0)),
    (0.0, 79.7, (0.0, 80.0, 20.0)),
    (0.0, 1.37, (0.0, 1.4, 0.2)),
    (12.3, 12.9, (12.3, 12.9, 0.1)),
])
def test_nice_ticks(lo, hi, expected):
    got = nice_ticks(lo, hi, target=5)
    assert got == pytest.approx(expected)


def test_ticks_are_round_multiples_and_cover_range():
    for lo, hi in ((-187.4, 53.2), (0.0, 148.6), (-3.3e-3, 7.7e-3), (1234.0, 98765.0)):
        a, b, step = nice_ticks(lo, hi)
        mantissa = step / 10 ** int(f"{step:e}".split("e")[1])
        assert mantissa == pytest.approx(1.0) or mantissa == pytest.approx(2.0) or mantissa == pytest.approx(5.0)
        assert a <= lo and b >= hi
        ticks = tick_values(a, b, step)
        assert 3 <= len(ticks) - 1 <= 12
        assert ticks[0] == pytest.approx(a) and ticks[-1] == pytest.approx(b)


def test_degenerate_range_and_labels():
    a, b, step = nice_ticks(5.0, 5.0)
    assert a < 5.0 < b
    assert tick_label(40.0, 20.0) == "40"
    assert tick_label(0.4, 0.2) == "0.4"
    assert tick_label(-0.0, 0.05) == "-0.00" or tick_label(0.0, 0.05) == "0.00"
