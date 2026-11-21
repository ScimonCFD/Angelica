"""
Regression tests for the segmented wall heat-loss discretisation and for the
tutorial import paths.

Wall heat loss
--------------
A pipe with ``n_thermal_segments = n`` has n-1 interior FV nodes; the junction
nodes carry no wall term.  The whole wall area U·π·D·L must therefore be shared
by the n-1 interior nodes.  (Before v1.10.1 each node used L/n, silently
dropping 1/n of the pipe area — 10 % less heat loss with 10 segments.)

For pure advection with upwind differencing the discrete solution is exactly

    θ_out = (T_out − T_amb)/(T_in − T_amb) = (1 + NTU/(n−1))^−(n−1)

which converges to exp(−NTU) as n → ∞.

Tutorial paths
--------------
Every tutorial inserts ``SRC_ROOT`` into ``sys.path``.  It must point at the
directory that contains the ``angelica`` package, otherwise the tutorial
silently runs against whatever Angelica version is installed.
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica import (  # noqa: E402
    NetworkCase,
    NonIsothermalSolverSettings,
    Pipe,
    PressureBoundary,
    SteadyNonIsothermalIncompressibleSolver,
    ThermalBoundary,
    ThermalFluid,
)
from angelica.closures import UpwindScheme  # noqa: E402

_CP = 2000.0
_U, _D, _L = 100.0, 0.01, 5.0
_T_IN, _T_AMB = 80.0, 20.0

_FLUID = ThermalFluid.from_constants(
    density_kg_per_m3=900.0,
    viscosity_pa_s=0.05,
    specific_heat_j_per_kg_k=_CP,
    thermal_conductivity_w_per_m_k=1e-9,   # pure advection
)


def _solve(n_segments: int):
    case = NetworkCase(
        name=f"wall-loss-n{n_segments}",
        fluid_model=_FLUID,
        components=(Pipe(
            start_node=1, end_node=2, diameter_m=_D, length_m=_L,
            absolute_roughness_m=0.0, heat_transfer_coefficient_w_per_m2k=_U,
            ambient_temperature_c=_T_AMB, n_thermal_segments=n_segments,
        ),),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=111_325.0),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=101_325.0),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=_T_IN),),
    )
    solver = SteadyNonIsothermalIncompressibleSolver(
        non_isothermal_settings=NonIsothermalSolverSettings(
            max_temperature_iterations=30, temperature_tolerance_k=1e-6,
        ),
        convection_scheme=UpwindScheme(),
    )
    result = solver.solve(case)
    assert result.converged
    mdot = result.component_flows[0].mass_flow_kg_per_s
    ntu = _U * math.pi * _D * _L / (mdot * _CP)
    theta = (result.node_temperatures_c[2] - _T_AMB) / (_T_IN - _T_AMB)
    return result, ntu, theta


@pytest.mark.parametrize("n", [3, 5, 10, 50])
def test_segmented_heat_loss_matches_discrete_upwind_solution(n):
    """θ_out equals (1 + NTU/(n−1))^−(n−1): the full wall area is applied."""
    _, ntu, theta = _solve(n)
    theta_discrete = (1.0 + ntu / (n - 1)) ** (-(n - 1))
    assert abs(theta - theta_discrete) < 1e-4, (
        f"n={n}: θ_solver={theta:.5f}, θ_discrete={theta_discrete:.5f} "
        f"(missing wall area would give {(1 + ntu * (n - 1) / n / (n - 1)) ** (-(n - 1)):.5f})"
    )


def test_segmented_heat_loss_converges_to_analytical():
    """With many segments the outlet approaches exp(−NTU) (within 0.1 K)."""
    _, ntu, theta = _solve(400)
    T_out = _T_AMB + theta * (_T_IN - _T_AMB)
    T_exact = _T_AMB + math.exp(-ntu) * (_T_IN - _T_AMB)
    assert abs(T_out - T_exact) < 0.1


@pytest.mark.parametrize("n", [5, 50])
def test_global_energy_balance_closes(n):
    """Reported wall loss matches the enthalpy drop of the stream."""
    result, _, _ = _solve(n)
    geb = result.global_energy_balance
    assert geb is not None
    assert geb.energy_error_pct < 0.1, f"energy error {geb.energy_error_pct:.3f} %"


# ── Tutorial import paths ──────────────────────────────────────────────────────

_TUTORIALS = sorted((SRC_ROOT / "angelica" / "tutorials").glob("*/*/run.py"))


def test_tutorials_found():
    assert len(_TUTORIALS) >= 40


@pytest.mark.parametrize("run_py", _TUTORIALS, ids=lambda p: f"{p.parent.parent.name}/{p.parent.name}")
def test_tutorial_src_root_points_at_package(run_py):
    """SRC_ROOT must contain angelica/ so the tutorial uses the source tree."""
    text = run_py.read_text(encoding="utf-8")
    m = re.search(r'SRC_ROOT\s*=\s*Path\(__file__\)\.resolve\(\)\.parents\[(\d+)\](\s*/\s*"src")?', text)
    assert m, "tutorial does not define SRC_ROOT"
    root = run_py.resolve().parents[int(m.group(1))]
    if m.group(2):
        root = root / "src"
    assert (root / "angelica" / "__init__.py").exists(), f"SRC_ROOT resolves to {root}"
