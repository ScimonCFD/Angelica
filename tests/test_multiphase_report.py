"""
Per-pipe two-phase report ``ComponentFlowResult.multiphase`` (v1.16.0):
holdup, flow pattern, superficial velocities, liquid inventory, and the
elevation / friction split of the pressure gradient.
"""
from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures import BeggsBrillCorrelation  # noqa: E402
from angelica.core.case import InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary  # noqa: E402
from angelica.core.components import Fitting, Pipe  # noqa: E402
from angelica.io.reporting import export_solve_result_csv, multiphase_rows  # noqa: E402
from angelica.properties.black_oil import BlackOilFluid  # noqa: E402
from angelica.solvers import SteadyBlackOilSolver, SteadyCompositionalSolver  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "tutorial_compositional_07",
    SRC_ROOT / "angelica/tutorials/steady_compositional/07_hilly_condensate_beggs_brill/run.py",
)
hilly = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(hilly)


@pytest.fixture(scope="module")
def compositional_runs():
    hom = SteadyCompositionalSolver().solve(hilly.build_case())
    bb = SteadyCompositionalSolver(turbulent_pipe_correlation=BeggsBrillCorrelation()).solve(hilly.build_case())
    return hom, bb


def _section_drop(res, i):
    return res.node_pressures_pa[i + 1] - res.node_pressures_pa[i + 2]


def test_gradients_add_up_to_the_section_pressure_drop(compositional_runs):
    for res in compositional_runs:
        for i, ((_, L, _), cf) in enumerate(zip(hilly.SECTIONS, res.component_flows)):
            mp = cf.multiphase
            total = (mp.elevation_gradient_pa_per_m + mp.friction_gradient_pa_per_m) * L
            assert total == pytest.approx(_section_drop(res, i), rel=3e-3)


def test_homogeneous_report_has_no_slip(compositional_runs):
    hom, _ = compositional_runs
    for cf in hom.component_flows:
        mp = cf.multiphase
        assert mp.flow_pattern == "no-slip"
        assert mp.liquid_holdup == mp.no_slip_liquid_fraction


def test_beggs_brill_report_kinematics(compositional_runs):
    _, bb = compositional_runs
    area = math.pi * hilly.D ** 2 / 4
    for (_, L, dz), cf in zip(hilly.SECTIONS, bb.component_flows):
        mp = cf.multiphase
        assert mp.flow_pattern in ("segregated", "transition", "intermittent", "distributed")
        assert mp.liquid_holdup > mp.no_slip_liquid_fraction        # slip at low velocity
        assert mp.superficial_liquid_velocity_m_per_s + mp.superficial_gas_velocity_m_per_s == \
            pytest.approx(mp.mixture_velocity_m_per_s, rel=1e-12)
        assert mp.superficial_liquid_velocity_m_per_s == \
            pytest.approx(mp.no_slip_liquid_fraction * mp.mixture_velocity_m_per_s, rel=1e-12)
        assert mp.liquid_volume_m3 == pytest.approx(mp.liquid_holdup * area * L, rel=1e-12)
        assert (mp.elevation_gradient_pa_per_m > 0) == (dz > 0)


def test_csv_export_and_rows(compositional_runs, tmp_path):
    _, bb = compositional_runs
    rows = multiphase_rows(bb)
    assert len(rows) == len(hilly.SECTIONS) + 1 and rows[-1][0] == "Total liquid in pipes"
    assert rows[-1][7] == pytest.approx(sum(cf.multiphase.liquid_volume_m3 for cf in bb.component_flows), abs=1e-3)
    out = tmp_path / "r.csv"
    export_solve_result_csv(bb, str(out))
    text = out.read_text(encoding="utf-8")
    assert "=== MULTIPHASE FLOW ===" in text and "H_L holdup (-)" in text


def test_black_oil_report_and_reverse_flow_sign():
    """Black-oil pipes get a report; the elevation gradient follows the flow
    direction (a pipe drawn uphill but flowing downhill reports a negative
    elevation gradient). Fittings get none."""
    fluid = BlackOilFluid(api_gravity=35.0, gas_gravity=0.7, gor_sc_m3_per_m3=100.0, wor_sc_m3_per_m3=0.2)
    case = NetworkCase(
        name="bo-report", fluid_model=fluid,
        pressure_inlets=(PressureBoundary(node_id=3, pressure_pa=40e5),),
        pressure_outlets=(PressureBoundary(node_id=1, pressure_pa=20e5),),
        thermal_inlets=(ThermalBoundary(node_id=3, temperature_c=60.0, bc_type="fixed_temperature"),),
        inlet_fluid_bcs=(InletFluidBC(node_id=3, api_gravity=35.0, gas_gravity=0.7,
                                      gor_sc_m3_per_m3=100.0, wor_sc_m3_per_m3=0.2),),
        components=(Pipe(start_node=1, end_node=2, diameter_m=0.1, length_m=2000.0,
                         height_change_m=50.0, absolute_roughness_m=4.6e-5),
                    Fitting(start_node=2, end_node=3, diameter_m=0.1, loss_coefficient=1.0)),
    )
    res = SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation()).solve(case)
    assert res.converged
    pipe, fitting = res.component_flows
    assert pipe.mass_flow_kg_per_s < 0.0                     # flows from node 2 to node 1: downhill
    assert pipe.multiphase is not None and fitting.multiphase is None
    assert pipe.multiphase.elevation_gradient_pa_per_m < 0.0
    assert 0.0 < pipe.multiphase.no_slip_liquid_fraction < 1.0
