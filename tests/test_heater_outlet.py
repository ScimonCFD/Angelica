"""
Heater/cooler set by its outlet temperature (roadmap item 5.3, v1.29.0).

Validation data and their standing
----------------------------------
Water enthalpy — READ IN SOURCE: NIST Chemistry WebBook, Thermophysical
Properties of Fluid Systems, water, isobar 0.5 MPa (retrieved 2026-11-27):
h(10 °C) = 42.5077, h(20 °C) = 84.3824, h(60 °C) = 251.583,
h(80 °C) = 335.373 kJ/kg.  Angelica's liquid water has a constant
cp = 4182 J/(kg·K), so this checks that the heater duty is the stream's
enthalpy change with that model.  Tolerance 0.51 %: the largest deviation of
NIST's cp from 4182 between 5 and 80 °C (4.20333 J/(g·K) at 5 °C, same
table) — a property of the data, not of these results.  Errors obtained:
+0.048 % (heating 20 → 60 °C), −0.044 % (cooling 80 → 10 °C).

Everything else is internal consistency: the duty equals ṁ·(h_out − h_in)
with the fluid's own enthalpy, the outlet equals the set point, the energy
balance closes, and the out-of-range sweeps end converged or diagnosed.
"""
from __future__ import annotations

import importlib.util
import math
import sys
from dataclasses import replace
from functools import lru_cache
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT / "src", ROOT / "tests"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from _two_loop_networks import solver_for, with_heater  # noqa: E402

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary  # noqa: E402
from angelica.core.components import HeatSource, Pipe  # noqa: E402
from angelica.properties.water_liquid import build_water_thermal_fluid  # noqa: E402
from angelica.solvers import SteadyNonIsothermalIncompressibleSolver  # noqa: E402

TUTORIALS = ROOT / "src" / "angelica" / "tutorials"
THERMAL_MODES = ["liquid_thermal", "gas", "black_oil", "compositional"]
NIST_H_KJ_PER_KG = {10.0: 42.5077, 20.0: 84.3824, 60.0: 251.583, 80.0: 335.373}   # 0.5 MPa
NIST_TOL = 0.0051


def water_line(t_in: float, t_set: float) -> NetworkCase:
    """Adiabatic water line at about 0.5 MPa with a heater/cooler."""
    def pipe(cid, a, b):
        return Pipe(component_id=cid, start_node=a, end_node=b, diameter_m=0.1, length_m=50.0,
                    absolute_roughness_m=46e-6, heat_transfer_coefficient_w_per_m2k=0.0)
    return NetworkCase(
        name="water heater", fluid_model=build_water_thermal_fluid(),
        pressure_inlets=(PressureBoundary(1, 5.1e5),), pressure_outlets=(PressureBoundary(4, 4.9e5),),
        thermal_inlets=(ThermalBoundary(1, t_in),),
        components=(pipe("in", 1, 2),
                    HeatSource(component_id="h", start_node=2, end_node=3, diameter_m=0.1,
                               pressure_drop_mode="fixed", pressure_drop_pa=0.0, outlet_temperature_c=t_set),
                    pipe("out", 3, 4)),
    )


@pytest.mark.parametrize("t_in, t_set", [(20.0, 60.0), (80.0, 10.0)])
def test_water_duty_against_nist(t_in, t_set):
    res = SteadyNonIsothermalIncompressibleSolver().solve(water_line(t_in, t_set))
    assert res.converged
    h = res.heaters[0]
    assert h.state == "at outlet temperature"
    assert h.outlet_temperature_c == pytest.approx(t_set, abs=1e-9)
    nist = h.mass_flow_kg_per_s * (NIST_H_KJ_PER_KG[t_set] - NIST_H_KJ_PER_KG[t_in]) * 1e3
    assert abs(h.duty_w / nist - 1.0) < NIST_TOL


@lru_cache(maxsize=None)
def heated(mode: str, **kw):
    case = with_heater(mode, **kw)
    return case, solver_for(mode).solve(case)


def enthalpy(case, res, mode, node, T):
    """Specific enthalpy (J/kg) of the export-line fluid at a node, or None for
    the liquid (cp model)."""
    from types import SimpleNamespace
    fluid = case.fluid_model
    P = res.node_pressures_pa[node]
    if mode == "liquid_thermal":
        return None
    if mode == "gas":
        return fluid.enthalpy_j_per_kg(P, T)
    if mode == "black_oil":
        return fluid.enthalpy_for_link(None, P, T)
    return fluid.enthalpy_for_link(SimpleNamespace(zs=res.component_flows[-1].zs), P, T)


# Outlet temperature set point, every thermal mode.  The heater's own outlet
# is the set point by construction (its exit node is held there); what is
# tested is the node after it: it mixes the heater's stream with the
# outflows of a hydraulic solution converged to 1e-3 in mass, so it can
# miss the set point by about 1e-3 × the temperature in °C (review: 0.04 K
# for a transparent black-oil heater).  Obtained here: ≤ 0.0040 K
# (tolerance 0.01 K, the outer loops' temperature criterion, which this
# network meets — not a general bound); duty = ṁ·Δh
# within 7.4e-6 (tolerance 1e-3: 0.01 K of inlet temperature against the
# 7.8–56.5 K changes here); energy balance ≤ 1.9e-4 % (roadmap figure 0.01 %).
@pytest.mark.parametrize("mode", THERMAL_MODES)
@pytest.mark.parametrize("t_set", [40.0, 5.0])
def test_outlet_temperature_set_point(mode, t_set):
    case, res = heated(mode, outlet_temperature_c=t_set)
    assert res.converged
    h = res.heaters[0]
    assert h.state == "at outlet temperature"
    assert h.outlet_temperature_c == pytest.approx(t_set, abs=1e-9)
    assert res.node_temperatures_c[8] == pytest.approx(t_set, abs=0.01)
    assert abs(res.global_energy_balance.energy_error_pct) < 0.01
    T_in = res.node_temperatures_c[5]
    h_in, h_out = enthalpy(case, res, mode, 5, T_in), enthalpy(case, res, mode, 8, t_set)
    if h_in is None:
        expected = h.mass_flow_kg_per_s * 4182.0 * (t_set - T_in)
    else:
        expected = h.mass_flow_kg_per_s * (h_out - h_in)
    assert h.duty_w == pytest.approx(expected, rel=1e-3)
    assert (h.duty_w > 0) == (t_set > T_in)          # heater or cooler as needed


@pytest.mark.parametrize("mode", THERMAL_MODES)
def test_duty_limited(mode):
    """A 200 kW limit below the duty the 40 °C set point needs (0.91 MW in
    black-oil, 1.76 MW gas, 2.04 MW compositional, 5.0 MW water): the heater
    works at ±200 kW and the outlet lies between inlet and set point."""
    case, res = heated(mode, outlet_temperature_c=40.0, max_duty_w=2e5)
    assert res.converged
    h = res.heaters[0]
    assert h.state == "duty limited"
    assert abs(h.duty_w) == pytest.approx(2e5)
    lo, hi = sorted((h.inlet_temperature_c, 40.0))
    assert lo < h.outlet_temperature_c < hi
    assert abs(res.global_energy_balance.energy_error_pct) < 0.01


def test_fixed_duty_in_a_gas_uses_the_real_enthalpy():
    """City gate (compositional tutorial 11), 500 kW: before v1.29.0 the
    heater outlet came from a constant cp and the energy balance was off by
    0.39 %; now h(P_out, T_out) = h_in + Q/ṁ and it closes (4.2e-4 %)."""
    spec = importlib.util.spec_from_file_location(
        "t11", TUTORIALS / "steady_compositional" / "11_city_gate_regulator" / "run.py")
    t11 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(t11)
    res = t11.solve(heater_w=5e5)
    assert res.converged
    assert abs(res.global_energy_balance.energy_error_pct) < 0.01
    assert res.heaters[0].state == "fixed duty" and res.heaters[0].duty_w == 5e5


# Out of range.  Obtained: every set point converged except black-oil below
# its property floor region (−20 °C; also −5 … −10 °C off-suite), where the
# cold crude's extrapolated viscosity (thousands of mPa·s) turns the export
# line laminar and the outer loop cycles — the diagnosis names temperature
# and density (and the laminar/turbulent switch from −6 °C).
OUT_OF_RANGE = {
    "liquid_thermal": (1.0, 99.0, 150.0),
    "gas": (-60.0, 0.0, 300.0),
    "black_oil": (-20.0, 0.0, 250.0),
    "compositional": (-40.0, 0.0, 200.0),
}


@pytest.mark.parametrize("mode", THERMAL_MODES)
def test_set_points_out_of_range(mode):
    for t_set in OUT_OF_RANGE[mode]:
        case = with_heater(mode, outlet_temperature_c=t_set)
        res = solver_for(mode).solve(case)
        h = res.heaters[0]
        if mode == "black_oil" and t_set < 0.0:
            assert not res.converged and {"temperature", "density"} <= set(res.diagnosis.failed)
            continue
        assert res.converged, (mode, t_set)
        assert h.outlet_temperature_c == pytest.approx(t_set, abs=1e-9)


@pytest.mark.parametrize("mode", THERMAL_MODES)
def test_zero_duty_limit_reversed_flow_and_rest(mode):
    _, res = heated(mode, outlet_temperature_c=40.0, max_duty_w=0.0)
    assert res.converged and res.heaters[0].duty_w == 0.0 and res.heaters[0].state == "duty limited"
    case = with_heater(mode, outlet_temperature_c=40.0)
    p_in = max(b.pressure_pa for b in case.pressure_inlets)
    rev = replace(case, pressure_outlets=(PressureBoundary(node_id=6, pressure_pa=1.1 * p_in),))
    res = solver_for(mode).solve(rev)
    assert res.converged and res.heaters[0].mass_flow_kg_per_s < 0.0
    assert res.heaters[0].outlet_temperature_c == pytest.approx(40.0, abs=1e-9)
    # At rest a square-law device (heater, fitting) does not settle: round-off
    # pressure differences of the stopping tolerance drive flow through it
    # (roadmap limitation).  Required here: no exception, and a stated reason.
    p_out = case.pressure_outlets[0].pressure_pa
    rest = replace(case, pressure_inlets=tuple(replace(b, pressure_pa=p_out) for b in case.pressure_inlets))
    res = solver_for(mode).solve(rest)
    assert res.converged or res.diagnosis is not None


def test_gas_network_at_rest_with_a_fitting_does_not_raise():
    """Fixed in v1.29.0: round-off flows starved a junction (inflow 1/47 of
    its outflow) and its mixing equation gave −310 °C, then the EOS raised."""
    from _two_loop_networks import build_case

    from angelica.core.components import Fitting
    case = build_case("gas")
    comps = []
    for p in case.components:
        if p.component_id == "export":
            comps += [Fitting(component_id="valve", start_node=5, end_node=8, diameter_m=p.diameter_m,
                              loss_coefficient=5.0), replace(p, start_node=8)]
        else:
            comps.append(p)
    p_out = case.pressure_outlets[0].pressure_pa
    rest = replace(case, components=tuple(comps),
                   pressure_inlets=tuple(replace(b, pressure_pa=p_out) for b in case.pressure_inlets))
    res = solver_for("gas").solve(rest)
    assert res.converged or res.diagnosis is not None
    # between the coldest and warmest boundary or ambient temperature (15–50 °C)
    assert all(10.0 <= t <= 55.0 for t in res.node_temperatures_c.values()), res.node_temperatures_c


def test_transparent_heater_in_a_water_line_converges():
    """Fixed in v1.29.0: a heater with no pressure drop in a water network
    carried 246 kg/s against the line's 64 (its coupling was capped)."""
    case = with_heater("liquid_thermal", power_w=0.0)
    case = replace(case, components=tuple(
        replace(c, pressure_drop_mode="fixed", pressure_drop_pa=0.0) if c.component_id == "heater" else c
        for c in case.components))
    res = solver_for("liquid_thermal").solve(case)
    assert res.converged
    flows = {f.label.split(":")[-1]: f.mass_flow_kg_per_s for f in res.component_flows}
    assert flows["heater"] == pytest.approx(flows["export"], rel=1e-3)


def test_invalid_inputs():
    kw = {"component_id": "h", "start_node": 1, "end_node": 2, "diameter_m": 0.1}
    for bad in ({"outlet_temperature_c": math.nan}, {"outlet_temperature_c": -300.0},
                {"max_duty_w": -1.0}, {"max_duty_w": math.nan}, {"power_w": math.inf}):
        with pytest.raises(ValueError):
            HeatSource(**kw, **bad)


def test_isothermal_solve_reports_the_set_point_as_ignored():
    from angelica.solvers import SteadyCompressibleSolver
    from angelica.solvers.steady_compressible import CompressibleSolverSettings
    res = SteadyCompressibleSolver(compressible_settings=CompressibleSolverSettings(isothermal_temperature_c=20.0)) \
        .solve(with_heater("gas", outlet_temperature_c=40.0))
    assert res.converged
    assert any("heater/cooler duty ignored" in w for w in res.warnings)
    assert res.heaters[0].duty_w == 0.0 and res.heaters[0].state.startswith("isothermal")


def test_gui_heater_properties_reach_the_solver(tmp_path):
    from angelica.gui.io import (
        build_network_case_from_scene,
        build_solver_from_scene,
        load_scene_from_file,
        save_scene_to_file,
    )
    scene = load_scene_from_file(TUTORIALS / "steady_black_oil" / "16_heater_station" / "16_heater_station.gui.json")
    link = scene.links[0]
    heater = next(c for c in link.components if c.component_type == "heat_source")
    scene.update_link_component_properties(link.link_id, heater.component_id,
                                           {**heater.properties, "max_duty_w": "2000000"})
    save_scene_to_file(scene, tmp_path / "s.gui.json")
    scene = load_scene_from_file(tmp_path / "s.gui.json")
    case = build_network_case_from_scene(scene)
    comp = next(c for c in case.components if isinstance(c, HeatSource))
    assert comp.outlet_temperature_c == 80.0 and comp.max_duty_w == 2e6
    res = build_solver_from_scene(scene).solve(case)
    assert res.converged and res.heaters[0].state == "duty limited"


def test_heaters_in_the_exports(tmp_path):
    from angelica.io.reporting import export_solve_result_csv, export_solve_result_workbook
    _, res = heated("gas", outlet_temperature_c=40.0)
    export_solve_result_csv(res, tmp_path / "r.csv")
    text = (tmp_path / "r.csv").read_text(encoding="utf-8")
    assert "Heaters / coolers" in text and "at outlet temperature" in text
    pytest.importorskip("openpyxl")
    export_solve_result_workbook(res, tmp_path / "r.xlsx")
    from openpyxl import load_workbook
    assert "Heaters" in load_workbook(tmp_path / "r.xlsx").sheetnames


@pytest.mark.parametrize("heater_kw", [{}, {"power_w": 1e4}, {"outlet_temperature_c": 60.0}])
def test_heater_after_a_dead_well(heater_kw):
    """Black-oil tutorial 13 without lift gas: the well is dead.  A heater
    between its wellhead and the flowline carries no flow either (found in
    review: with a 0.01 m/s coupling floor it carried a phantom 0.025 kg/s
    and the solve did not converge; v1.28.0 converged)."""
    from angelica.closures import BeggsBrillCorrelation
    from angelica.solvers import SteadyBlackOilSolver
    spec = importlib.util.spec_from_file_location("t13", TUTORIALS / "steady_black_oil" / "13_gas_lift" / "run.py")
    t13 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(t13)
    case = t13.build_case(0.0)
    comps = []
    for c in case.components:
        if c.component_id == "flowline":
            comps += [HeatSource(component_id="htr", start_node=2, end_node=20, diameter_m=c.diameter_m,
                                 pressure_drop_mode="rated", pressure_drop_pa=0.2e5, rated_mass_flow_kg_per_s=10.0,
                                 **heater_kw), replace(c, start_node=20)]
        else:
            comps.append(c)
    res = SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation()).solve(
        replace(case, components=tuple(comps)))
    assert res.converged
    assert res.well_states[0].state != "producing"
    assert res.heaters[0].state == "no flow" and res.heaters[0].duty_w == 0.0
