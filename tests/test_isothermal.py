"""
Isothermal solves of the gas, black-oil and compositional solvers
(``isothermal_temperature_c`` in their settings; roadmap goal part 1b).

Internal consistency only: an isothermal solve is a modelling choice, with
nothing to compare against but itself.  Checked: every node and the well
tubing at the given temperature; heater duties and other thermal
boundaries reported as ignored; out of range (−40 to 200 °C, network at
rest, reversed flow) every case converged or not converged with a stated
reason — never an exception.
"""
from __future__ import annotations

import importlib.util
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT / "src", ROOT / "tests"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from _two_loop_networks import build_case, solver_for  # noqa: E402

from angelica.closures.well_model import well_tubing_state  # noqa: E402
from angelica.core.case import PressureBoundary  # noqa: E402

TUTORIALS = ROOT / "src" / "angelica" / "tutorials"


def tutorial(path: str):
    spec = importlib.util.spec_from_file_location(path.replace("/", "_"), TUTORIALS / path / "run.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_black_oil_well_tubing_is_isothermal():
    """Gas-lifted black-oil well (tutorial 13): isothermal tubing ends at the
    temperature it starts at; with the Ramey profile it cools."""
    t13 = tutorial("steady_black_oil/13_gas_lift")
    from angelica.properties.black_oil import BlackOilFluid
    well, fluid = t13.make_well(60_000.0), BlackOilFluid(**t13.CRUDE)
    iso = well_tubing_state(well, fluid, t13.P_RES, 60.0, 3.0, link_state=SimpleNamespace(isothermal=True))
    ramey = well_tubing_state(well, fluid, t13.P_RES, 60.0, 3.0)
    assert iso.wellhead_temperature_c == pytest.approx(60.0, abs=1e-9)
    assert iso.bottomhole_temperature_c == pytest.approx(60.0, abs=1e-9)
    assert ramey.wellhead_temperature_c < 59.0
    assert iso.wellhead_pressure_pa != pytest.approx(ramey.wellhead_pressure_pa, rel=1e-4)


def test_compositional_well_tubing_is_isothermal():
    t09 = tutorial("steady_compositional/09_gas_condensate_wells")
    well, fluid = t09.make_well(1), t09.fluid_of(1)
    p_r = t09.WELLS[1][2] * t09.BAR
    iso = well_tubing_state(well, fluid, p_r, 80.0, 5.0, link_state=SimpleNamespace(isothermal=True, zs=()))
    assert iso.wellhead_temperature_c == pytest.approx(80.0, abs=1e-9)


@pytest.mark.parametrize("T", [20.0, 60.0])
def test_isothermal_network_with_wells_and_separator(T):
    """Black-oil tutorial 14 (wells, separator) solved isothermally."""
    from angelica.closures import BeggsBrillCorrelation
    from angelica.solvers import SteadyBlackOilSolver
    from angelica.solvers.steady_black_oil import BlackOilSolverSettings
    t14 = tutorial("steady_black_oil/14_separator")
    res = SteadyBlackOilSolver(black_oil_settings=BlackOilSolverSettings(isothermal_temperature_c=T),
                               turbulent_pipe_correlation=BeggsBrillCorrelation()).solve(t14.build_case())
    assert res.converged
    assert all(t == T for t in res.node_temperatures_c.values())
    assert res.separators[0].temperature_c == T
    assert all(ws.state == "producing" for ws in res.well_states)


def test_heater_duty_is_reported_as_ignored():
    """Compositional tutorial 11 (city gate) with its heater on."""
    from angelica.solvers import SteadyCompositionalSolver
    from angelica.solvers.steady_compositional import CompositionalSolverSettings
    t11 = tutorial("steady_compositional/11_city_gate_regulator")
    res = SteadyCompositionalSolver(
        compositional_settings=CompositionalSolverSettings(isothermal_temperature_c=20.0),
    ).solve(t11.build_case(heater_w=2e5))
    assert res.converged
    assert any("heater/cooler duty ignored" in w for w in res.warnings)


# Out of range.  Results obtained (reviewed in v1.28.0): every case converged
# except black-oil at 5, 5.5 and 6 °C, where a pipe sits at Re ≈ 2300 and the
# friction factor jumps from 64/Re to Colebrook (roadmap limitation,
# maintainer's choice); the diagnosis names the pipe and both regimes.
NOT_CONVERGING = {("black_oil", 5.0), ("black_oil", 5.5), ("black_oil", 6.0)}


@pytest.mark.parametrize("mode, temps", [
    ("gas", (-150.0, -40.0, -6.0, 0.0, 100.0, 200.0, 1000.0)),
    ("black_oil", (-150.0, -40.0, -6.0, 0.0, 4.5, 5.0, 5.5, 6.0, 7.0, 100.0, 200.0)),
    ("compositional", (-60.0, -20.0, 0.0, 150.0, 250.0)),
])
def test_isothermal_out_of_range(mode, temps):
    for T in temps:
        res = solver_for(mode, T).solve(build_case(mode))
        assert all(t == T for t in res.node_temperatures_c.values()), T
        if (mode, T) in NOT_CONVERGING:
            assert not res.converged, T
            switches = res.diagnosis.pattern_switches
            assert switches and set(switches[0].patterns) == {"laminar", "turbulent"}, (T, res.diagnosis)
        else:
            assert res.converged, (mode, T)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -273.15, -300.0])
def test_isothermal_temperature_must_be_physical(bad):
    from angelica.core.case import ThermalBoundary
    from angelica.solvers import BlackOilSolverSettings, CompositionalSolverSettings, CompressibleSolverSettings
    for settings in (CompressibleSolverSettings, BlackOilSolverSettings, CompositionalSolverSettings):
        with pytest.raises(ValueError):
            settings(isothermal_temperature_c=bad)
    with pytest.raises(ValueError):
        ThermalBoundary(node_id=1, temperature_c=bad)
    ThermalBoundary(node_id=1, temperature_c=bad, bc_type="zero_gradient")   # temperature unused there


def test_devices_whose_temperature_change_is_dropped_are_named():
    """Compressible tutorial 06 (compressor heats 20 → 69 °C) and 05 (valve
    Joule-Thomson 40 → 31 °C): isothermal, the warnings name them."""
    from angelica.solvers import CompressibleSolverSettings, SteadyCompressibleSolver
    solver = SteadyCompressibleSolver(compressible_settings=CompressibleSolverSettings(isothermal_temperature_c=20.0))
    from angelica.gui.io import build_network_case_from_scene, load_scene_from_file
    for folder, word in (("06_compressor_heating", "pump/compressor"), ("05_valve_jt_effect", "Joule-Thomson")):
        scene_file = next((TUTORIALS / "steady_compressible" / folder).glob("*.gui.json"))
        res = solver.solve(build_network_case_from_scene(load_scene_from_file(scene_file)))
        assert any(word in w for w in res.warnings), (folder, res.warnings)


@pytest.mark.parametrize("mode", ["gas", "black_oil", "compositional"])
def test_isothermal_at_rest_and_reversed(mode):
    case = build_case(mode)
    p_out = case.pressure_outlets[0].pressure_pa
    rest = replace(case, pressure_inlets=tuple(replace(b, pressure_pa=p_out) for b in case.pressure_inlets))
    res = solver_for(mode, 40.0).solve(rest)
    assert res.converged
    assert max(abs(f.mass_flow_kg_per_s) for f in res.component_flows) < 1e-6
    p_in = max(b.pressure_pa for b in case.pressure_inlets)
    rev = replace(case, pressure_outlets=(PressureBoundary(node_id=6, pressure_pa=1.1 * p_in),))
    res = solver_for(mode, 40.0).solve(rev)
    assert res.converged
    assert res.component_flows[-1].mass_flow_kg_per_s < 0.0        # flows back in from the outlet


@pytest.mark.parametrize("scene_file, T", [
    ("steady_black_oil/02_looped_gathering_network/looped_gathering_network.gui.json", 40.0),
    ("steady_compositional/03_looped_network/03_looped_network.gui.json", 40.0),
    ("steady_compressible/03_looped_gas_pipeline_heat_loss/looped_gas_pipeline_heat_loss.gui.json", 15.0),
])
def test_gui_scene_isothermal_round_trip(tmp_path, scene_file, T):
    """The GUI's "Isothermal at" setting is saved, loaded and reaches the
    solver; with it a scene needs no fixed-temperature boundary."""
    from angelica.gui.io import (
        build_network_case_from_scene,
        build_solver_from_scene,
        load_scene_from_file,
        save_scene_to_file,
    )
    scene = load_scene_from_file(TUTORIALS / scene_file)
    scene.update_solver_settings({"isothermal_temperature_c": T})
    for node in scene.nodes:
        if node.node_type in ("source", "sink"):
            node.properties["thermal_bc_type"] = "zero_gradient"
    path = tmp_path / "iso.gui.json"
    save_scene_to_file(scene, path)
    scene = load_scene_from_file(path)
    res = build_solver_from_scene(scene).solve(build_network_case_from_scene(scene))
    assert res.converged
    assert all(t == T for t in res.node_temperatures_c.values())
    scene.update_solver_settings({"isothermal_temperature_c": None})      # blank: energy solved
    with pytest.raises(ValueError, match="fixed temperature"):
        build_network_case_from_scene(scene)          # then a fixed-temperature boundary is needed
    scene = load_scene_from_file(TUTORIALS / scene_file)
    res = build_solver_from_scene(scene).solve(build_network_case_from_scene(scene))
    assert res.global_energy_balance is not None
