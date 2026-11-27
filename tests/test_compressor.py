"""
Gas compressor (roadmap item 5.2, v1.30.0).

Validation and its standing
---------------------------
1. ANALYTICAL — ideal gas with constant cp: the isentrope is
   T2s = T1·r^(R/cp) and W_s = cp·(T2s − T1) exactly (R per unit mass).
   R is the equation of state's own (8.314 J/(mol·K)); with 8.3144626 the
   reference itself was off by 5.6e-5 (found in review).  Errors obtained
   with the default 20 steps: −3.9e-6 (ratio 1.5) to −8.8e-5 (ratio 6) in
   W_s, −0.0001 to −0.013 K in T2s.
2. EOS ENTROPY (the roadmap's isentropic-limit check) — the march uses only
   h(P, T) and v(P, T); its end point must lie on the isentrope of the
   equation of state, s(P2, T2s) = s(P1, T1), with the entropy computed
   independently:
     - Peng-Robinson methane: s = cp·ln T − (R/M)·ln P + S_dep(thermo PR)/M,
       the entropy belonging to Angelica's h = cp·T + H_dep;
     - compositional natural gas: thermo's flash entropy (the same flash
       whose H Angelica uses).
   Obtained, T2s·|Δs| / W_s: 6.2e-5 … 1.3e-4 (PR, ratio 1.5 … 6),
   4.1e-6 … 4.7e-5 (compositional).
   Tolerance for 1 and 2: 5e-4 of W_s, the accuracy asked of the march
   (with 20 steps in ln P its error is ≤ 6e-5 of W_s against a 320-step
   march up to ratio 6, measured for PR methane, the natural gas and the
   black-oil gas); the tolerance is a numerical bound, not a property of
   data.
3. Network consistency (gas, black-oil, compositional): ratio, rise or
   power met at the node pressures within 1e-4 (the outer-loop criterion);
   the reported power equals ṁ·(h(P2, T2) − h(P1, T1)) recomputed here from
   the fluid at the reported state; mass through the compressor equals the
   pipes' within 2e-3 (the hydraulic criterion, per node); energy balance
   ≤ 0.01 %;
   discharge above the isentropic one.

Black-oil has no entropy function: its compression rests on its h and v
(consistency with the energy balance only).  No published compressor case
was read; GPSA's examples are not in the repository.
"""
from __future__ import annotations

import math
import sys
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from angelica.closures.compressor import CompressorFluid, compress  # noqa: E402
from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary  # noqa: E402
from angelica.core.components import Compressor, Pipe  # noqa: E402
from angelica.properties.black_oil import BlackOilFluid  # noqa: E402
from angelica.properties.compositional_fluid import CompositionalFluid  # noqa: E402
from angelica.properties.compressible_fluid import CompressibleFluid  # noqa: E402
from angelica.properties.eos import _R as R  # noqa: E402  (the EOS's gas constant, 8.314)
from angelica.properties.eos import IdealGasEOS, PengRobinsonEOS  # noqa: E402
from angelica.solvers import SteadyBlackOilSolver, SteadyCompositionalSolver, SteadyCompressibleSolver  # noqa: E402

M_CH4, CP = 0.016043, 2220.0
NAMES = ["methane", "ethane", "propane", "n-butane", "nitrogen", "carbon dioxide"]
ZS = [0.90, 0.05, 0.02, 0.01, 0.01, 0.01]


def methane(eos) -> CompressibleFluid:
    return CompressibleFluid.from_constants(eos=eos, viscosity_pa_s=1.1e-5, specific_heat_j_per_kg_k=CP,
                                            thermal_conductivity_w_per_m_k=0.033, reference_pressure_pa=30e5,
                                            reference_temperature_c=20.0)


PR_CH4 = PengRobinsonEOS(molecular_weight_kg_per_mol=M_CH4, critical_temperature_k=190.56,
                         critical_pressure_pa=4.599e6, acentric_factor=0.011)   # methane constants (inputs)
FLUIDS = {
    "gas": methane(PR_CH4),
    "black_oil": BlackOilFluid(api_gravity=35.0, gas_gravity=0.7, gor_sc_m3_per_m3=1e6, wor_sc_m3_per_m3=0.0),
    "compositional": CompositionalFluid(NAMES, ZS),
}
SOLVERS = {"gas": SteadyCompressibleSolver, "black_oil": SteadyBlackOilSolver,
           "compositional": SteadyCompositionalSolver}
LINK = SimpleNamespace(zs=tuple(ZS), component=None)


@pytest.mark.parametrize("ratio", [1.5, 2.0, 4.0, 6.0])
def test_ideal_gas_isentrope_is_analytical(ratio):
    fl = CompressorFluid(methane(IdealGasEOS(molecular_weight_kg_per_mol=M_CH4)), LINK)
    res = compress(fl, 30e5, 20.0, ratio * 30e5, 1.0)
    T2s = 293.15 * ratio ** (R / M_CH4 / CP) - 273.15
    assert res.isentropic_work_j_per_kg == pytest.approx(CP * (T2s - 20.0), rel=5e-4)
    assert res.t_out_isentropic_c == pytest.approx(T2s, abs=5e-4 * (T2s - 20.0))


def pr_entropy(P, T_c):
    pr = PR_CH4._pr_object(P, T_c)
    s_dep = getattr(pr, "S_dep_g", None) or getattr(pr, "S_dep_l")
    return CP * math.log(T_c + 273.15) - R / M_CH4 * math.log(P) + s_dep / M_CH4


def flash_entropy(P, T_c):
    from angelica.properties.compositional_fluid import _tp_flash
    st = _tp_flash(tuple(NAMES), P, T_c, tuple(ZS), "PR")
    return st.S() / st.MW() * 1e3                     # J/(kg·K)


@pytest.mark.parametrize("ratio", [1.5, 3.0, 6.0])
@pytest.mark.parametrize("fluid, entropy", [("gas", pr_entropy), ("compositional", flash_entropy)])
def test_isentrope_against_eos_entropy(fluid, entropy, ratio):
    fl = CompressorFluid(FLUIDS[fluid], LINK)
    res = compress(fl, 30e5, 20.0, ratio * 30e5, 1.0)
    ds = entropy(ratio * 30e5, res.t_out_isentropic_c) - entropy(30e5, 20.0)
    assert abs(ds) * (res.t_out_isentropic_c + 273.15) <= 5e-4 * res.isentropic_work_j_per_kg


@pytest.mark.parametrize("fluid", ["gas", "black_oil", "compositional"])
def test_efficiencies_order_the_discharge(fluid):
    """Work: ideal isentropic < isentropic at η < polytropic at the same η.
    The same number means a less efficient machine when it is polytropic:
    in compression the isentropic efficiency of a machine is below its
    polytropic one (each small step's losses heat the gas the next step
    compresses)."""
    fl = CompressorFluid(FLUIDS[fluid], LINK)
    s = compress(fl, 30e5, 20.0, 90e5, 1.0)
    i = compress(fl, 30e5, 20.0, 90e5, 0.75, "isentropic")
    p = compress(fl, 30e5, 20.0, 90e5, 0.75, "polytropic")
    assert i.actual_work_j_per_kg == pytest.approx(s.isentropic_work_j_per_kg / 0.75, rel=1e-9)
    assert s.isentropic_work_j_per_kg < i.actual_work_j_per_kg < p.actual_work_j_per_kg
    assert s.t_out_c < i.t_out_c < p.t_out_c


def case(fluid: str, diameter: float = 0.3, **kw) -> NetworkCase:
    def pipe(cid, a, b, length):
        return Pipe(component_id=cid, start_node=a, end_node=b, diameter_m=0.3, length_m=length,
                    absolute_roughness_m=46e-6, heat_transfer_coefficient_w_per_m2k=2.0, ambient_temperature_c=10.0)
    return NetworkCase(
        name="compressor station", fluid_model=FLUIDS[fluid],
        pressure_inlets=(PressureBoundary(1, 30e5),), pressure_outlets=(PressureBoundary(4, 50e5),),
        thermal_inlets=(ThermalBoundary(1, 20.0),),
        components=(pipe("suction", 1, 2, 2000.0),
                    Compressor(component_id="k", start_node=2, end_node=3, diameter_m=diameter, **kw),
                    pipe("discharge", 3, 4, 20000.0)),
    )


@lru_cache(maxsize=None)
def run(fluid: str, spec: str, value: float, eff_type: str = "isentropic", diameter: float = 0.3):
    return SOLVERS[fluid]().solve(case(fluid, diameter=diameter, **{spec: value, "efficiency_type": eff_type}))


def check_spec(c, spec, value):
    """Specification met at the node pressures (outer-loop criterion 1e-4)."""
    if spec == "pressure_ratio":
        assert c.pressure_ratio == pytest.approx(value, rel=1e-4)
    elif spec == "pressure_rise_pa":
        assert c.outlet_pressure_pa / c.inlet_pressure_pa == pytest.approx(1 + value / c.inlet_pressure_pa, rel=1e-4)
    else:
        # The criterion is on the ratio, |ln r − ln target| < 1e-4; the power
        # varies as r^~1.4 at ratio 2, so it is met within about 1.4e-4
        # (obtained: 1.2e-4).
        assert c.power_w == pytest.approx(value, rel=3e-4)


SPECS = [("pressure_ratio", 2.0), ("pressure_rise_pa", 30e5), ("power_w", 3e6)]


@pytest.mark.parametrize("fluid", ["gas", "black_oil", "compositional"])
@pytest.mark.parametrize("spec, value", SPECS)
def test_network_compressor(fluid, spec, value):
    res = run(fluid, spec, value)
    assert res.converged
    c = res.compressors[0]
    assert c.state == "running" and c.mass_flow_kg_per_s > 0.0
    check_spec(c, spec, value)
    assert c.outlet_temperature_c > c.outlet_temperature_isentropic_c > c.inlet_temperature_c
    fl = CompressorFluid(FLUIDS[fluid], LINK)
    dh = fl.h(c.outlet_pressure_pa, c.outlet_temperature_c) - fl.h(c.inlet_pressure_pa, c.inlet_temperature_c)
    assert c.power_w == pytest.approx(c.mass_flow_kg_per_s * dh, rel=1e-3)
    flows = {f.label.split(":")[-1]: f.mass_flow_kg_per_s for f in res.component_flows}
    # The hydraulic criterion is |in − out|/(in + out) < 1e-3 per node, so two
    # neighbouring flows may differ by up to 2e-3 (obtained: 1.3e-3, gas,
    # ratio 2 — the compressor link carries the rounding of a stiff law).
    assert flows["k"] == pytest.approx(flows["suction"], rel=2e-3)
    assert flows["discharge"] == pytest.approx(flows["suction"], rel=2e-3)
    assert abs(res.global_energy_balance.energy_error_pct) < 0.01


@pytest.mark.parametrize("ratio", [1.0, 1.5])
def test_too_small_a_ratio_passes_no_flow(ratio):
    """30 bar × 1.5 = 45 < 50 bar at the outlet: the compressor cannot push,
    and it does not run backwards."""
    res = run("gas", "pressure_ratio", ratio)
    assert res.converged
    assert res.compressors[0].state == "no flow"
    assert max(abs(f.mass_flow_kg_per_s) for f in res.component_flows) < 1e-6


# Out of range (gas).  Obtained: ratios 1.7–6 and powers 0.1–10 MW converge;
# ratio 10 and 30 MW (discharge ≈ 280 °C, the suction line nearly emptied)
# and powers of 1–30 kW (under 1 % of the design flow) end not converged with
# "compressor power" named — never an exception (roadmap limitation).
@pytest.mark.parametrize("spec, value", [("pressure_ratio", 1.7), ("pressure_ratio", 6.0),
                                         ("power_w", 1e5), ("power_w", 1e7),
                                         ("pressure_ratio", 10.0), ("power_w", 3e7), ("power_w", 0.0)])
def test_out_of_range(spec, value):
    res = run("gas", spec, value)
    if not res.converged:
        assert value in (10.0, 3e7)
        assert any("compressor" in f for f in res.diagnosis.failed)
        return
    c = res.compressors[0]
    if value == 0.0:
        assert c.state == "no flow"
    else:
        check_spec(c, spec, value)


@pytest.mark.parametrize("spec, value", SPECS)
@pytest.mark.parametrize("diameter", [0.1, 0.05])
def test_small_compressor_nozzle_still_meets_the_specification(spec, value, diameter):
    """Review of v1.30.0: the link's own loss (K·ρv²/2) was taken off the rise
    and the specification missed by up to 4.4e-3 at D = 0.05 m; the rise is
    now adjusted until the node pressures meet it."""
    res = run("gas", spec, value, diameter=diameter)
    assert res.converged
    check_spec(res.compressors[0], spec, value)


def gas_case(components, inlets=((1, 30e5),), outlets=((9, 50e5),)):
    return NetworkCase(name="g", fluid_model=FLUIDS["gas"],
                       pressure_inlets=tuple(PressureBoundary(n, p) for n, p in inlets),
                       pressure_outlets=tuple(PressureBoundary(n, p) for n, p in outlets),
                       thermal_inlets=tuple(ThermalBoundary(n, 20.0) for n, _ in inlets), components=components)


def gas_pipe(cid, a, b, length=2000.0):
    return Pipe(component_id=cid, start_node=a, end_node=b, diameter_m=0.3, length_m=length,
                absolute_roughness_m=46e-6, heat_transfer_coefficient_w_per_m2k=2.0, ambient_temperature_c=10.0)


def test_both_ends_fixed_cannot_meet_a_ratio():
    """30 bar → 50 bar boundaries with no pipe: the ratio is 1.667 whatever the
    machine; ratio 2 has no steady state.  v1.30.0 before review reported it
    converged (5 370 kg/s, 553 MW); now not converged, named."""
    res = SteadyCompressibleSolver().solve(gas_case(
        (Compressor(component_id="k", start_node=1, end_node=9, diameter_m=0.3, pressure_ratio=2.0),)))
    assert not res.converged
    assert any("compressor" in f for f in res.diagnosis.failed)


@pytest.mark.parametrize("where", ["inlet at a boundary", "outlet at a boundary"])
def test_compressor_attached_to_a_boundary(where):
    k = Compressor(component_id="k", start_node=1 if where.startswith("inlet") else 2,
                   end_node=2 if where.startswith("inlet") else 9, diameter_m=0.3, pressure_ratio=2.0)
    p = gas_pipe("p", 2, 9, 20000.0) if where.startswith("inlet") else gas_pipe("p", 1, 2, 20000.0)
    res = SteadyCompressibleSolver().solve(gas_case((k, p)))
    assert res.converged
    c = res.compressors[0]
    check_spec(c, "pressure_ratio", 2.0)
    gb = res.global_balance
    assert gb.mass_outlet_kg_per_s == pytest.approx(gb.mass_inlet_kg_per_s, rel=2e-3)


def test_identical_compressors_in_parallel_share_the_flow():
    comps = (gas_pipe("s", 1, 2), gas_pipe("d", 3, 9, 20000.0),
             Compressor(component_id="k1", start_node=2, end_node=3, diameter_m=0.3, power_w=1.5e6),
             Compressor(component_id="k2", start_node=2, end_node=3, diameter_m=0.3, power_w=1.5e6))
    res = SteadyCompressibleSolver().solve(gas_case(comps))
    assert res.converged
    a, b = res.compressors
    assert a.mass_flow_kg_per_s == pytest.approx(b.mass_flow_kg_per_s, rel=1e-3)
    check_spec(a, "power_w", 1.5e6)


def test_flow_without_a_pressure_rise_is_labelled():
    """Power spec with max_pressure_ratio 1: flow passes, no compression
    (review: it was labelled "no flow")."""
    comps = (gas_pipe("s", 1, 2), gas_pipe("d", 3, 9, 20000.0),
             Compressor(component_id="k", start_node=2, end_node=3, diameter_m=0.3, power_w=1e6,
                        max_pressure_ratio=1.0))
    res = SteadyCompressibleSolver().solve(gas_case(comps, outlets=((9, 20e5),)))
    c = res.compressors[0]
    assert c.mass_flow_kg_per_s > 0.0 and abs(c.power_w) < 1.0
    assert c.state in ("no pressure rise", "at max pressure ratio")


def test_liquid_feed_is_reported():
    """Black-oil with GOR 10: the stream is mostly liquid at suction."""
    fluid = BlackOilFluid(api_gravity=35.0, gas_gravity=0.7, gor_sc_m3_per_m3=10.0, wor_sc_m3_per_m3=0.0)
    c = case("black_oil", pressure_ratio=2.0)
    from dataclasses import replace
    res = SteadyBlackOilSolver().solve(replace(c, fluid_model=fluid))
    assert any("not all gas at suction" in w for w in res.warnings), res.warnings


def test_compressors_in_the_exports(tmp_path):
    from angelica.io.reporting import export_solve_result_csv
    res = run("gas", "pressure_ratio", 2.0)
    export_solve_result_csv(res, tmp_path / "r.csv")
    assert "Compressors" in (tmp_path / "r.csv").read_text(encoding="utf-8")


def test_gui_compressor_round_trip(tmp_path):
    import json

    from angelica.gui.io import (
        build_network_case_from_scene,
        build_solver_from_scene,
        load_scene_from_file,
        save_scene_to_file,
    )
    data = json.loads((ROOT / "src/angelica/tutorials/steady_compressible/06_compressor_heating/"
                       "06_compressor_heating.gui.json").read_text(encoding="utf-8"))
    for link in data["links"]:
        for comp in link["components"]:
            if comp["component_type"] == "pump":
                comp["component_type"] = "compressor"
                comp["properties"] = {"diameter_m": "0.3", "specification": "pressure_ratio", "value": "2.0",
                                      "efficiency": "0.75", "efficiency_type": "polytropic",
                                      "max_pressure_ratio": "10"}
    (tmp_path / "in.gui.json").write_text(json.dumps(data), encoding="utf-8")
    scene = load_scene_from_file(tmp_path / "in.gui.json")
    save_scene_to_file(scene, tmp_path / "s.gui.json")
    scene = load_scene_from_file(tmp_path / "s.gui.json")
    res = build_solver_from_scene(scene).solve(build_network_case_from_scene(scene))
    assert res.converged
    check_spec(res.compressors[0], "pressure_ratio", 2.0)


def test_polytropic_in_a_network():
    i, p = run("gas", "pressure_ratio", 2.0), run("gas", "pressure_ratio", 2.0, "polytropic")
    assert p.converged
    assert p.compressors[0].pressure_ratio == pytest.approx(2.0, rel=1e-4)
    assert p.compressors[0].outlet_temperature_c > i.compressors[0].outlet_temperature_c


def test_invalid_inputs_and_liquids():
    kw = {"component_id": "k", "start_node": 1, "end_node": 2, "diameter_m": 0.3}
    for bad in ({}, {"pressure_ratio": 2.0, "power_w": 1e6}, {"pressure_ratio": 0.5}, {"pressure_rise_pa": -1.0},
                {"power_w": math.inf}, {"pressure_ratio": 2.0, "efficiency": 0.0},
                {"pressure_ratio": 2.0, "efficiency": 1.2}, {"pressure_ratio": 2.0, "efficiency_type": "adiabatic"},
                {"pressure_ratio": 2.0, "max_pressure_ratio": 0.5}, {"pressure_ratio": 2.0, "diameter_m": 0.0}):
        with pytest.raises(ValueError):
            Compressor(**{**kw, **bad})
    from angelica.properties.single_component import SingleComponentFluid
    from angelica.solvers import SteadyIsothermalIncompressibleSolver
    water = NetworkCase(name="w", fluid_model=SingleComponentFluid(density_kg_per_m3=998.0, viscosity_pa_s=1e-3),
                        pressure_inlets=(PressureBoundary(1, 2e5),), pressure_outlets=(PressureBoundary(2, 3e5),),
                        components=(Compressor(**kw, pressure_ratio=2.0),))
    with pytest.raises(ValueError, match="Pump for liquids"):
        SteadyIsothermalIncompressibleSolver().solve(water)


def test_isothermal_solve_names_the_compressor():
    from angelica.solvers.steady_compressible import CompressibleSolverSettings
    res = SteadyCompressibleSolver(compressible_settings=CompressibleSolverSettings(isothermal_temperature_c=20.0)) \
        .solve(case("gas", pressure_ratio=2.0))
    assert res.converged
    assert any("pump/compressor" in w for w in res.warnings)


def test_flow_specified_inlet_with_a_ratio():
    """The inlet node's pressure is free (FlowBoundary): the ratio fixes it.
    Review of v1.30.0: an exception (the pressure floor did not hold on a
    flow boundary, −72 kPa), then a stalled bracket; now converged."""
    from angelica.core.case import FlowBoundary
    c = NetworkCase(name="f", fluid_model=FLUIDS["gas"], pressure_inlets=(),
                    pressure_outlets=(PressureBoundary(9, 50e5),), flow_inlets=(FlowBoundary(1, 20.0),),
                    thermal_inlets=(ThermalBoundary(1, 20.0),),
                    components=(gas_pipe("a", 1, 2), Compressor(component_id="k", start_node=2, end_node=3,
                                                                diameter_m=0.3, pressure_ratio=2.0),
                                gas_pipe("c", 3, 9, 20000.0)))
    res = SteadyCompressibleSolver().solve(c)
    assert res.converged
    check_spec(res.compressors[0], "pressure_ratio", 2.0)
    assert res.compressors[0].mass_flow_kg_per_s == pytest.approx(20.0, rel=1e-3)
