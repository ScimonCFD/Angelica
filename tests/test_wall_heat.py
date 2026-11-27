"""
Wall heat transfer: cylindrical layers, buried pipe, computed inner film
(roadmap item 6.1, v1.38.0).

Validation
----------
- Cylindrical conduction through concentric layers, R' = ln(r_o/r_i)/(2πk):
  the exact solution of the radial heat equation (ANALYTICAL); a single
  layer is checked against it to 1e-12, and the thin-wall limit against the
  old flat-wall sum (t/r → 0).
- Buried pipe, S/L = 2π/arccosh(2H/D): checked ANALYTICALLY by the method
  of images (a line source and its image make the cylinder an isotherm;
  their temperature difference gives the same S to 1e-12), and against the
  worked value of the EES heat-transfer library, function SF_2 (D = 0.5 m,
  depth 2 m, L = 5 m: 11.33 m; READ IN A SECONDARY SOURCE, citing
  Incropera & DeWitt, 5th ed., pp. 193–194): Angelica 11.347 m (+0.15 %;
  the library value matches the ln(4z/D) form of the same table entry,
  11.331 m — the arccosh form is the exact one).  Tolerance 0.2 %: the
  difference between the two forms of that table at this depth.
- Inner film: Gnielinski's correlation and its range, Nu = 3.66 laminar —
  READ IN A SECONDARY SOURCE (Wikipedia "Nusselt number", citing Incropera
  & DeWitt).  NOT validated against measurements; by consistency it agrees
  with Gnielinski's own simplified power laws (same page) within 13.1 % and
  with Dittus–Boelter within 10.4 % over Re 1e4–1e6 (guards 15 % and 12 %).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT / "src", ROOT / "tests"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from angelica.closures.wall_heat import (  # noqa: E402
    NU_LAMINAR,
    RE_LAMINAR,
    RE_TURBULENT,
    gnielinski_nusselt,
    inner_film,
    nusselt,
)
from angelica.core.components import Pipe  # noqa: E402


def pipe(**kw):
    base = {"start_node": 1, "end_node": 2, "diameter_m": 0.2, "length_m": 1000.0}
    return Pipe(**{**base, **kw})


def test_single_cylindrical_layer_is_exact():
    p = pipe(insulation_thickness_m=0.05, insulation_thermal_conductivity_w_per_mk=0.04)
    r_i, r_o = 0.1, 0.15
    exact = 0.04 / (r_i * math.log(r_o / r_i))                 # U on the bore: k/(r_i·ln(r_o/r_i))
    assert p.u_overall_w_per_m2k == pytest.approx(exact, rel=1e-12)


def test_thin_wall_limit_is_the_flat_sum():
    kw = {"inner_film_coefficient_w_per_m2k": 500.0, "wall_thickness_m": 1e-5,
          "wall_thermal_conductivity_w_per_mk": 50.0, "outer_film_coefficient_w_per_m2k": 10.0}
    flat = 1.0 / (1 / 500 + 1e-5 / 50 + 1 / 10)
    assert pipe(diameter_m=10.0, **kw).u_overall_w_per_m2k == pytest.approx(flat, rel=1e-5)


def test_buried_pipe_by_the_method_of_images():
    """Line source at depth a = √(H² − r²) and its image above the surface:
    on the cylinder (axis depth H, radius r) the distance ratio is constant,
    ρ_image/ρ_source = (H + a)/r, and S/L = 2π/ln of that ratio."""
    H, r = 1.2, 0.25
    a = math.sqrt(H * H - r * r)
    ratios = []
    for k in range(12):
        th = 2 * math.pi * k / 12
        x, y = r * math.cos(th), -H + r * math.sin(th)          # surface at y = 0
        ratios.append(math.hypot(x, y - a) / math.hypot(x, y + a))
    assert max(ratios) - min(ratios) < 1e-12
    s_images = 2 * math.pi / math.log(ratios[0])
    s_formula = 2 * math.pi / math.acosh(H / r)
    assert s_images == pytest.approx(s_formula, rel=1e-12)
    p = pipe(diameter_m=2 * r, burial_depth_m=H, soil_thermal_conductivity_w_per_mk=1.5)
    assert p.u_overall_w_per_m2k * math.pi * 2 * r == pytest.approx(1.5 * s_formula, rel=1e-12)  # U·πD = k·S/L


def test_buried_pipe_against_the_ees_example():
    s = 2 * math.pi * 5.0 / math.acosh(2 * 2.0 / 0.5)
    assert s == pytest.approx(11.33, rel=2e-3)
    p = pipe(diameter_m=0.5, length_m=5.0, burial_depth_m=2.0, soil_thermal_conductivity_w_per_mk=1.0)
    assert p.u_overall_w_per_m2k * math.pi * 0.5 * 5.0 == pytest.approx(s, rel=1e-12)


def test_layer_inputs_are_checked():
    for bad in ({"wall_thickness_m": 0.01},                                     # no conductivity
                {"burial_depth_m": 1.0},                                        # no soil k
                {"burial_depth_m": 0.05, "soil_thermal_conductivity_w_per_mk": 1.0},   # above the surface
                {"burial_depth_m": 1.0, "soil_thermal_conductivity_w_per_mk": 1.0,
                 "outer_film_coefficient_w_per_m2k": 10.0},
                {"insulation_thickness_m": -0.01, "insulation_thermal_conductivity_w_per_mk": 0.04}):
        with pytest.raises(ValueError):
            pipe(**bad)


def test_inner_film_correlation():
    assert nusselt(1000.0, 5.0) == NU_LAMINAR
    assert nusselt(RE_LAMINAR, 5.0) == NU_LAMINAR
    assert nusselt(RE_TURBULENT, 5.0) == pytest.approx(gnielinski_nusselt(RE_TURBULENT, 5.0), rel=1e-12)
    mid = nusselt(2650.0, 5.0)
    assert NU_LAMINAR < mid < gnielinski_nusselt(RE_TURBULENT, 5.0)            # continuous across
    worst_simplified = worst_db = 0.0
    for re in (1e4, 3e4, 1e5, 1e6):
        for pr in (0.7, 1.2):
            g = gnielinski_nusselt(re, pr)
            worst_simplified = max(worst_simplified, abs(g / (0.0214 * (re**0.8 - 100) * pr**0.4) - 1))
            worst_db = max(worst_db, abs(g / (0.023 * re**0.8 * pr**0.4) - 1))
        for pr in (3.0, 10.0, 100.0):
            g = gnielinski_nusselt(re, pr)
            worst_simplified = max(worst_simplified, abs(g / (0.012 * (re**0.87 - 280) * pr**0.4) - 1))
    assert worst_simplified < 0.15 and worst_db < 0.12
    film = inner_film(20.0, 0.2, 1e-3, 4180.0, 0.6)          # water, Re 1.27e5, Pr 6.97
    assert film.in_range and film.reynolds == pytest.approx(4 * 20.0 / (math.pi * 0.2 * 1e-3))
    assert film.h_w_per_m2k == pytest.approx(gnielinski_nusselt(film.reynolds, film.prandtl) * 0.6 / 0.2)
    assert not inner_film(20.0, 0.2, 1e-3, 4180.0, 1e-4).in_range             # Pr 4.2e4 > 2000


def _hot_water_line(**pipe_kw):
    from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
    from angelica.properties.thermal_fluid import ThermalFluid
    fluid = ThermalFluid.from_constants(density_kg_per_m3=980.0, viscosity_pa_s=4.5e-4,
                                        specific_heat_j_per_kg_k=4190.0, thermal_conductivity_w_per_m_k=0.66)
    return NetworkCase(name="w", fluid_model=fluid, pressure_inlets=(PressureBoundary(1, 4e5),),
                       pressure_outlets=(PressureBoundary(2, 1e5),), thermal_inlets=(ThermalBoundary(1, 80.0),),
                       components=(pipe(diameter_m=0.15, length_m=2000.0, absolute_roughness_m=46e-6,
                                        ambient_temperature_c=5.0, n_thermal_segments=20, **pipe_kw),))


def test_solver_uses_the_computed_film_and_reports_it():
    """The reported U is the layers' U with the film the correlation gives at
    the pipe's flow and properties; the global energy balance closes as for
    a fixed U (< 0.01 %, the suite's criterion)."""
    from angelica.solvers import SteadyNonIsothermalIncompressibleSolver
    kw = {"compute_inner_film": True, "wall_thickness_m": 0.007, "wall_thermal_conductivity_w_per_mk": 45.0,
          "insulation_thickness_m": 0.04, "insulation_thermal_conductivity_w_per_mk": 0.035,
          "outer_film_coefficient_w_per_m2k": 12.0}
    case = _hot_water_line(**kw)
    res = SteadyNonIsothermalIncompressibleSolver().solve(case)
    f = res.component_flows[0]
    assert res.converged and f.inner_film_coefficient_w_per_m2k is not None
    film = inner_film(f.mass_flow_kg_per_s, 0.15, 4.5e-4, 4190.0, 0.66)
    assert f.inner_film_coefficient_w_per_m2k == pytest.approx(film.h_w_per_m2k, rel=1e-9)
    assert f.heat_transfer_coefficient_w_per_m2k == pytest.approx(case.components[0].overall_u(film.h_w_per_m2k),
                                                                rel=1e-12)
    assert abs(res.global_energy_balance.energy_error_pct) < 0.01


def test_buried_line_cools_less_than_bare_in_air():
    from angelica.solvers import SteadyNonIsothermalIncompressibleSolver
    s = SteadyNonIsothermalIncompressibleSolver()
    air = s.solve(_hot_water_line(wall_thickness_m=0.007, wall_thermal_conductivity_w_per_mk=45.0,
                                  outer_film_coefficient_w_per_m2k=25.0))
    soil = s.solve(_hot_water_line(wall_thickness_m=0.007, wall_thermal_conductivity_w_per_mk=45.0,
                                   burial_depth_m=1.0, soil_thermal_conductivity_w_per_mk=1.2))
    assert soil.node_temperatures_c[2] > air.node_temperatures_c[2]
    assert soil.component_flows[0].heat_transfer_coefficient_w_per_m2k < air.component_flows[0].heat_transfer_coefficient_w_per_m2k


@pytest.mark.parametrize("mode", ["gas", "black_oil", "compositional"])
def test_computed_film_in_every_thermal_mode(mode):
    """The film is computed and reported in the gas, black-oil and
    compositional solvers too (homogeneous mixture properties)."""
    from _two_loop_networks import fluid_for, solver_for

    from angelica.core.case import FlowBoundary, NetworkCase, PressureBoundary, ThermalBoundary
    m, p_out = {"gas": (10.0, 40e5), "black_oil": (10.0, 10e5), "compositional": (10.0, 40e5)}[mode]
    case = NetworkCase(name="p", fluid_model=fluid_for(mode), pressure_inlets=(), flow_inlets=(FlowBoundary(1, m),),
                       pressure_outlets=(PressureBoundary(2, p_out),), thermal_inlets=(ThermalBoundary(1, 40.0),),
                       components=(pipe(length_m=3000.0, absolute_roughness_m=46e-6, ambient_temperature_c=5.0,
                                        compute_inner_film=True, wall_thickness_m=0.01,
                                        wall_thermal_conductivity_w_per_mk=45.0,
                                        burial_depth_m=1.5, soil_thermal_conductivity_w_per_mk=1.0),))
    res = solver_for(mode).solve(case)
    f = res.component_flows[0]
    assert res.converged
    assert f.inner_film_coefficient_w_per_m2k is not None and f.inner_film_coefficient_w_per_m2k > 0.0
    assert f.heat_transfer_coefficient_w_per_m2k == pytest.approx(
        case.components[0].overall_u(f.inner_film_coefficient_w_per_m2k), rel=1e-12)
    assert res.node_temperatures_c[2] < 40.0


def test_gui_auto_inner_film_and_burial(tmp_path):
    from angelica.gui.io import build_network_case_from_scene, load_scene_from_file
    scene = load_scene_from_file(ROOT / "src/angelica/tutorials/steady_black_oil/07_hot_oil_pipeline/"
                                 "07_hot_oil_pipeline.gui.json")
    comp = scene.links[0].components[0]
    scene.update_link_component_properties(scene.links[0].link_id, comp.component_id,
                                           {**comp.properties, "inner_film_coefficient_w_per_m2k": "auto",
                                            "burial_depth_m": "1.2", "soil_thermal_conductivity_w_per_mk": "1.4"})
    p = next(c for c in build_network_case_from_scene(scene).components if isinstance(c, Pipe))
    assert p.compute_inner_film and p.burial_depth_m == 1.2 and p.soil_thermal_conductivity_w_per_mk == 1.4


def test_buried_insulated_pipe_by_hand():
    """Burial with wall and insulation: the shape factor takes the outer
    diameter over the insulation (D_o = 0.2 + 2·(0.01 + 0.05) = 0.32 m)."""
    p = pipe(wall_thickness_m=0.01, wall_thermal_conductivity_w_per_mk=45.0, insulation_thickness_m=0.05,
             insulation_thermal_conductivity_w_per_mk=0.04, burial_depth_m=1.0, soil_thermal_conductivity_w_per_mk=1.2)
    r_i = 0.1
    R = r_i * math.log(0.11 / 0.10) / 45.0 + r_i * math.log(0.16 / 0.11) / 0.04 + r_i * math.acosh(1.0 / 0.16) / 1.2
    assert p.u_overall_w_per_m2k == pytest.approx(1.0 / R, rel=1e-12)


def test_computed_film_alone_goes_in_series_with_the_given_u():
    """"auto" with only an overall U given: 1/U = 1/h_i + 1/U_given (v1.38.0
    review: the given U was dropped and U became h_i)."""
    p = pipe(heat_transfer_coefficient_w_per_m2k=5.0, compute_inner_film=True)
    assert p.overall_u(100.0) == pytest.approx(1.0 / (1.0 / 100.0 + 1.0 / 5.0), rel=1e-12)
    assert p.overall_u() == 5.0                       # before the first film (no film yet)
    with pytest.raises(ValueError, match="not both"):
        pipe(compute_inner_film=True, inner_film_coefficient_w_per_m2k=100.0)
    with pytest.raises(ValueError, match="burial_depth_m > 0"):
        pipe(soil_thermal_conductivity_w_per_mk=1.0)
    with pytest.raises(ValueError, match="Pipe p7"):
        pipe(component_id="p7", wall_thickness_m=0.01)
    with pytest.raises(ValueError):
        inner_film(float("nan"), 0.2, 1e-3, 4180.0, 0.6)


def test_film_range_flags_and_warning():
    """Transition (Re 2300–3000) and laminar flow still developing over the
    pipe length are flagged, and the solve says so in its warnings."""
    film_t = inner_film(0.8, 0.2, 2e-3, 2000.0, 0.13)                             # Re 2546: transition
    assert 2300.0 < film_t.reynolds < 3000.0 and not film_t.in_range
    lam = inner_film(0.5, 0.2, 2e-3, 2000.0, 0.13, length_m=1000.0)               # Re 1592, Pr 30.8: L_e 490 m
    assert lam.in_range
    assert not inner_film(0.5, 0.2, 2e-3, 2000.0, 0.13, length_m=100.0).in_range  # developing over 100 m
    from angelica.core.case import FlowBoundary, NetworkCase, PressureBoundary, ThermalBoundary
    from angelica.properties.thermal_fluid import ThermalFluid
    from angelica.solvers import SteadyNonIsothermalIncompressibleSolver
    fluid = ThermalFluid.from_constants(density_kg_per_m3=880.0, viscosity_pa_s=2e-3,
                                        specific_heat_j_per_kg_k=2000.0, thermal_conductivity_w_per_m_k=0.13)
    case = NetworkCase(name="t", fluid_model=fluid, pressure_inlets=(), flow_inlets=(FlowBoundary(1, 0.8),),
                       pressure_outlets=(PressureBoundary(2, 1e5),), thermal_inlets=(ThermalBoundary(1, 60.0),),
                       components=(pipe(component_id="tr", ambient_temperature_c=5.0, compute_inner_film=True,
                                        wall_thickness_m=0.01, wall_thermal_conductivity_w_per_mk=45.0,
                                        outer_film_coefficient_w_per_m2k=20.0),))
    res = SteadyNonIsothermalIncompressibleSolver().solve(case)
    assert res.component_flows[0].inner_film_in_range is False
    assert any("inner film" in w and "tr" in w for w in res.warnings)


@pytest.mark.parametrize("m", [4.0, 6.0, 8.0, 10.0])
@pytest.mark.parametrize("outside", ["air", "buried"])
@pytest.mark.parametrize("n", [1, 20])
def test_viscous_oil_with_computed_film_converges(m, outside, n):
    """Heavy oil, μ = 0.05·e^(−0.05(T−20)) Pa·s: the film taken at the pipe's
    mean temperature flipped between laminar and turbulent and did not
    converge at 4–10 kg/s (v1.38.0 review).  At the inlet temperature (and
    under-relaxed) every case converges; the balance closes (< 0.01 %)."""
    from angelica.core.case import FlowBoundary, NetworkCase, PressureBoundary, ThermalBoundary
    from angelica.properties.thermal_fluid import ThermalFluid
    from angelica.solvers import SteadyNonIsothermalIncompressibleSolver
    fluid = ThermalFluid.from_functions(density_fn=lambda T: 900.0,
                                        viscosity_fn=lambda T: 0.05 * math.exp(-0.05 * (T - 20.0)),
                                        specific_heat_fn=lambda T: 2000.0, thermal_conductivity_fn=lambda T: 0.13)
    extra = ({"outer_film_coefficient_w_per_m2k": 25.0} if outside == "air"
             else {"burial_depth_m": 1.0, "soil_thermal_conductivity_w_per_mk": 1.5})
    case = NetworkCase(name="h", fluid_model=fluid, pressure_inlets=(), flow_inlets=(FlowBoundary(1, m),),
                       pressure_outlets=(PressureBoundary(2, 1e5),), thermal_inlets=(ThermalBoundary(1, 70.0),),
                       components=(pipe(length_m=5000.0, absolute_roughness_m=46e-6, ambient_temperature_c=5.0,
                                        n_thermal_segments=n, compute_inner_film=True, wall_thickness_m=0.01,
                                        wall_thermal_conductivity_w_per_mk=45.0, **extra),))
    res = SteadyNonIsothermalIncompressibleSolver().solve(case)
    assert res.converged
    assert abs(res.global_energy_balance.energy_error_pct) < 0.01
    f = res.component_flows[0]
    T_in = res.node_temperatures_c[1]
    mu = 0.05 * math.exp(-0.05 * (T_in - 20.0))
    film = inner_film(m, 0.2, mu, 2000.0, 0.13, 5000.0)
    assert f.inner_film_coefficient_w_per_m2k == pytest.approx(film.h_w_per_m2k, rel=1e-6)   # converged film


def test_gui_round_trip_keeps_the_wall_fields(tmp_path):
    from angelica.gui.io import build_network_case_from_scene, load_scene_from_file, save_scene_to_file
    scene = load_scene_from_file(ROOT / "src/angelica/tutorials/steady_black_oil/07_hot_oil_pipeline/"
                                 "07_hot_oil_pipeline.gui.json")
    comp = scene.links[0].components[0]
    scene.update_link_component_properties(scene.links[0].link_id, comp.component_id,
                                           {**comp.properties, "inner_film_coefficient_w_per_m2k": "Auto ",
                                            "wall_thickness_m": "0.01", "wall_thermal_conductivity_w_per_mk": "45",
                                            "burial_depth_m": "1.2", "soil_thermal_conductivity_w_per_mk": "1.4"})
    save_scene_to_file(scene, tmp_path / "s.gui.json")
    p = next(c for c in build_network_case_from_scene(load_scene_from_file(tmp_path / "s.gui.json")).components
             if isinstance(c, Pipe))
    assert p.compute_inner_film and p.burial_depth_m == 1.2 and p.wall_thickness_m == 0.01
