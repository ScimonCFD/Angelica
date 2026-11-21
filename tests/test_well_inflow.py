"""
Tests for well inflow (IPR), check valves, composition propagation through
devices, and the hybrid convection scheme with reverse flow (v1.12.0).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures.well_inflow import (  # noqa: E402
    mass_per_oil_sc_m3,
    oil_rate_derivative,
    oil_rate_sc_m3_per_day,
)
from angelica.core.case import InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary  # noqa: E402
from angelica.core.components import CheckValve, Pipe, WellInflow  # noqa: E402
from angelica.properties.black_oil import BlackOilFluid  # noqa: E402
from angelica.properties.single_component import SingleComponentFluid  # noqa: E402
from angelica.solvers import (  # noqa: E402
    SteadyBlackOilSolver,
    SteadyIsothermalIncompressibleSolver,
)

BAR = 1e5


# ── IPR formulas ─────────────────────────────────────────────────────────────

class TestIprFormulas:
    def test_linear(self):
        assert oil_rate_sc_m3_per_day(80 * BAR, 50 * BAR, 10.0, "linear", 200 * BAR) == pytest.approx(300.0)

    def test_vogel_saturated_open_flow(self):
        """P_b ≥ P_r: q_max = J·P_r/1.8 (Vogel's absolute open flow)."""
        q = oil_rate_sc_m3_per_day(80 * BAR, 0.0, 10.0, "vogel", 200 * BAR)
        assert q == pytest.approx(10.0 * 80 / 1.8)

    def test_vogel_equals_linear_above_bubble_point(self):
        q_v = oil_rate_sc_m3_per_day(200 * BAR, 150 * BAR, 5.0, "vogel", 100 * BAR)
        assert q_v == pytest.approx(5.0 * 50)

    def test_composite_vogel_continuous_at_bubble_point(self):
        args = (200 * BAR, 5.0, "vogel", 100 * BAR)
        above = oil_rate_sc_m3_per_day(args[0], 100 * BAR + 1.0, *args[1:])
        below = oil_rate_sc_m3_per_day(args[0], 100 * BAR - 1.0, *args[1:])
        assert above == pytest.approx(below, rel=1e-6)

    def test_vogel_below_linear(self):
        for p_wf in (10, 40, 70):
            assert (oil_rate_sc_m3_per_day(80 * BAR, p_wf * BAR, 10.0, "vogel", 100 * BAR)
                    < oil_rate_sc_m3_per_day(80 * BAR, p_wf * BAR, 10.0, "linear", 100 * BAR))

    def test_injection_when_flowing_pressure_exceeds_reservoir(self):
        assert oil_rate_sc_m3_per_day(80 * BAR, 90 * BAR, 10.0, "vogel", 200 * BAR) == pytest.approx(-100.0)

    def test_derivative_matches_finite_difference(self):
        p_r, J, pb = 80 * BAR, 10.0, 200 * BAR
        for p_wf in (10 * BAR, 40 * BAR, 70 * BAR):
            h = 10.0
            fd = (oil_rate_sc_m3_per_day(p_r, p_wf - h, J, "vogel", pb)
                  - oil_rate_sc_m3_per_day(p_r, p_wf + h, J, "vogel", pb)) / (2 * h)
            assert oil_rate_derivative(p_r, p_wf, J, "vogel", pb) == pytest.approx(fd, rel=1e-6)

    def test_mass_per_stock_tank_m3(self):
        fluid = BlackOilFluid(api_gravity=35.0, gas_gravity=0.65, gor_sc_m3_per_m3=100.0,
                              wor_sc_m3_per_m3=0.5)
        rho_o = 141.5 / (131.5 + 35.0) * 999.0
        expected = rho_o + 100.0 * 0.65 * 1.225 + 0.5 * 1025.0
        assert mass_per_oil_sc_m3(fluid) == pytest.approx(expected, rel=2e-3)

    def test_component_validation(self):
        with pytest.raises(ValueError):
            WellInflow(start_node=1, end_node=2, diameter_m=0.1)
        with pytest.raises(ValueError):
            WellInflow(start_node=1, end_node=2, diameter_m=0.1,
                       productivity_index_m3_per_day_per_bar=5.0, ipr_model="fetkovich")


# ── IPR inside the black-oil solver ──────────────────────────────────────────

def _single_well_case(model: str, p_res=12e6, p_wf=5e6) -> tuple[NetworkCase, BlackOilFluid]:
    fluid = BlackOilFluid(api_gravity=32.0, gas_gravity=0.65, gor_sc_m3_per_m3=60.0,
                          wor_sc_m3_per_m3=0.2, reference_pressure_pa=8e6, reference_temperature_c=70.0)
    case = NetworkCase(
        name="single-well",
        fluid_model=fluid,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=p_res),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=p_wf),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=70.0, bc_type="fixed_temperature"),),
        components=(WellInflow(start_node=1, end_node=2, diameter_m=0.1016,
                               productivity_index_m3_per_day_per_bar=8.0, ipr_model=model),),
    )
    return case, fluid


@pytest.mark.parametrize("model", ["linear", "vogel"])
def test_solver_reproduces_ipr(model):
    case, fluid = _single_well_case(model)
    res = SteadyBlackOilSolver().solve(case)
    assert res.converged
    p_b = fluid.pvt(12e6, 70.0).bubble_point_pa
    q = oil_rate_sc_m3_per_day(12e6, 5e6, 8.0, model, p_b)
    m_expected = q * mass_per_oil_sc_m3(fluid) / 86_400.0
    assert res.component_flows[0].mass_flow_kg_per_s == pytest.approx(m_expected, rel=1e-6)


def test_well_inflow_rejected_outside_black_oil():
    case = NetworkCase(
        name="ipr-water",
        fluid_model=SingleComponentFluid(density_kg_per_m3=1000.0, viscosity_pa_s=1e-3),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=3e5),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=1e5),),
        components=(WellInflow(start_node=1, end_node=2, diameter_m=0.1,
                               productivity_index_m3_per_day_per_bar=8.0),),
    )
    with pytest.raises(ValueError, match="black-oil"):
        SteadyIsothermalIncompressibleSolver().solve(case)


# ── Check valve ──────────────────────────────────────────────────────────────

def _check_valve_case(p_branch: float) -> NetworkCase:
    """Source 1 (3 bar) → junction 2 → sink 3 (1 bar), plus a branch
    2 → check valve → 4 → pipe → sink 5 at ``p_branch``."""
    water = SingleComponentFluid(density_kg_per_m3=1000.0, viscosity_pa_s=1e-3)
    return NetworkCase(
        name="check-valve",
        fluid_model=water,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=3e5),),
        pressure_outlets=(PressureBoundary(node_id=3, pressure_pa=1e5),
                          PressureBoundary(node_id=5, pressure_pa=p_branch)),
        components=(
            Pipe(start_node=1, end_node=2, diameter_m=0.1, length_m=100.0, absolute_roughness_m=4.5e-5),
            Pipe(start_node=2, end_node=3, diameter_m=0.1, length_m=100.0, absolute_roughness_m=4.5e-5),
            CheckValve(component_id="cv", start_node=2, end_node=4, diameter_m=0.08),
            Pipe(start_node=4, end_node=5, diameter_m=0.08, length_m=50.0, absolute_roughness_m=4.5e-5),
        ),
    )


def _cv_flow(res) -> float:
    return next(c.mass_flow_kg_per_s for c in res.component_flows if c.label.endswith("cv"))


def test_check_valve_passes_forward_flow():
    res = SteadyIsothermalIncompressibleSolver().solve(_check_valve_case(1e5))
    assert res.converged
    assert _cv_flow(res) > 1.0


def test_check_valve_blocks_reverse_flow():
    """The branch sink at 4 bar would push water back into the junction; the
    check valve stops it while the main path keeps flowing."""
    res = SteadyIsothermalIncompressibleSolver().solve(_check_valve_case(4e5))
    assert res.converged
    assert abs(_cv_flow(res)) < 1e-6
    assert res.component_flows[1].mass_flow_kg_per_s > 1.0


# ── Wells with check valves in the black-oil solver (tutorial 11 layout) ─────

def _two_well_network(check: bool):
    wells = ((10, 8e6, 40.0, 0.70, 150.0, 0.1), (20, 2.4e6, 25.0, 0.60, 15.0, 0.8))
    comps, inlets, thermal, fluids = [], [], [], []
    for res, p, api, gg, gor, wor in wells:
        wh, down = res + 1, res + 2
        inlets.append(PressureBoundary(node_id=res, pressure_pa=p))
        thermal.append(ThermalBoundary(node_id=res, temperature_c=60.0, bc_type="fixed_temperature"))
        fluids.append(InletFluidBC(node_id=res, api_gravity=api, gas_gravity=gg,
                                   gor_sc_m3_per_m3=gor, wor_sc_m3_per_m3=wor))
        comps.append(WellInflow(component_id=f"ipr{res}", start_node=res, end_node=wh, diameter_m=0.1016,
                                productivity_index_m3_per_day_per_bar=20.0))
        start = wh
        if check:
            comps.append(CheckValve(component_id=f"cv{res}", start_node=wh, end_node=down, diameter_m=0.1016))
            start = down
        comps.append(Pipe(component_id=f"fl{res}", start_node=start, end_node=100, diameter_m=0.1016,
                          length_m=2000.0, absolute_roughness_m=46e-6))
    comps.append(Pipe(component_id="trunk", start_node=100, end_node=101, diameter_m=0.1016,
                      length_m=6000.0, absolute_roughness_m=46e-6))
    return NetworkCase(
        name="two-wells", fluid_model=BlackOilFluid(api_gravity=32.0, gas_gravity=0.65,
                                                    gor_sc_m3_per_m3=50.0, wor_sc_m3_per_m3=0.3),
        pressure_inlets=tuple(inlets), pressure_outlets=(PressureBoundary(node_id=101, pressure_pa=1.5e6),),
        thermal_inlets=tuple(thermal), inlet_fluid_bcs=tuple(fluids), components=tuple(comps),
    )


def test_weak_well_takes_backflow_without_check_valve():
    res = SteadyBlackOilSolver().solve(_two_well_network(check=False))
    assert res.converged
    flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in res.component_flows}
    assert flows["ipr20"] < -1e-3


def test_check_valve_shuts_in_weak_well_and_balances_close():
    res = SteadyBlackOilSolver().solve(_two_well_network(check=True))
    assert res.converged
    flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in res.component_flows}
    assert abs(flows["ipr20"]) < 1e-6
    assert res.node_pressures_pa[21] == pytest.approx(2.4e6, rel=1e-6)   # wellhead at reservoir pressure
    assert res.global_energy_balance.energy_error_pct < 1e-3


def test_composition_propagates_through_devices():
    """The flowline downstream of well 10's IPR and check valve must carry
    well 10's crude (40 °API, GOR 150), not the case default."""
    res = SteadyBlackOilSolver().solve(_two_well_network(check=True))
    cf = next(c for c in res.component_flows if c.label.endswith("fl10"))
    rho_link = cf.mass_flow_kg_per_s * 3600.0 / cf.volumetric_flow_m3_per_h
    P = 0.5 * (res.node_pressures_pa[12] + res.node_pressures_pa[100])
    T = res.node_temperatures_c[12]
    well = BlackOilFluid(api_gravity=40.0, gas_gravity=0.70, gor_sc_m3_per_m3=150.0, wor_sc_m3_per_m3=0.1)
    default = BlackOilFluid(api_gravity=32.0, gas_gravity=0.65, gor_sc_m3_per_m3=50.0, wor_sc_m3_per_m3=0.3)
    rho_well = well.pvt(P, T).mixture_density_kg_per_m3
    rho_default = default.pvt(P, T).mixture_density_kg_per_m3
    assert abs(rho_link - rho_well) < 0.25 * abs(rho_link - rho_default)


# ── Hybrid scheme with reverse flow ──────────────────────────────────────────

def test_hybrid_scheme_reverse_flow_matches_forward():
    """A pipe defined against the flow must give the same outlet temperature
    and a closed energy balance (the hybrid coefficients were wrong for F < 0)."""
    from angelica import (
        NonIsothermalSolverSettings,
        SteadyNonIsothermalIncompressibleSolver,
        ThermalFluid,
    )
    from angelica.closures import HybridScheme

    fluid = ThermalFluid.from_constants(density_kg_per_m3=900.0, viscosity_pa_s=0.01,
                                        specific_heat_j_per_kg_k=2000.0,
                                        thermal_conductivity_w_per_m_k=0.13)

    def solve(reverse: bool):
        a, b = (2, 1) if reverse else (1, 2)
        case = NetworkCase(
            name="hybrid", fluid_model=fluid,
            pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=3e5),),
            pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=1e5),),
            thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=80.0),),
            components=(Pipe(start_node=a, end_node=b, diameter_m=0.05, length_m=500.0,
                             absolute_roughness_m=4.5e-5, heat_transfer_coefficient_w_per_m2k=20.0,
                             ambient_temperature_c=10.0, n_thermal_segments=10),),
        )
        solver = SteadyNonIsothermalIncompressibleSolver(
            non_isothermal_settings=NonIsothermalSolverSettings(temperature_tolerance_k=1e-8,
                                                                max_temperature_iterations=50),
            convection_scheme=HybridScheme())
        return solver.solve(case)

    fwd, rev = solve(False), solve(True)
    assert rev.node_temperatures_c[2] == pytest.approx(fwd.node_temperatures_c[2], abs=1e-6)
    assert rev.global_energy_balance.energy_error_pct < 1e-6
    assert math.isfinite(rev.node_temperatures_c[2])
