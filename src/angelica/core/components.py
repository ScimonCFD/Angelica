from __future__ import annotations

import math
from dataclasses import dataclass

FITTING_PRESET_LIBRARY = {
    "regular_90_flanged": {"name": "90 elbow, regular, flanged", "loss_coefficient": 0.3},
    "regular_90_threaded": {"name": "90 elbow, regular, threaded", "loss_coefficient": 1.5},
    "long_radius_90_flanged": {"name": "90 elbow, long radius, flanged", "loss_coefficient": 0.2},
    "long_radius_90_threaded": {"name": "90 elbow, long radius, threaded", "loss_coefficient": 0.7},
    "long_radius_45_flanged": {"name": "45 elbow, long radius, flanged", "loss_coefficient": 0.2},
    "regular_45_threaded": {"name": "45 elbow, regular, threaded", "loss_coefficient": 0.4},
    "return_bend_180_flanged": {"name": "180 return bend, flanged", "loss_coefficient": 0.2},
    "return_bend_180_threaded": {"name": "180 return bend, threaded", "loss_coefficient": 1.5},
    "tee_line_flow_flanged": {"name": "Tee line flow, flanged", "loss_coefficient": 0.2},
    "tee_line_flow_threaded": {"name": "Tee line flow, threaded", "loss_coefficient": 0.9},
    "tee_branch_flow_flanged": {"name": "Tee branch flow, flanged", "loss_coefficient": 1.0},
    "tee_branch_flow_threaded": {"name": "Tee branch flow, threaded", "loss_coefficient": 2.0},
    "union_threaded": {"name": "Union, threaded", "loss_coefficient": 0.08},
    "globe_valve_fully_open": {"name": "Globe valve, fully open", "loss_coefficient": 10.0},
    "angle_valve_fully_open": {"name": "Angle valve, fully open", "loss_coefficient": 2.0},
    "gate_valve_fully_open": {"name": "Gate valve, fully open", "loss_coefficient": 0.15},
    "gate_valve_quarter_closed": {"name": "Gate valve, 1/4 closed", "loss_coefficient": 0.26},
    "gate_valve_half_closed": {"name": "Gate valve, 1/2 closed", "loss_coefficient": 2.1},
    "gate_valve_three_quarters_closed": {
        "name": "Gate valve, 3/4 closed",
        "loss_coefficient": 17.0,
    },
    "swing_check_forward_flow": {"name": "Swing check, forward flow", "loss_coefficient": 2.0},
    "swing_check_backward_flow": {
        "name": "Swing check, backward flow",
        "loss_coefficient": float("inf"),
    },
    "ball_valve_fully_open": {"name": "Ball valve, fully open", "loss_coefficient": 0.05},
    "ball_valve_one_third_closed": {"name": "Ball valve, 1/3 closed", "loss_coefficient": 5.5},
    "ball_valve_two_thirds_closed": {"name": "Ball valve, 2/3 closed", "loss_coefficient": 210.0},
}


@dataclass(frozen=True)
class PressureChanger:
    start_node: int
    end_node: int
    diameter_m: float
    component_id: str = ""

    @property
    def area_m2(self) -> float:
        return math.pi * self.diameter_m**2 / 4.0


@dataclass(frozen=True)
class Pipe(PressureChanger):
    length_m: float = 0.0
    absolute_roughness_m: float = 0.0
    hazen_williams_c: float = 130.0
    height_change_m: float = 0.0
    heat_transfer_coefficient_w_per_m2k: float = 0.0
    ambient_temperature_c: float = 20.0
    n_thermal_segments: int = 10
    # Multi-layer composite U (resistances in series).
    # When any layer parameter is set, these override heat_transfer_coefficient_w_per_m2k.
    inner_film_coefficient_w_per_m2k: float = 0.0
    wall_thickness_m: float = 0.0
    wall_thermal_conductivity_w_per_mk: float = 0.0
    insulation_thickness_m: float = 0.0
    insulation_thermal_conductivity_w_per_mk: float = 0.0
    outer_film_coefficient_w_per_m2k: float = 0.0

    def __post_init__(self) -> None:
        if self.length_m <= 0.0:
            raise ValueError(f"Pipe length_m must be positive (got {self.length_m})")
        if abs(self.height_change_m) > self.length_m:
            raise ValueError(
                f"abs(height_change_m)={abs(self.height_change_m):.4f} exceeds "
                f"length_m={self.length_m:.4f} — pipe cannot rise or fall more than its own length"
            )

    @property
    def u_overall_w_per_m2k(self) -> float:
        """Effective overall heat-transfer coefficient (W/m²K).

        Computes the composite U from inner-film + wall + insulation + outer-film
        resistances in series (flat-wall approximation, valid for thin walls).
        Falls back to heat_transfer_coefficient_w_per_m2k when no layer parameters are set.
        """
        R = 0.0
        if self.inner_film_coefficient_w_per_m2k > 0.0:
            R += 1.0 / self.inner_film_coefficient_w_per_m2k
        if self.wall_thickness_m > 0.0 and self.wall_thermal_conductivity_w_per_mk > 0.0:
            R += self.wall_thickness_m / self.wall_thermal_conductivity_w_per_mk
        if self.insulation_thickness_m > 0.0 and self.insulation_thermal_conductivity_w_per_mk > 0.0:
            R += self.insulation_thickness_m / self.insulation_thermal_conductivity_w_per_mk
        if self.outer_film_coefficient_w_per_m2k > 0.0:
            R += 1.0 / self.outer_film_coefficient_w_per_m2k
        if R > 0.0:
            return 1.0 / R
        return self.heat_transfer_coefficient_w_per_m2k


@dataclass(frozen=True)
class Fitting(PressureChanger):
    loss_coefficient: float = float("nan")

    def __post_init__(self) -> None:
        if math.isnan(self.loss_coefficient):
            raise ValueError("Fitting requires loss_coefficient (e.g. loss_coefficient=2.0)")
        if self.loss_coefficient <= 0.0:
            raise ValueError(
                f"Fitting loss_coefficient must be positive (got {self.loss_coefficient})"
            )


@dataclass(frozen=True)
class Choke(Fitting):
    """Production choke (bean) with critical / subcritical two-phase flow.

    With ``flow_model="sachdeva"`` (default) and a black-oil fluid,
    ``SteadyBlackOilSolver`` computes the flow with the Sachdeva et al. (1986)
    model (see ``angelica.closures.choke``): when P2/P1 drops below the
    critical ratio the flow is choked and no longer depends on the downstream
    pressure.  With pure liquid the model reduces to the orifice equation below.

    With ``flow_model="orifice"``, or with fluids that have no black-oil PVT
    (water, dead oil, gas solvers), the choke behaves as a K-factor fitting
    derived from the incompressible orifice equation:

        ṁ = Cd · A_orifice · √(2·ρ·ΔP)

    which maps to a loss coefficient:

        K = (A_pipe / (Cd · A_orifice))²  =  (D_pipe / D_orifice)⁴ / Cd²

    so that the pressure drop is ΔP = K · ρ · v² / 2, exactly like a Fitting.

    Parameters
    ----------
    start_node, end_node, diameter_m : inherited (pipe cross-section)
    orifice_diameter_m : inner diameter of the choke bean/nozzle (m)
    discharge_coefficient : orifice discharge coefficient Cd (default 0.62,
        appropriate for a sharp-edged orifice plate per ISO 5167; Sachdeva
        et al. report Cd ≈ 0.75 for liquid-dominated and ≈ 0.85 for
        gas-dominated flow through production chokes)
    flow_model : "sachdeva" (default, critical-flow aware) or "orifice"

    Typical oilfield usage
    ----------------------
    In the US oilfield, choke size is quoted in 64ths of an inch ("beans").
    A 16/64" bean is orifice_diameter_m ≈ 0.00635 m.
    A 32/64" (½") bean is orifice_diameter_m ≈ 0.0127 m.
    """

    orifice_diameter_m: float = float("nan")
    discharge_coefficient: float = 0.62
    flow_model: str = "sachdeva"

    def __post_init__(self) -> None:
        if math.isnan(self.orifice_diameter_m):
            raise ValueError("Choke requires orifice_diameter_m (inner diameter of the nozzle/bean).")
        if self.orifice_diameter_m <= 0.0:
            raise ValueError(
                f"Choke orifice_diameter_m must be positive (got {self.orifice_diameter_m})"
            )
        if self.orifice_diameter_m >= self.diameter_m:
            raise ValueError(
                f"Choke orifice_diameter_m ({self.orifice_diameter_m:.4f} m) must be "
                f"less than pipe diameter_m ({self.diameter_m:.4f} m)"
            )
        if self.flow_model not in ("sachdeva", "orifice"):
            raise ValueError(
                f"Choke flow_model must be 'sachdeva' or 'orifice' (got {self.flow_model!r})"
            )
        if not (0.0 < self.discharge_coefficient <= 1.0):
            raise ValueError(
                f"discharge_coefficient must be in (0, 1] (got {self.discharge_coefficient})"
            )
        K = (self.diameter_m / self.orifice_diameter_m) ** 4 / self.discharge_coefficient ** 2
        object.__setattr__(self, "loss_coefficient", K)


@dataclass(frozen=True)
class CheckValve(Fitting):
    """Non-return (check) valve: a K-factor fitting that only lets flow pass
    from ``start_node`` to ``end_node``.

    Forward flow sees the loss coefficient ``loss_coefficient`` (default 2.0, a
    swing check valve fully open).  When the pressure at ``end_node`` exceeds
    the pressure at ``start_node`` the valve closes and carries no flow.  Works
    in every solver.

    Typical use: on each well of a gathering network, so that a low-pressure
    well is shut in instead of taking back-flow from the manifold.
    """

    loss_coefficient: float = 2.0
    blocks_reverse_flow = True   # class attribute read by the minor-loss model


@dataclass(frozen=True)
class WellInflow(Fitting):
    """Well inflow performance relationship (IPR) between a reservoir node and
    the bottom-hole node.

    ``start_node`` is the reservoir (a pressure boundary at the average
    reservoir pressure P_r, with the reservoir temperature), ``end_node`` the
    bottom hole (flowing pressure P_wf).  The stock-tank oil rate is

      linear:  q_o = J · (P_r − P_wf)
      vogel:   q_o = J · (P_r − P_wf)                        for P_wf ≥ P_b
               q_o = J · (P_r − P_b) + J·P_b/1.8 · (1 − 0.2·x − 0.8·x²),
                     x = P_wf / P_b                          for P_wf < P_b
               (composite Vogel; P_b is the bubble point of the well fluid at
               reservoir temperature, capped at P_r)

      fetkovich:     q_o = C · (P_r² − P_wf²)ⁿ        (pressures in bar)
      backpressure:  q_g = C · (P_r² − P_wf²)ⁿ        (gas rate, gas wells)

    with J = ``productivity_index_m3_per_day_per_bar`` (stock-tank m³ of oil
    per day per bar), C = ``ipr_coefficient`` (standard m³/d per bar^(2n), of
    oil or gas) and n = ``ipr_exponent`` (0.5–1).  The mass rate includes the
    produced gas and water, ṁ = q_o·(ρ_o,sc + GOR·ρ_g,sc + WOR·ρ_w,sc), with
    q_o = q_g/GOR for the back-pressure equation.  If P_wf > P_r the well
    takes injection by the same law (use a ``CheckValve`` to prevent it).

    Requires ``SteadyBlackOilSolver``.  ``diameter_m`` is the tubing bore used
    only to report a velocity.
    """

    productivity_index_m3_per_day_per_bar: float = float("nan")
    ipr_model: str = "vogel"
    ipr_coefficient: float = float("nan")
    ipr_exponent: float = 1.0
    loss_coefficient: float = 1.0   # unused; keeps the Fitting interface
    is_well_inflow = True           # class attribute read by the flow models

    def __post_init__(self) -> None:
        from angelica.closures.well_inflow import validate_ipr_fields
        validate_ipr_fields("WellInflow", self.ipr_model, self.productivity_index_m3_per_day_per_bar,
                            self.ipr_coefficient, self.ipr_exponent)


@dataclass(frozen=True)
class Well(Fitting):
    """Producing oil well: reservoir inflow (IPR) plus production tubing.

    ``start_node`` is the reservoir (a pressure boundary at the average
    reservoir pressure, with the reservoir temperature), ``end_node`` the
    wellhead.  For each oil rate the flowing bottom-hole pressure follows the
    IPR (see ``WellInflow``) and the tubing pressure drop is integrated upward
    with Beggs-Brill and the local black-oil PVT, giving the well
    deliverability curve P_wh(q).  The network sees the stable branch of that
    curve (dP_wh/dq < 0, the right-hand intersection of nodal analysis); if
    the wellhead pressure exceeds the curve's maximum the well is dead and
    produces nothing (it never takes injection).

    Parameters
    ----------
    productivity_index_m3_per_day_per_bar : J, stock-tank m³ oil / d / bar
        (linear and Vogel)
    ipr_model : "vogel" (default), "linear", "fetkovich" or "backpressure"
        (see ``WellInflow``)
    ipr_coefficient, ipr_exponent : C and n of the Fetkovich and
        back-pressure laws
    diameter_m : tubing inner diameter
    tubing_length_m : measured length of the tubing
    tubing_height_change_m : vertical rise from bottom hole to wellhead
        (defaults to ``tubing_length_m``, a vertical well)
    tubing_roughness_m : absolute roughness (default 4.6e-5 m)
    tubing_segments : integration steps along the tubing (default 25)
    tubing_correlation : "beggs_brill" (default), "hagedorn_brown" (oil wells)
        or "gray" (gas and gas-condensate wells)

    Tubing temperature (see ``angelica.closures.well_model``): the fluid
    enters the bottom hole at its reservoir enthalpy and loses heat to the
    formation on the way up (Ramey 1962); the well delivers its wellhead
    temperature to the network.

    geothermal_gradient_k_per_m : fall of the undisturbed formation
        temperature per metre of vertical depth above the reservoir
        (default 0.03 K/m; the formation is at the reservoir temperature at
        the bottom hole)
    formation_conductivity_w_per_mk : k_e (default 2.4; 0 = adiabatic tubing)
    formation_diffusivity_m2_per_s : α (default 1.0e-6)
    production_time_days : time on production, for the transient formation
        conduction f(t) (default 365)
    wellbore_diameter_m : diameter at the cement/formation interface
        (default: the larger of 0.2159 m — an 8½-in. hole — and 1.5 × tubing
        diameter)
    tubing_u_w_per_m2k : overall coefficient of the wellbore (tubing wall,
        annulus, casing, cement) referred to the tubing inner diameter;
        default inf = no wellbore resistance (formation conduction only)

    Requires ``SteadyBlackOilSolver``.
    """

    productivity_index_m3_per_day_per_bar: float = float("nan")
    ipr_model: str = "vogel"
    ipr_coefficient: float = float("nan")
    ipr_exponent: float = 1.0
    tubing_length_m: float = float("nan")
    tubing_height_change_m: float = float("nan")
    tubing_roughness_m: float = 4.6e-5
    tubing_segments: int = 25
    tubing_correlation: str = "beggs_brill"
    geothermal_gradient_k_per_m: float = 0.03
    formation_conductivity_w_per_mk: float = 2.4
    formation_diffusivity_m2_per_s: float = 1.0e-6
    production_time_days: float = 365.0
    wellbore_diameter_m: float = float("nan")
    tubing_u_w_per_m2k: float = float("inf")
    loss_coefficient: float = 1.0   # unused; keeps the Fitting interface
    is_well_model = True            # class attribute read by the flow models

    def __post_init__(self) -> None:
        from angelica.closures.well_inflow import validate_ipr_fields
        validate_ipr_fields("Well", self.ipr_model, self.productivity_index_m3_per_day_per_bar,
                            self.ipr_coefficient, self.ipr_exponent)
        if math.isnan(self.tubing_length_m) or self.tubing_length_m <= 0.0:
            raise ValueError("Well requires a positive tubing_length_m.")
        if math.isnan(self.tubing_height_change_m):
            object.__setattr__(self, "tubing_height_change_m", self.tubing_length_m)
        if abs(self.tubing_height_change_m) > self.tubing_length_m + 1e-9:
            raise ValueError("Well tubing_height_change_m cannot exceed tubing_length_m.")
        if self.tubing_segments < 2:
            raise ValueError("Well tubing_segments must be at least 2.")
        if self.tubing_correlation not in ("beggs_brill", "hagedorn_brown", "gray"):
            raise ValueError(
                "Well tubing_correlation must be 'beggs_brill', 'hagedorn_brown' or 'gray' "
                f"(got {self.tubing_correlation!r})")
        if math.isnan(self.wellbore_diameter_m):
            object.__setattr__(self, "wellbore_diameter_m", max(0.2159, 1.5 * self.diameter_m))
        if self.wellbore_diameter_m < self.diameter_m:
            raise ValueError("Well wellbore_diameter_m cannot be smaller than the tubing diameter.")
        if self.geothermal_gradient_k_per_m < 0.0:
            raise ValueError("Well geothermal_gradient_k_per_m must be non-negative.")
        if self.formation_conductivity_w_per_mk < 0.0:
            raise ValueError("Well formation_conductivity_w_per_mk must be non-negative.")
        if self.formation_diffusivity_m2_per_s <= 0.0 or self.production_time_days <= 0.0:
            raise ValueError("Well formation_diffusivity_m2_per_s and production_time_days must be positive.")
        if math.isnan(self.tubing_u_w_per_m2k) or self.tubing_u_w_per_m2k < 0.0:
            raise ValueError("Well tubing_u_w_per_m2k must be non-negative (inf = no wellbore resistance).")


@dataclass(frozen=True)
class Pump(PressureChanger):
    """Centrifugal pump (or compressor) defined by a head curve.

    ``curve_points_q_head`` holds ``(Q, H)`` pairs with Q the volumetric flow
    in **m³/h** (not m³/s) and H the head in metres of fluid.  A single point
    is treated as the design point of a standard three-term curve.
    """

    curve_points_q_head: tuple[tuple[float, float], ...] = ()

    def __post_init__(self) -> None:
        if not self.curve_points_q_head:
            raise ValueError("Pump curve_points_q_head must contain at least one (Q, head) point")
        points = tuple(sorted((float(q), float(head)) for q, head in self.curve_points_q_head))
        if any(q < 0.0 for q, _ in points):
            raise ValueError("Pump curve flow values must be non-negative.")
        if any(head < 0.0 for _, head in points):
            raise ValueError("Pump curve head values must be non-negative.")
        for i in range(1, len(points)):
            prev_q, prev_head = points[i - 1]
            curr_q, curr_head = points[i]
            if curr_q <= prev_q:
                raise ValueError("Pump curve flow values must be strictly increasing.")
            if curr_head >= prev_head:
                raise ValueError("Pump curve head values must decrease as flow increases.")
        object.__setattr__(self, "curve_points_q_head", points)


@dataclass(frozen=True)
class HeatSource(PressureChanger):
    """Inline heater (power_w > 0) or cooler (power_w < 0) with fixed thermal power.

    Hydraulic pressure drop is controlled by pressure_drop_mode:
      "rated"  — ΔP = ΔP_rated · (ṁ/ṁ_rated)²  (quadratic, derived from K factor)
      "fixed"  — ΔP = pressure_drop_pa (lagged coefficient, updated each iteration)

    Set pressure_drop_pa = 0 for a hydraulically transparent device.
    """

    power_w: float = 0.0
    pressure_drop_mode: str = "rated"
    pressure_drop_pa: float = 0.0
    rated_mass_flow_kg_per_s: float = 1.0
    n_thermal_segments: int = 10
    length_m: float = 1.0

    def __post_init__(self) -> None:
        if self.length_m <= 0.0:
            raise ValueError(f"HeatSource length_m must be positive (got {self.length_m})")
