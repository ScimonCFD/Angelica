# Validation manual

This manual lists what each part of Angelica has been checked against, with
the reference, the kind of evidence, the error obtained and the test that
reproduces it.  It is generated from the repository itself — test docstrings
and assertions, module docstrings, tutorial READMEs, `ROADMAP.md`,
`STATUS.md` and `CHANGELOG.md` — and quotes them; nothing in it comes from
outside the repository.  Every row names the test to rerun:

```bash
pip install -e .[dev]
pytest tests/test_hydrate.py::test_methane_against_semenov_2024 -q    # one row
pytest -q                                                              # everything (~17 min)
```

**Read this first.**  Most of Angelica is verified by *internal
consistency* (an independent code path, an exact limit, a conservation law)
or *analytically*, not against measured data.  Comparisons with data read in
a published source are few and listed as such; where the repository only
cites a value, or says it was recalled from memory, the manual says so.
Section [Not validated and known deviations](#not-validated-and-known-deviations)
lists what has no external check and every recorded discrepancy.

## Kinds of evidence

`ROADMAP.md` §4 "Standing rules":

- "No parameter is tuned to make a validation pass. A mismatch is explained or recorded as a limitation."
- "Validation data are labelled in the tests: read in the source, read in a secondary source, or from memory. Memory does not count as validation."
- "Each validation test states the error obtained and where its tolerance comes from (data uncertainty, published accuracy of the method); a tolerance is never relaxed to make a test pass."
- "At the end of each stage the log summarises: validated against data read in a source / validated only by internal consistency / not validated."
- "(Review 2026-11-25) Every new component is also run outside its tutorial's range before tagging: zero flow, reversed flow, dead wells, set-points that cannot be met. Each case must end converged or with converged=False and a stated reason, never with an exception."

Further wording used in the repo and kept below: "ANALYTICAL", "CONSISTENCY",
"regression guard" (tolerance chosen after the results; not an accuracy
claim), "COMPARED, NOT VALIDATED", "CROSS-CHECK, not validation".

Labels used in the tables (column "Label"):

| Label | Meaning (repo wording) |
|---|---|
| read in source | data "READ IN SOURCE" / "read in the source" |
| secondary source | "READ IN A SECONDARY SOURCE" / "secondary source" |
| analytical | exact closed-form solution ("ANALYTICAL") |
| consistency | internal consistency / independent code path ("CONSISTENCY") |
| from memory | repo says the reference value was recalled or "WRITTEN FROM MEMORY" (does not count) |
| not labelled | the test cites a value but does not say whether it was read in a source |
| regression | stored reference produced by Angelica itself |

Overall standing stated by the repo (`ROADMAP.md` §1, "Validation debt (2026-11-25)"):
"Most models since v1.19 are validated by internal consistency only: of the last
ten log summaries, eight say 'validated against data read in a source: none'.
Data read in a source cover only the black-oil enthalpy (∂h/∂P)_T of gas and
water and the methane solubility in water. Without a source: Hagedorn-Brown,
Gray, Ramey, Fetkovich and gas back-pressure; gas lift and compositional wells;
the separator." Later releases added data read in source for: water heater duty
(NIST), unit factors (NIST SP 811), water saturation pressure (NIST), methane
hydrate and methanol-inhibitor data (Semenov et al. 2023/2024). "The worked
examples that would close this (Brill & Mukherjee, Economides et al., GPSA)
are not in the repository's reference folder (`refs/` does not exist).
Recalling them from memory does not count."

Row count by primary label in the tables below (193 rows): consistency 100,
analytical 33, not labelled 19, read in source 15 (one of them published
output of another solver, EPANET; several are comparisons with
same-laboratory data or parameters read in source), secondary source 12,
regression 6 (incl. Hanoi), not validated / not done 6, from memory 2.

---

## 1. Single-phase friction and pipe velocity loop

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Laminar flow in three parallel smooth pipes (oil 880 kg/m³, 0.1 Pa·s, ΔP 50 kPa, Re 58–179) | Hagen-Poiseuille Q = πD⁴ΔP/(128μL): 0.7069, 0.2796, 1.4381 m³/h | analytical | Test tolerance 1e-6 m³/h on 0.706858347, 0.279567712, 1.438106675 m³/h; README: "Max relative error vs analytical < 0.001 %". Tolerance source not stated. | `tests/test_tutorial_suite.py::TutorialSuiteSmokeTests::test_laminar_parallel_pipes_match_poiseuille`; `src/angelica/tutorials/steady_isothermal_incompressible/10_laminar_poiseuille_benchmark/benchmark_tables.py` |
| Stagnant pipe velocity (ΔP → 0) | Hagen-Poiseuille v = D²·ΔP/(32·μ·L) | analytical | "agreement is to round-off (checked to 6 decimals)" | `tests/test_friction.py::ColebrookFrictionTests::test_stagnant_pipe_carries_no_phantom_flow` |
| Velocity loop solution satisfies Darcy-Weisbach (water, 200 mm, 50 kPa/km) | ΔP ≈ f·(L/D)·(ρ/2)·v² + elevation | consistency | Tolerance 1 % of ΔP; fixed-point and secant agree within 1 %; v within 0.5–10 m/s (docstring expects 1–5 m/s). Obtained not stated. | `tests/test_velocity_loops.py::TestVelocityLoopPhysics::test_darcy_weisbach_satisfied_at_solution`, `TestVelocityLoopMethods::test_fixed_point_and_secant_agree`, `TestVelocityLoopMethods::test_velocity_is_physically_reasonable` |
| Hazen-Williams single pipe (D 0.05 m, L 100 m, C 130, 200→100 kPa) | Head-loss form h_f = 10.67·L·Q^1.852/(C^1.852·D^4.871) | analytical | Relative flow error < 1e-4 ("exceeds 0.01%" message). Obtained not stated. | `tests/test_gui_io.py::GuiIoTests::test_hazen_williams_single_pipe_matches_analytical`; `tests/test_friction.py::HazenWilliamsPipeCorrelationTests` (sign/symmetry/coupling only) |
| Calibration multiplier = friction length scaling | multiplier c ≡ length L·c | consistency | Tolerance 1e-9 ("review measured ≤ 3e-16") | `tests/test_calibration.py::test_multiplier_equals_scaling_the_friction_length` |
| Laminar–turbulent jump at Re 2300 | 64/Re = 0.0278 vs Colebrook 0.0475 (ε/D = 3e-4) | consistency (diagnosis only) | "Not fixed (a modelling choice)"; a viscous line near Re 2300 can have no solution; Re cycles 2184 → 2306 → 3219 → 2334 → 2566 (f 0.029 ↔ 0.047) | `tests/test_black_oil_range.py::test_laminar_turbulent_cycle_is_diagnosed`, `test_ambient_sweep_of_the_regulated_line` |

## 2. Network hydraulics (isothermal incompressible)

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Hanoi network: 32 trunk-node heads and 34 trunk-link flows (Hazen-Williams, C = 130) | Values generated by Angelica itself (tutorial README); topology and demands from Fujiwara & Khang (1990), which publish no heads or flows | regression | ±0.05 m head (the baseline's 2 printed decimals), ±1 m³/h flow | `tests/test_tutorial_suite.py::HanoiBenchmarkTests::test_hanoi_benchmark_matches_regression_baseline` |
| EPANET 2.2 tutorial network with pump (one-point curve 1000 gpm / 200 ft, Hazen-Williams): node pressures, heads, 7 link flows | EPANET 2.2 user manual, Appendix C: Listing 1 (input, pump curve "1 1000 200") and Listing 2 (report at 0:00) — read 2026-11-28; another solver's published output, not measurements | read in source | Pressures within 0.212 psi (node 2; others ≤ 0.047), heads ≤ 0.022 ft, flows ≤ 0.075 gpm; tolerances 0.25 psi, 0.2 ft, 1.5 gpm (no stated source) | `tests/test_pump_benchmark.py::PumpBenchmarkTests::test_epanet_pump_benchmark_matches_reference` |
| Three-reservoir junction (heads 30/20/5 m, Colebrook) | Independent analytical Colebrook solve in `benchmark_tables.py` ("Analytical helpers"); README values H_J 25.26 m, 112.35 / 62.36 / 49.99 m³/h | analytical (independent solve) | Obtained ≤ 0.00074 m³/h (benchmark_tables.py); test tolerances 0.01 m and 0.1 m³/h; mass balance error at J < 1e-3 | `tests/test_tutorial_suite.py::TutorialSuiteSmokeTests::test_three_reservoir_junction_converges`; `src/angelica/tutorials/steady_isothermal_incompressible/09_three_reservoir_junction/benchmark_tables.py` |
| Flow-boundary vs pressure-boundary cases give the same solution (tutorials 04, 05, 06 vs 02/03) | Converged flows of tutorials 02/03 | consistency | Not quantified in README | `src/angelica/tutorials/steady_isothermal_incompressible/04_inlet_flow/README.md` (and 05, 06) |
| Two-loop network, every fluid mode: junction mass balance | Solver criterion `nodal_mass_imbalance_rel_tolerance` = 1e-3 | consistency | Obtained ≤ 1.07e-5 | `tests/test_two_loop_networks.py::test_converged_with_mass_balance_at_every_junction` |
| Two-loop network: each pipe solved alone reproduces the network downstream pressure (Kirchhoff loop check) | Same correlation, independent single-pipe solve | consistency | Tolerance 1e-4 of the pipe ΔP (velocity loop tolerance, outer density criterion 1e-4). Obtained ≤ 3.02e-5 (black-oil thermal; compositional isothermal 3.00e-5), ≤ 6e-8 for liquids | `tests/test_two_loop_networks.py::test_each_pipe_alone_reproduces_the_network_pressure` |
| Crude-oil pipeline (32 °API, 65 °C) | README "Expected Results" 43.00 / 26.69 / 16.31 m³/h, 2.547 bar | regression | Test asserts 43.0 / 26.7 / 16.3 ± 0.1 m³/h | `tests/test_tutorial_suite.py` (crude-oil tutorial checks); `src/angelica/tutorials/steady_isothermal_incompressible/11_crude_oil_pipeline/README.md` |

## 3. Pumps

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| EPANET one-point curve completion (shutoff 133 % of design head, zero head at 2× design flow) | EPANET rule as stated in README | consistency | Design point 227.124707 m³/h (=1000 gpm) to 3 decimals; shut-off above shutoff head | `tests/test_pump.py::PumpCurveModelTests::test_single_point_epanet_curve_hits_design_point`, `test_single_point_epanet_curve_shuts_off_above_shutoff_head` |
| Pump heating, incompressible: ΔT = ΔP/(ρ·cp) | Shaft-work formula | analytical | Tolerance 1e-4 K; tutorial 08 "matches the solver result exactly" | `tests/test_non_isothermal.py::test_pump_heating_matches_shaft_work_incompressible`; `src/angelica/tutorials/steady_non_isothermal_incompressible/08_pump_heating/run.py` |
| Pump in compressible PR solver: enthalpy rise bounded by ΔP/ρ_out ≤ Δh ≤ ΔP/ρ_in | Bounds | consistency | Δh within [0.90·ΔP/ρ_out, 1.10·ΔP/ρ_in]; heating > 30 K | `tests/test_compressible.py::QuantitativeBenchmarks::test_pump_heating_compressible_pr` |
| Head-curve "pump as compressor" vs isentropic `Compressor` (tutorial 06) | Gas isentrope | consistency | "92.4 kJ/kg where the gas isentrope gives 108.5 (+17.5 %) and the gas leaves 6.5 K hotter" | `STATUS.md` v1.30.0 |

## 4. Compressible gas pipelines and gas correlations

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Single pipe mass flow | P₁²−P₂² momentum equation + Colebrook-White | analytical | Within 1 % | `tests/test_compressible.py::QuantitativeBenchmarks::test_single_pipe_mass_flow_vs_p_squared` |
| Two identical isothermal pipes in series, midpoint pressure | √((P₁²+P₂²)/2) | analytical | Tolerance 500 Pa | `tests/test_compressible.py::QuantitativeBenchmarks::test_two_pipes_series_intermediate_pressure` |
| Flow BC vs pressure BC (tutorial 02) | Same network with pressure BCs | consistency | Outlets 500.04 kPa vs 500 kPa; test tolerance 1000 Pa | `tests/test_tutorial_suite.py::CompressibleTutorialTests::test_t02_flow_bc_cross_validation`; tutorial README 02 |
| Ideal gas density | ρ = PM/RT | analytical | 6 decimal places | `tests/test_compressible.py::IdealGasEOSTests::test_density_at_standard_conditions` |
| PR EOS at low pressure vs ideal gas (methane, 1 atm, 20 °C) | Ideal gas | consistency | Within 1 % | `tests/test_compressible.py::PengRobinsonEOSTests::test_low_pressure_converges_to_ideal_gas` |
| Lee-Gonzalez-Eakin methane viscosity, 1 atm 20 °C | "accepted value ≈ 1.1e-5 Pa·s" | not labelled | Tolerance ±0.1e-5 Pa·s | `tests/test_compressible.py::LeeGonzalezEakinTests::test_atmospheric_value_close_to_known` |
| Lee-Gonzalez-Eakin accuracy | Docstring: "Typical error < 2 % for natural gas in the range 38–171 °C, SG 0.52–0.80" | not labelled (no test) | — | `src/angelica/properties/gas_correlations.py` |
| Weymouth, Panhandle A/B friction factors | Docstrings: Panhandle A "Empirical for Re ≈ 5 × 10⁶ – 11 × 10⁶; efficiency factor E = 1 assumed" | not labelled (no test found) | — | `src/angelica/closures/friction.py` |
| Methane critical constants used in tests | "Methane critical properties (NIST)" Tc 190.564 K, Pc 4 599 200 Pa, ω 0.01141 | not labelled | — | `tests/test_compressible.py` (module constants) |

## 5. Liquid property correlations (water, dead oil)

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Liquid water density, viscosity, cp, k (0–100 °C) | Docstring: Kell (1975) "error < 0.1 kg/m³"; Vogel "matches IAPWS 2008 within 1 %"; cp constant 4182 "(varies < 1 % between 15 and 80 °C)"; k "error < 0.5 %" | not labelled (no test cited) | — | `src/angelica/properties/water_liquid.py` |
| Water cp vs NIST | NIST cp: largest deviation from 4182 between 5 and 80 °C is 0.51 % (4.20333 J/(g·K) at 5 °C) | read in source | Used as the heater tolerance (see §20) | `tests/test_heater_outlet.py` (module docstring) |
| Black-oil water viscosity at 20 °C | "Known value: ≈ 1 cP" | not labelled | Tolerance ±0.15 cP | `tests/test_black_oil_pvt.py::TestViscosities::test_water_viscosity_at_20c` |
| Dead-oil density from API | SG = 141.5/(API+131.5), ρ = SG·999.064 | analytical | 875.5 and 824.5 kg/m³, ±0.5 | `tests/test_dead_oil.py::DeadOilDensityTests::test_30api_density`, `test_40api_density` |
| Beggs-Robinson dead-oil viscosity, 30 °API at 60 °F | Comment "Reference: ~180 cP" | not labelled | Accepted range 50–500 cP | `tests/test_dead_oil.py::DeadOilViscosityTests::test_beggs_robinson_30api_60f` |
| Beggs-Robinson data range: 16–58 °API, 70–295 °F (21.1–146.1 °C), 460 observations | "pengtools wiki; IHS Harmony reference"; original SPE-5434-PA "not accessible"/"not read" | secondary source | "Not a validation of the correlation: the tests check that the solver says when it is outside that range" | `tests/test_black_oil_range.py::test_warnings_outside_the_data_range_and_none_inside`; `src/angelica/properties/dead_oil.py` |
| Watson-Nelson cp (30 °API, 60 °C) and Cragoe k (30 °API, 20 °C) | "typical crude oil cp is 1600–2500 J/(kg·K)"; "typical crude oil k is 0.10–0.16 W/(m·K)" | not labelled | Range checks (k tested 0.08–0.20) | `tests/test_dead_oil.py::DeadOilSpecificHeatTests::test_30api_at_60c_in_expected_range`, `DeadOilThermalConductivityTests::test_30api_at_20c_in_expected_range` |

## 6. Black-oil PVT correlations

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Hall-Yarborough z at T_pr = 1.807, P_pr = 2.984 (γ_g 0.65) | Papay (1968) formula, "≈ 0.905"; "agreement within 1–2 % is expected" | not labelled (formula cross-check) | Tolerance z = 0.90 ± 0.05 | `tests/test_black_oil_pvt.py::TestZFactorHallYarborough::test_z_against_papay_reference` (a formula cross-check, not chart readings — class docstring corrected in v1.39.0) |
| Standing bubble point: API 35, γ_g 0.65, 150 °F, GOR 200 scf/STB | "Standing (1947) ... Expected ≈ 1030 psia" | not labelled | Tolerance ±60 psia ("±5 %") | `tests/test_black_oil_pvt.py::TestBubblePoint::test_pb_standing_example` |
| Gas FVF | Bg = P_sc·T/(T_sc·P) for z = 1 | analytical | 1e-10 | `tests/test_black_oil_pvt.py::TestGasFVF::test_bg_ideal_gas_relation` |
| Rs, Bo, Bw, viscosity trends | Physical ordering | consistency | — | `tests/test_black_oil_pvt.py` |
| Black-oil vs compositional, same volatile oil, linear IPR, 2 500 m tubing | Black-oil GOR/API/γ_g from the EOS stream's standard flash | consistency (comparison) | Mass per stock-tank m³ 0.1 %; wellhead pressures: oil 1 within 2–5 %, oil 2 EOS 13–22 % higher ("Nothing adjusted") | `ROADMAP.md` log v1.24.0 |
| Nodal analysis with Vasquez-Beggs | Case correlation inherited | consistency | Review: "+3.9 % on rate with Standing in its place" (fixed) | `tests/test_studies.py::test_nodal_keeps_the_case_pvt_correlation` |

## 7. Black-oil energy balance (enthalpy model, Joule-Thomson, Katz hydrate chart)

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Gas (∂h/∂P)_T = −μ_JT·cp, methane-gravity gas, 5 and 10 MPa, 280–320 K | NIST Chemistry WebBook isobaric tables (methane: Setzmann & Wagner 1991, Cp uncertainty ~1 %) | read in source | Errors (v1.20.0): −2.0 to −2.5 % at 5 MPa, −3.9 to −6.1 % at 10 MPa. Tolerance 10 %: "set with the results known and is documented as such" (gravity-based Z correlation; "a few-percent error is expected") | `tests/test_black_oil_energy.py::test_gas_isothermal_throttling_coefficient_vs_nist` |
| Water (∂h/∂P)_T, 10 MPa, 300–360 K | NIST WebBook (IAPWS-95, 0.1 %) | read in source | Errors −1.9 to −2.8 %; tolerance 10 % (as above) | `tests/test_black_oil_energy.py::test_water_isothermal_throttling_coefficient_vs_nist` |
| Adiabatic black-oil pipe h_in = h_out + g·Δz | Model's own enthalpy (exact first law) | consistency | Errors ≤ 0.002 K; tolerance 0.01 K = outer-loop tolerance; balance residual < ṁ·cp·0.01 K | `tests/test_black_oil_energy.py::test_adiabatic_pipe_conserves_h_plus_gz` |
| Adiabatic liquid friction heating | ΔT = ΔP·(1 − T·β)/(ρ·cp), bound ΔP/(ρ·cp) ("physics bound, not a fit") | analytical bound | Obtained 2.03 K against a 2.93 K bound (T·β ≈ 0.31) | `tests/test_black_oil_energy.py::test_dead_oil_warms_gassy_oil_cools` |
| Black-oil choke is isenthalpic | Fluid's own h(P, T) | consistency | Tolerance 50 J/kg | `tests/test_black_oil_energy.py::test_choke_is_isenthalpic_in_black_oil` |
| Katz gas-gravity hydrate chart reproduces digitised nodes | "Katz 1945, GPSA", READ IN A SECONDARY SOURCE (Janna 2015, UTP dissertation, Tables 4.1a–4.5a) | secondary source | Exact at nodes (rel 1e-9; 1e-6 °F) | `tests/test_black_oil_energy.py::test_katz_chart_reproduces_its_digitised_values` |
| Black-oil enthalpy continuous through the −6.7 °C numerical floor | Continuity | consistency | Let-down 7.97 → 2 MPa below the floor: +0.12 K at −6.2, −8, −10 °C (was +16.6 K, −11 K); tolerance 0.02 K ("ten times the inversion's resolution") | `tests/test_black_oil_range.py::test_enthalpy_is_continuous_and_increasing_through_the_floor`, `test_let_down_below_the_floor_continues_the_one_above_it` |
| Cold regulated line | Nodes between ambient and inlet T | consistency | Let-down changes T by < 0.5 K (+0.1 K obtained); was node at −22.7 °C and +16.6 K | `tests/test_black_oil_range.py::test_cold_regulated_line_has_physical_temperatures` |
| Black-oil enthalpy inversion | Bracket root | consistency | Fixed v1.27.0: returned 61.64 °C instead of 59.15 °C | `tests/test_regulator.py::test_black_oil_enthalpy_inversion_returns_the_root` |
| Hot-oil pipeline T_out vs NTU (tutorial black-oil 07, 20 km, U 2) | T(L) = T_amb + (T_in − T_amb)·exp(−U·π·D·L/(ṁ·Cp)) | analytical | "Actual error is ~0.27 °C; 2 °C tolerance covers variable-Cp effects" (run with `pipe_enthalpy_terms=False`) | `tests/test_tutorial_suite.py::NonIsothermalBlackOilBenchmarks::test_temperature_profile_matches_analytical`; `src/angelica/tutorials/steady_black_oil/07_hot_oil_pipeline/run.py` |
| Hot vs cold oil flow; viscosity ratio | Qualitative | consistency | Flow ratio ~1.35 (threshold 1.20); viscosity ratio ~70× (threshold 10×) | `tests/test_tutorial_suite.py::NonIsothermalBlackOilBenchmarks::test_hot_oil_flows_more_than_cold`, `test_viscosity_decreases_with_temperature` |
| Black-oil standard-volume mixing at junctions | Exact standard-volume conservation | consistency | Enthalpy equals mass average within 0.01 K·cp (old mass-averaged parameters: 2–4 K·cp) | `tests/test_gas_lift.py::test_black_oil_mixing_conserves_standard_volumes_and_enthalpy` |

## 8. Beggs-Brill two-phase correlation and surface tension

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Inclination correction, 6-in. tubing data set (v_sL 3.97 ft/s, v_sg 3.86 ft/s, ρ_L 47.61, ρ_g 5.88 lbm/ft³, μ_L 0.97 cP, μ_g 0.016 cP, σ 8.41 dyn/cm) | Hand calculation from published equations; coefficients "cross-checked against the whitson wiki tables — READ IN A SECONDARY SOURCE" (pyResToolbox also confirmed, ROADMAP v1.21.0) | secondary source / consistency ("a consistency check of the implementation, not a validation against measurements") | λ_L = 0.507, N_LV = 11.87, intermittent, C clamped to 0, H_L(90°) = 0.574, 0.215 psi/ft; H_L(−45°) = 0.304 (all to 3 decimals) | `tests/test_beggs_brill.py::TestInclinationCorrection::*` |
| H_L(90°) = 0.574 vs Brill & Mukherjee (1999) Example 4.2 | "recalled FROM MEMORY ... not counted as validation" | from memory | Agrees | `tests/test_beggs_brill.py::TestInclinationCorrection` (class docstring) |
| Horizontal gradient, hand case (v_m 3 m/s, C_L 0.333, D 0.1016 m) | Hand calculation: intermittent, N_Fr ≈ 9.03, H_L(0) ≈ 0.459, y ≈ 1.58 | analytical (hand calculation) | Regime/holdup/friction ordering checks | `tests/test_beggs_brill.py::TestBeggsBrillGradient::*` |
| Low-GOR limit: BB ≈ Colebrook | Single-phase limit | consistency | Flows agree within 10 % | `tests/test_beggs_brill.py::TestBeggsBrillIntegration::test_single_phase_bb_close_to_colebrook` |
| Payne et al. (1979) holdup factors 0.924 / 0.685 | "Neither factor could be read in a source"; pyResToolbox notes 0.924 "does not appear in the paper" | not validated ("UNVERIFIED") | Option applies factors exactly (12 places) | `src/angelica/closures/beggs_brill.py`; `tests/test_beggs_brill.py::TestInclinationCorrection::test_payne_factors` |
| Surface tension: Baker & Swerdloff (dead/live oil), Hough gas-water | Equations "Brill & Mukherjee ... Appendix B; Beggs ... 1991" as written in module | consistency (formula reproduction) | 6 / 4 / 3 decimal places; floor 1 dyn/cm | `tests/test_beggs_brill.py::TestSurfaceTension::*`; `src/angelica/properties/surface_tension.py` |
| Tubing correlations spread (tutorial black-oil 12 Part D) | — | consistency | Well 1 max P_wh 11.04 (BB) / 11.30 (HB) / 11.96 MPa (Gray); test case peaks 11.43 / 11.80 / 12.69 MPa | `src/angelica/tutorials/steady_black_oil/12_wells_with_tubing/README.md`; `tests/test_well_correlations.py::test_well_with_each_correlation` |
| Calibration: friction multiplier scales ΔP exactly | f·(L/D)·ρv²/2 | analytical | Tolerance 1e-6 (outer/velocity loop convergence); obtained ≤ 5.9e-9 | `tests/test_calibration.py::test_friction_multiplier_scales_the_pressure_drop_exactly` |
| Calibration: network response (two-loop black-oil, BB horizontal) | Monotonicity | consistency | 25.71 → 25.05 → 22.35 kg/s at 1, 1.05, 1.3; holdup ×1.3 gives 31.18 kg/s (more flow, explained) | `tests/test_calibration.py::test_network_response_to_the_multipliers`, `test_beggs_brill_gradient_multipliers` |

## 9. Multiphase report and erosional velocity

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| API RP 14E erosional velocity V_e = C/√ρ_m | API RP 14E (5th ed., 1991, §2.5) formula; hand calculation | analytical (hand calculation) | ρ_m = 50 lbm/ft³: 14.142 ft/s (C = 100), 17.678 ft/s (C = 125) | `tests/test_multiphase_report.py::test_api_14e_erosional_velocity_hand_calculation` |
| Elevation + friction gradients add up to the section ΔP | Sum | consistency | — | `tests/test_multiphase_report.py::test_gradients_add_up_to_the_section_pressure_drop` |

## 10. Compositional EOS, flash and transport properties

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Methane gas viscosity at 300 K, 1/50/100/150 bar | "CoolProp 6.6 reference correlations (methane: Quiñones-Cisneros & Huber)": 11.255 / 12.173 / 13.815 / 16.163 μPa·s | not labelled (reference implementation) | Tolerances 2 %, 2 %, 2 %, 4 % (source not stated). Obtained (CHANGELOG): −0.6 % at 50 bar, +0.1 % at 100 bar, +3.1 % at 150 bar | `tests/test_compositional.py::TestGasViscosityPressureDependence::test_methane_300K_vs_reference` |
| CH₄/C₂H₆ 90/10, 300 K, 100 bar | CoolProp 14.225 μPa·s | not labelled | Tolerance 2 %; obtained −0.3 % (CHANGELOG) | `...::test_methane_ethane_mixture_100bar` |
| Dense supercritical ethane, 350 K, 100 bar | CoolProp (Friend et al.) 24.42 μPa·s | not labelled | Tolerance 5 %; obtained −2.7 % (CHANGELOG) | `...::test_dense_supercritical_ethane_not_placeholder` |
| n-hexane liquid viscosity, 25 °C | thermo's DIPPR-based value 0.2993 mPa·s | regression | Tolerance 2 % | `...::test_liquid_viscosity_unchanged` |
| McAdams two-phase viscosity blend uses gas mass fraction | Independent formula | consistency | 1e-9 | `...::test_two_phase_mcadams_uses_gas_mass_fraction` |
| Tutorial compositional 01 flash (CH₄/C₂H₆ 80/20, 60 °C) | "CoolProp (GERG-2008 type Helmholtz mixture model with its reference viscosity correlations)" | not labelled | ρ 81.9 (80.1) at 100 bar, 14.2 (14.1) at 20 bar; μ 1.46 (1.48), 1.22 (1.22) ×1e-5 Pa·s; Cp 2869 (2880), 2294 (2290) J/(kg·K). No test. | `src/angelica/tutorials/steady_compositional/01_single_pipe/README.md` |
| Lean-gas cp | CoolProp | not labelled | 2 678 vs 2 680 J/(kg·K) (CHANGELOG, energy-assembly fix) | `CHANGELOG.md` |
| Surface tension of pure n-alkanes at 20 °C (Macleod-Sugden, Peneloux) | "CRC Handbook values WRITTEN FROM MEMORY — not read in a source, so this is a regression guard, NOT a validation" | from memory | Pentane 15.37 (−3.9 %), hexane 17.83 (−3.1 %), heptane 19.49 (−3.0 %), octane 20.51 (−5.0 %), decane 23.62 (−0.8 %). 6 % tolerance "set after seeing these results" | `tests/test_compositional.py::TestGasLiquidSurfaceTension::test_pure_n_alkanes_within_6_percent` |
| Three-phase bulk density with trace water | Dry EOS density | consistency | rel 1e-6; mass/volume rel 1e-9 | `tests/test_compositional.py::TestThreePhaseDensity::*` |
| Flash-cache rounding step | Outer-loop density tolerance 1e-4 | consistency | Old 0.1 K step moved density ~5e-4 | `tests/test_compositional.py::TestFlashCacheResolution` |
| EPPR78 k_ij | Bundled JSON database: CH₄–C₂H₆ 0.0058, CH₄–C₃H₈ 0.0189, C₂H₆–C₃H₈ −0.0029 | consistency | 3 decimal places | `tests/test_phase_envelope.py::TestBuildKijMatrix::test_known_pairs_from_eppr78` |

## 11. Compositional Beggs-Brill

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| No-slip density = bulk flash density; C_L = flash phase volumes | Flash | consistency | rel 1e-6 | `tests/test_compositional.py::TestCompositionalBeggsBrillProperties::*` |
| Network ΔP = Beggs-Brill gradient of the flash phases (dz 0, +100, −100 m) | `beggs_brill_gradient` | consistency | rel 2e-3 on 10 bar | `tests/test_compositional.py::TestCompositionalBeggsBrillSolver::test_pressure_drop_is_beggs_brill_of_flash_phases` |
| Roadmap check "black-oil vs compositional with matched phase properties within 1 %" | — | not done | "was not done as such: a black-oil PVT cannot be matched to an EOS fluid without tuning" | `ROADMAP.md` log v1.15.0 |

## 12. Phase envelope

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Pure methane saturation pressure, 130–180 K | "NIST REFPROP (Lemmon et al. 2009, J. Phys. Chem. Ref. Data 38:721-748)"; values 3.6732e5 … 32.8518e5 Pa | not labelled | "systematic over-prediction of ~0.7%"; tolerance 2 % ("conservative bound for PR EOS on Psat"; source not stated) | `tests/test_phase_envelope.py::TestComputePhaseEnvelope::test_pure_methane_saturation_vs_nist` |
| Pure methane critical point | "Tc ≈ 190.6 K, Pc ≈ 4.60 MPa (NIST)" | not labelled | Tolerance 5 K and 10 %; CHANGELOG v1.8.10: "matches NIST exactly via the PR EOS parameters" | `...::test_pure_methane_critical_near_literature` |
| CH₄/C₂H₆/C₃H₈ 0.70/0.20/0.10 critical point | "HYSYS SRK01 (session data, 2026-08-29)": critical ≈ −17 °C, ≈ 83 bar, cricondenbar ≈ 87 bar | not labelled (session data, not in repo) | Test bounds −35 to 5 °C and 70–100 bar; CHANGELOG v1.6.75: "matches HYSYS within ~4 °C and < 1 % in pressure" (before k_ij were added) | `...::test_critical_T_in_expected_range_PR`, `_SRK`, `test_critical_P_in_expected_range_PR`, `_SRK` |
| PR vs SRK critical point | Each other | consistency | Within 5 K and 10 % | `...::test_pr_srk_critical_T_within_5C`, `test_pr_srk_critical_P_within_10pct` |
| Quality lines, rich condensate, absent components, cricondentherm vs flash | Independent TP flashes | consistency | — | `tests/test_phase_envelope.py::TestQualityLineContinuation`, `TestRichCondensateEnvelope`, `TestAbsentComponents` |

## 13. Hydrate formation conditions (vdW-P, Munck constants) and Katz chart

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Lower quadruple points Q1: CH₄ 272.9 K 2.563 MPa; C₂H₆ 273.1 K 0.530 MPa; C₃H₈ 273.1 K 0.172 MPa; CO₂ 273 K 1.256 MPa | "Sloan & Koh, *Clathrate Hydrates of Natural Gases* — SECONDARY SOURCE, seen only through a search-engine summary of a ScienceDirect topic page (the table itself was not read)" | secondary source | Methane 2.751 MPa (+7.3 %, "part of it is the 0.25 K offset"), ethane 0.496 (−6.5 %), propane 0.166 (−3.6 %), CO₂ 1.240 (−1.2 %). Tolerance 11 %: accuracy reported by Antunes et al. (2018, Braz. J. Chem. Eng. 35, 265 — READ IN SOURCE) for single-guest vdW-P + Munck; "not chosen from the results" | `tests/test_hydrate.py::test_lower_quadruple_points` |
| Upper quadruple points Q2: C₂H₆ 287.8 K at 3.390 MPa; C₃H₈ 278.8 K at 0.556 MPa | Same secondary source | secondary source | Ethane 286.98 K (−0.82 K), propane 278.52 K (−0.28 K). Tolerance 1 K ("≈ 11 % in pressure on this part of the curve; the roadmap's criterion") | `tests/test_hydrate.py::test_upper_quadruple_points` |
| Methane solubility in water (Henry term), 283 K and 298 K, 1.2–10.4 MPa, 13 points | Böttger, Kamps & Maurer (2016), Fluid Phase Equilib. 407, 209, NIST ThermoML archive | read in source | Errors −4.1 … +3.5 % (largest absolute +0.0040 mol/kg at 10.36 MPa); all inside the expanded (95 %) uncertainty of each point (the tolerance) | `tests/test_hydrate.py::test_methane_solubility_against_measurements` |
| Pure methane hydrate, 274–289 K, 2.9–14.2 MPa (12 points) | Semenov et al. (2024), Data in Brief 53, 110138, Table 1, system 1 — READ IN SOURCE (Europe PMC, PMC10877679) | read in source | Model colder by 0.60 K at 2.9 MPa growing to 0.91 K at 14.2 MPa. Guard 1.0 K "set after the results, not an accuracy claim" (test asserts −1.0 < e < −0.55) | `tests/test_hydrate.py::test_methane_against_semenov_2024` |
| Hydrate T at 10–100 MPa: methane 40 MPa 297 K; CO₂ 10 MPa 283.6 K, 40 MPa 286.2 K | "Sloan & Koh values as quoted by Algaba et al. (2024), arXiv:2408.01819, Table III" | secondary source | Methane 40 MPa −1.56 K; CO₂ 10 MPa +0.36 K, 40 MPa −1.35 K. Tolerance 2 K "a regression guard chosen after seeing the results ... NOT a validation" | `tests/test_hydrate.py::test_high_pressure_hydrate_temperature_guard` |
| CO₂ hydrate at 100 MPa | Same source, 289.7 K | secondary source | Model 280.8 K (xfail strict, "Known limitation") | `tests/test_hydrate.py::test_co2_hydrate_at_100_mpa_known_limitation` |
| Katz chart vs Towler-Mokhatab at 500 psia, γ 0.6–1.0 | Janna (2015), Tables 4.1b–4.5b (read in the dissertation) — "CROSS-CHECK, not validation" | secondary source | Differences −1.45, −0.13, −0.77, −1.11, −1.81 °F; tolerance 2 °F "regression guard chosen after the results"; at 2000 psia 2.4 °F (γ 0.6) to 8.5 °F (γ 1.0) | `tests/test_hydrate.py::test_katz_chart_against_towler_mokhatab` |
| Dissolved gas lowers water activity (CO₂ vs CH₄) | Ordering | consistency ("not validation") | 0.97 < x_w(CO₂) < 0.985; x_w(CH₄) > 0.999 | `tests/test_hydrate.py::test_dissolved_gas_lowers_water_activity` |
| Propane promotes sII; curve monotonic; inverse consistent | Ordering / inverse | consistency | Inverse within 0.02 K; NG pressure < 0.5 × methane | `tests/test_hydrate.py::test_propane_promotes_structure_ii_in_natural_gas`, `test_curve_monotonic_and_inverse_consistent` |
| Hydrate margin per pipe = stand-alone function | Same function | consistency | 1e-6 K (tutorial 08); two-loop networks | `tests/test_hydrate.py::test_compositional_results_report_hydrate_margin`; `tests/test_two_loop_networks.py::test_hydrate_margin_per_pipe` |
| Langmuir constants | Munck et al. (1988) "as tabulated ... by Antunes et al. (2018) ... Table 2" ("verified verbatim ... via pdftotext of the open PDF", ROADMAP v1.19.0) | read in source (parameters) | — | `src/angelica/properties/hydrate.py` |

## 14. Hydrate as a phase (amount formed)

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Methane hydration number | Sum, Burruss & Sloan (1997), J. Phys. Chem. B 101, 7371: Raman 6.04 ± 0.03 — "Read in a secondary source (a search-engine excerpt of the USGS abstract ... HTTP 403)"; measurement conditions not read | secondary source — "COMPARED, NOT VALIDATED (a recorded discrepancy)" | Model 5.895 at 2 °C and 70 bar (−2.4 %, outside ±0.03); test pins the value and the sI bound 46/8 = 5.75 | `tests/test_hydrate_flash.py::test_pure_methane_is_limited_by_the_gas` |
| Incipient boundary = vdW-P curve | vdW-P curve | consistency | Amount zero at the curve, grows continuously (0.198 mol/K, STATUS) | `tests/test_hydrate_flash.py::test_incipient_boundary_is_the_vdwp_curve` |
| Mass balance of water and each former | Sums | consistency | Tolerance 1e-12 | `tests/test_hydrate_flash.py` (module docstring) |
| Leftover gas at its own hydrate condition / hydrate-vapour curve | Same functions | consistency | 1e-5 K (formation solve stops at 1e-4; obtained ≤ 1.5e-6) | `tests/test_hydrate_flash.py::test_equilibrium_with_free_water_left`, `test_water_limited_gas_ends_on_the_hydrate_vapour_curve` |
| Equilibrium case: 95 % water, 12 °C | — | consistency | 0.0888 mol water in hydrate, hydration number 6.156, propane 3.0 % → 0.27 % | `tests/test_hydrate_flash.py::test_equilibrium_with_free_water_left` |
| Oil-gas-water feed, 8 °C, 70 bar | — | consistency | 0.499977 mol water in hydrate, 2.273e-5 vapour, G/RT = −7.7e-13; tolerance 1e-9 | `tests/test_hydrate_flash.py::test_oil_gas_water_feed` |
| Little gas consumed | — | consistency | ξ/z_gas = 6.0677 at 1e-5 and 1e-3; leftover ≤ 6e-6 of feed gas | `tests/test_hydrate_flash.py::test_little_gas_is_consumed_to_the_stated_accuracy` |

## 15. Water content of gas and water dew point

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Water saturation pressure (Wagner form), 1–91 °C | NIST Chemistry WebBook, water, saturation properties (IAPWS-95) (retrieved 2026-11-27) | read in source | Deviations −7.2e-5 … +4.9e-5; tolerance 1e-4 "a regression guard set after the results (no published accuracy of this Wagner form was read)" | `tests/test_water_content.py::test_saturation_pressure_against_nist` |
| Water dew point inverts the curve | NIST pressures → NIST temperatures | read in source | Within 0.002 K | `tests/test_water_content.py::test_dew_point_inverts_the_saturation_curve` |
| Saturated gas reproduces the water-saturated hydrate curve | Exact | consistency | — | `tests/test_water_content.py::test_saturated_gas_reproduces_the_water_saturated_curve` |
| Undersaturated gas hydrate T between dew point and saturated curve (70 bar) | Ordering | consistency | Dew 15 → 15.43 °C, 5 → 7.22, 0 → 3.06 (saturated curve 17.03) | `tests/test_water_content.py::test_undersaturated_gas_forms_hydrate_between_its_dew_point_and_the_saturated_curve` |
| Compositional pipe with 60 ppm water | Direct call | consistency | 0.420 °C at the 80 bar end (subcooling −11.58 K); tolerance 1e-9 K | `tests/test_water_content.py::test_compositional_pipe_uses_the_stream_water_content` |
| GUI curve at water content (60 ppm) | Solver rule | consistency | 77.26 bar at 0 °C, 116.45 at 5 °C; at 10 °C saturated point 25.51 bar | `tests/test_water_content.py::test_curve_at_water_content_follows_the_solver_rule` |

## 16. Hydrate inhibitors (methanol, MEG)

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Methanol depression ΔT_h, 5–50 mass %, 3–13 MPa (new samples) | [S24] Semenov, Mendgaziev, Istomin, Sergeeva, Vinokurov, Gong, Li, Stoporev (2024), Data in Brief 53, 110138, Tables 1 and 10 (PMC10877679) — READ IN SOURCE; "new samples FROM THE SAME LABORATORY" as the fit data (not independent) | read in source (same lab) | |ΔT_h error| ≤ 0.07 K at 5 %, 0.19 at 10 %, 0.39 at 20 %, 0.78 at 30 %, 0.98 at 40 %, 0.94 at 44.7 %, 1.39 K at 50 %. Under-predicted at 3 MPa (−1.39 K, safe side); over-predicted above ~7 MPa (+0.3 to +1.02 K at 9.7–12.8 MPa, 20–50 %, unsafe side). Tolerance "a regression guard set after the results (≈ 1.2× the obtained), not a claimed accuracy; [S24] gives no uncertainty for ΔT_h" | `tests/test_inhibitor.py::test_methanol_against_independent_samples` |
| Methanol vs its fitting data | [S23] Semenov, Gong, Medvedev, Stoporev, Istomin, Vinokurov, Li (2023), Data in Brief 46, 108892, Table 1 (PMC9876827) — READ IN SOURCE | consistency (fitted data) | |error| ≤ 0.41 K (fit at 6.5 MPa; grows toward 5 and 8 MPa). Guard 0.5 K | `tests/test_inhibitor.py::test_methanol_reproduces_its_fitting_data` |
| MEG 20 % and 40 % (left out of the fit, same dataset) | [S23] | consistency ("a check within the same dataset, not independent validation") | Max 0.13 and 0.41 K; guards ≈ 1.2× obtained | `tests/test_inhibitor.py::test_meg_against_its_dataset` |
| Ice points of methanol solutions | [S24] Table 10 | read in source (comparison, recorded) | Within 0.07 K up to 10 % (module docstring: "within 0.1 K"); too high by 1.0, 3.0, 5.3, 8.9 K at 20, 30, 40, 50 %; not used as a limit | `tests/test_inhibitor.py::test_ice_points_of_methanol_solutions`; `src/angelica/properties/inhibitor.py` |
| Absolute hydrate T with inhibitor at ≥ 9.7 MPa | [S24] | read in source (recorded, not a tolerance claim) | Too cold by up to 1.89 K (50 %, 12.76 MPa): pure-water offset (−0.9 K) + depression over-prediction (+1.0 K) | `tests/test_inhibitor.py::test_absolute_error_against_the_2024_data_at_high_pressure` |
| Ice-point Δc_p constants | Liquid water 76.014 J/(mol·K) at 273.16 K (NIST WebBook, IAPWS-95, read 2026-11-27); ice Ih 2096.71 J/(kg·K) (IAPWS R10-06(2009), Table 6, read 2026-11-27) | read in source (parameters) | — | `src/angelica/properties/inhibitor.py` |
| Solver hydrate checks with inhibitor | Direct call | consistency | Tolerance 1e-9 K; black-oil depressions at 70 bar, 30 % methanol: 16.36 and 14.91 K | `tests/test_inhibitor.py::test_compositional_solver_uses_the_inhibitor`, `test_black_oil_solver_uses_the_inhibitor`, `test_black_oil_depression_is_the_smaller_of_two_gases` |
| Partition with inhibitor | — | consistency | 0.1069 mol water in hydrate, methanol 20 % → 24.1 %; no free water case 4.3755e-6 mol | `tests/test_inhibitor.py::test_hydrate_partition_with_an_inhibitor`, `test_hydrate_partition_without_free_water_ignores_the_inhibitor` |

## 17. Heat transfer and thermal solver

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Oil pipeline through an icy lake (2 m/s, D 30 cm, L 200 m, h 16.3 W/m²K, ṁ 125.6 kg/s) | "Cengel, Y.A. & Ghajar, A.J., *Heat and Mass Transfer*, 5th ed., McGraw-Hill, 2015, Example 8-3, pp. 486-488": T_e = 19.74 °C | not labelled (textbook citation with pages; read status not stated) | Within 0.05 °C (tolerance source not stated) | `tests/test_non_isothermal.py::test_cengel_example_8_3_oil_pipeline_through_lake` |
| Single pipe with heat loss, upwind and hybrid/power-law schemes | Exponential-decay (NTU) solution | analytical | Upwind within 1 °C; other schemes within 1.5 °C | `tests/test_non_isothermal.py::test_outlet_temperature_upwind`, `test_outlet_temperature_alternative_schemes` |
| Mesh convergence to NTU | NTU | analytical | Error decreasing; < 0.5 °C at n = 100 | `tests/test_non_isothermal.py::test_mesh_convergence_approaches_ntu_analytical` |
| Two pipes in series | Chained NTU | analytical | Within 1 °C | `tests/test_non_isothermal.py::test_series_pipes_chained_ntu` |
| Adiabatic mixing junction, inline heater fixed flow, symmetric loops | Exact mixing T; ΔT = Q/(ṁcp); symmetric NTU | analytical | 0.1 K (mixing), 0.2 K (heater), 0.1 K (adiabatic loop), 0.3 K (heat-loss loop); flows 1e-4–1e-6 | `tests/test_non_isothermal.py::test_adiabatic_mixing_junction_matches_exact_mixing_temperature`, `test_inline_heater_fixed_flow_matches_exact_delta_t`, `test_symmetric_adiabatic_loop_stays_isothermal_and_splits_flow_evenly`, `test_symmetric_heat_loss_loop_matches_exact_branch_temperature` |
| Global energy conservation, adiabatic network | ΣH_in = ΣH_out | consistency | Within 1 % | `tests/test_non_isothermal.py::test_global_energy_conservation_adiabatic_network` |
| Segmented wall heat loss, pure advection upwind | Discrete solution θ_out = (1 + NTU/(n−1))^−(n−1) | analytical | 1e-4; many segments → exp(−NTU) within 0.1 K; energy error < 0.1 % | `tests/test_wall_heat_loss.py::test_segmented_heat_loss_matches_discrete_upwind_solution`, `test_segmented_heat_loss_converges_to_analytical`, `test_global_energy_balance_closes` |
| Compressible single pipe heat loss (1 segment) | T_amb + (T_in − T_amb)·exp(−NTU) | analytical | 1e-6 K ("confirms the analytical path is active, not that ṁ is exact") | `tests/test_compressible.py::QuantitativeBenchmarks::test_single_pipe_heat_loss_ntu` |
| Tutorial non-isothermal 01 (water, 500 m, D 25 mm) | NTU, ~38.8 °C | analytical | "exact regardless of the number of pipe segments" | `src/angelica/tutorials/steady_non_isothermal_incompressible/01_single_pipe_heat_loss/README.md` |
| Tutorial non-isothermal 04 inline heater | ΔT = Q/(ṁ·cp) ≈ 0.86 K | analytical | Within 0.5 K | `tests/test_tutorial_suite.py::TutorialSuiteSmokeTests::test_inline_heater_converges`; README 04 |
| Cylindrical conduction through layers | R' = ln(r_o/r_i)/(2πk), exact radial solution | analytical | Single layer to 1e-12; thin-wall limit = flat-wall sum | `tests/test_wall_heat.py::test_single_cylindrical_layer_is_exact`, `test_thin_wall_limit_is_the_flat_sum` |
| Buried pipe shape factor S/L = 2π/arccosh(2H/D) | Method of images | analytical | Same S to 1e-12 | `tests/test_wall_heat.py::test_buried_pipe_by_the_method_of_images` |
| Buried pipe, D 0.5 m, depth 2 m, L 5 m | EES heat-transfer library, function SF_2: 11.33 m — READ IN A SECONDARY SOURCE, citing Incropera & DeWitt, 5th ed., pp. 193–194 | secondary source | Angelica 11.347 m (+0.15 %); library value matches ln(4z/D) form 11.331 m. Tolerance 0.2 %: "the difference between the two forms of that table at this depth" | `tests/test_wall_heat.py::test_buried_pipe_against_the_ees_example` |
| Inner film: Gnielinski (Re ≥ 3000) and Nu = 3.66 laminar | READ IN A SECONDARY SOURCE (Wikipedia "Nusselt number", citing Incropera & DeWitt p. 490) | secondary source / consistency; "NOT validated against measurements" | Agrees with Gnielinski's simplified power laws within 13.1 % and Dittus–Boelter within 10.4 % over Re 1e4–1e6 (guards 15 % and 12 %) | `tests/test_wall_heat.py::test_inner_film_correlation` |
| Computed film in solvers | Energy balance | consistency | < 0.01 % | `tests/test_wall_heat.py::test_solver_uses_the_computed_film_and_reports_it`, `test_viscous_oil_with_computed_film_converges` |

## 18. Energy balance and first law (gas and compositional pipes, valves)

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Adiabatic PR-methane pipes (80 → 50 bar, level and 800 m rise, 1 and 20 segments): h_in = h_out + g·Δz | Exact first law with the fluid's own enthalpy | consistency | Errors (v1.18.0) ≤ 0.001 K on 11–14 K of cooling; tolerance 0.01 K = outer-loop temperature tolerance; residual < ṁ·cp·0.01 K | `tests/test_compressible.py::PipeEnthalpyBalanceTests::*` |
| Ideal-gas horizontal pipe has no JT | Ideal gas h(T) | analytical | — | `tests/test_compressible.py::PipeEnthalpyBalanceTests::test_ideal_gas_horizontal_pipe_has_no_jt` |
| Valve (K = 50) isenthalpic, methane PR 70 bar / 40 °C | h_in = h_out | consistency | Relative error < 1e-3 (outer loop 0.01 K); energy error < 1e-3 % | `tests/test_compressible.py::QuantitativeBenchmarks::test_fitting_isenthalpic_enthalpy_conserved` |
| Adiabatic compositional pipes (single-phase 100 → 20 bar; rich gas condensing to VF 0.85–0.95 at 80 → 40 bar, level and 500 m rise) | "thermo's PH flash" — "CONSISTENCY check against an independent calculation" | consistency | ≤ 0.001 K; tolerance 0.01 K; energy error < 0.01 %; potential energy term rel 1e-9 | `tests/test_compositional.py::TestCompositionalPipeFirstLaw::*` |
| Compositional valve isenthalpic and carries composition | PH flash of inlet state | consistency | — | `tests/test_compositional.py::TestCompositionalPipeFirstLaw::test_valve_is_isenthalpic_and_carries_composition` |
| Fitting adiabatic in incompressible solver | T_out = T_in | analytical | 1e-4 K | `tests/test_non_isothermal.py::test_fitting_is_adiabatic_incompressible` |
| Global mass and energy balance of every tutorial GUI file | Boundary fluxes, heat sources, wall loss (independent of assembly), shaft work, JT of valves | consistency | < 0.01 % each | `tests/test_tutorials.py::test_tutorial_gui_file` |
| Two-loop thermal runs | Energy balance | consistency | Obtained ≤ 2.77e-4 % (roadmap figure 0.01 %) | `tests/test_two_loop_networks.py::test_thermal_runs_close_the_energy_balance` |

## 19. Compressor

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Ideal gas, constant cp: T2s = T1·r^(R/cp), W_s = cp·(T2s − T1) | Exact | analytical | Errors (20 steps): −3.9e-6 (ratio 1.5) to −8.8e-5 (ratio 6) in W_s; −0.0001 to −0.013 K in T2s. Tolerance 5e-4 of W_s ("a numerical bound, not a property of data") | `tests/test_compressor.py::test_ideal_gas_isentrope_is_analytical` |
| Isentrope vs EOS entropy (PR methane; compositional natural gas with thermo flash entropy) | s(P2, T2s) = s(P1, T1), entropy computed independently | analytical / against the EOS | T2s·|Δs|/W_s: 6.2e-5 … 1.3e-4 (PR, ratio 1.5 … 6); 4.1e-6 … 4.7e-5 (compositional). Tolerance 5e-4 of W_s | `tests/test_compressor.py::test_isentrope_against_eos_entropy` |
| Network consistency (gas, black-oil, compositional) | Specification at node pressures; power = ṁ·Δh | consistency | Ratio/rise/power within 1e-4 (outer-loop criterion); mass within 2e-3; energy balance ≤ 0.01 % (STATUS: power 3e-4, balance ≤ 1.5e-3 %) | `tests/test_compressor.py::test_network_compressor`, `test_polytropic_in_a_network` |
| Small nozzle meets specification | Node pressures | consistency | Before review missed by up to 4.4e-3 at D = 0.05 m | `tests/test_compressor.py::test_small_compressor_nozzle_still_meets_the_specification` |

## 20. Heater / cooler

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Water heater duty (heating 20 → 60 °C, cooling 80 → 10 °C) | NIST Chemistry WebBook, water, isobar 0.5 MPa (retrieved 2026-11-27): h(10 °C) 42.5077, h(20 °C) 84.3824, h(60 °C) 251.583, h(80 °C) 335.373 kJ/kg | read in source | +0.048 % (heating), −0.044 % (cooling). Tolerance 0.51 %: largest deviation of NIST cp from 4182 between 5 and 80 °C, "a property of the data, not of these results" | `tests/test_heater_outlet.py::test_water_duty_against_nist` |
| Outlet = set point; duty = ṁ·(h_out − h_in) | Fluid's own enthalpy | consistency | Node after heater ≤ 0.0040 K (this network; ≈ 1e-3 × T °C general); duty 7.4e-6; balance ≤ 1.9e-4 % (STATUS) | `tests/test_heater_outlet.py::test_outlet_temperature_set_point`, `test_duty_limited` |
| Fixed duty in a gas (city gate, 500 kW) | h(P_out, T_out) = h_in + Q/ṁ | consistency | Balance now 4.2e-4 % (was 0.39 %; 7.1 % methane, 5.9 % condensate in the test network) | `tests/test_heater_outlet.py::test_fixed_duty_in_a_gas_uses_the_real_enthalpy` |
| Tutorial black-oil 16 heater station | Energy balance | consistency | Residual 0.034 kW of 2 668 kW (2.9e-4 %) | `src/angelica/tutorials/steady_black_oil/16_heater_station/README.md` |

## 21. Chokes

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Orifice K-factor | K = (D_pipe/D_orifice)⁴/Cd² ("standard orifice equation"); default Cd 0.62 "(ISO 5167 sharp-edged)" | analytical | — | `tests/test_choke.py::TestChokeComponent::test_k_factor_formula`, `test_default_cd_is_sharp_edged` |
| Sachdeva pure-gas limit | y_c = (2/(k+1))^(k/(k−1)); G(y) peaks at y_c | analytical | — | `tests/test_choke.py::TestSachdevaLimits::test_pure_gas_matches_isentropic_ratio`, `test_pure_gas_mass_flux_maximum_at_critical_ratio` |
| Sachdeva pure-liquid limit | G = Cd·√(2·ρ_L·ΔP) | analytical | — | `tests/test_choke.py::TestSachdevaLimits::test_pure_liquid_mass_flux_is_orifice_equation` |
| Solver rate = Cd·A·G at solved upstream pressure (critical) | Closed form | consistency | rel 1e-6 | `tests/test_choke.py::TestCriticalChokeIntegration::test_choke_flow_matches_closed_form` |
| Orifice model vs Sachdeva at critical flow | — | consistency | Orifice > 1.1 × Sachdeva | `tests/test_choke.py::TestCriticalChokeIntegration::test_orifice_model_overpredicts_critical_choke` |

## 22. Separators

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Black-oil gas mass fraction = in-situ free-gas mass fraction | Black-oil PVT | consistency (identity) | < 1e-12 | `tests/test_separator.py::test_black_oil_gas_fraction_is_the_in_situ_free_gas` |
| Compositional split | Independent `_hc_phase_split` + immiscible water; species conserved | consistency | < 1e-9 | `tests/test_separator.py::test_compositional_split_conserves_species_and_matches_phase_split` |
| Enthalpy additivity of the split | Feed enthalpy | consistency | Tolerance 0.0015 K·cp black-oil (obtained 2e-4 to 1.5e-3 K at 0.5–8 MPa); 1e-5 K·cp compositional (obtained 6e-6 K) | `tests/test_separator.py::test_black_oil_split_conserves_enthalpy`, `test_compositional_split_conserves_enthalpy` |
| Boundary separator node temperature | Same separator 1 m away | consistency | Tolerance 0.01 K (obtained 2e-4 K) | `tests/test_separator.py::test_boundary_node_with_throughflow_is_at_the_mixing_temperature` |
| Stage-separation optimum (tutorial compositional 10) | — | consistency (model result) | Optimum ≈ 18 bar one stage, ≈ 30 bar two stages; LP stage gain +2.5 / +4.9 / +4.1 % | `src/angelica/tutorials/steady_compositional/10_separator_train/README.md` |

## 23. Pressure regulators

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Controlling: flow = downstream network solved alone at the set point | Independent solve | consistency | Obtained 1.2e-7; tolerance solver nodal tolerance | `tests/test_regulator.py::test_reducing_regulator_holds_its_outlet` |
| Wide open = check valve of same K | Independent network | consistency | Obtained 1e-7 | `tests/test_regulator.py::test_reducing_regulator_wide_open_is_a_check_valve` |
| Isenthalpic let-down (64 → 20 bar, ~27 K cooling) | Fluid's own h(P, T) | consistency | Tolerance cp·0.01 K | `tests/test_regulator.py::test_compositional_reducing_regulator_is_isenthalpic` |
| Sweeps monotonic | — | consistency | Tolerance 1e-3 (flow may vary 1e-5) | `tests/test_regulator.py::test_regulator_sweeps_converge_and_are_monotonic` |
| BPR then reducer | Downstream pipe solved alone | consistency | Obtained 1.3e-7 | `tests/test_regulator.py::test_back_pressure_then_reducing_regulator_solution` |
| Regulator pairs off-suite sweep | — | consistency | 725 solves converged and consistent (ROADMAP v1.27.1) | `ROADMAP.md` log v1.27.1 |

## 24. Wells

### 24.1 Inflow (IPR)

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Linear, composite Vogel, Fetkovich, back-pressure formulas; inverses; derivatives | Formulas (Vogel 1968; Brown 1984; Fetkovich 1973; Rawlins & Schellhardt 1936 as cited in module) | consistency | Vogel q_max = J·P_r/1.8; continuity at P_b rel 1e-6; derivatives vs FD rel 1e-6 / 1e-5; formula rel 1e-12; n = 1 productivity 2·C·P_r rel 1e-3 | `tests/test_well_inflow.py::TestIprFormulas`, `TestFetkovichAndBackPressure` |
| Solver reproduces IPR | Formula | consistency | rel 1e-6 | `tests/test_well_inflow.py::test_solver_reproduces_ipr`, `test_solver_reproduces_fetkovich_and_backpressure` |
| Against a multirate-test worked example | Fetkovich 1973, Beggs's examples | not validated | "the maintainer could supply one" | `ROADMAP.md` log v1.23.0 |

### 24.2 Well tubing (deliverability)

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Dead-oil tubing vs independent integration | Hydrostatic + Haaland friction | consistency | Tolerance rel 5e-3 ("to 0.5 %", CHANGELOG) | `tests/test_well_model.py::test_single_phase_tubing_matches_hydrostatic_plus_darcy` |
| Hagedorn-Brown (modified) and Gray hand calculations | Equations "from the whitson wiki and chart-fit coefficients from pyResToolbox — both read in secondary sources"; HB coefficients also in rNodal (credits Guo, Lyons & Ghalambor 2007; "The book itself was not read") | secondary source / consistency; "NOT VALIDATED against measured well data" | rel 1e-12 on holdup and gradient; dry-gas limit rel 1e-9; same static column at zero rate rel 1e-6 | `tests/test_well_correlations.py::test_hagedorn_brown_hand_calculation`, `test_hagedorn_brown_bubble_flow_griffith`, `test_gray_hand_calculation`, `test_single_phase_gas_limit`, `test_well_with_each_correlation` |
| Bistable / dead / unstable well states | Logic | consistency | e.g. weak well static 3.04 MPa, maximum 5.49 MPa at 0.94 kg/s, manifold 5.42 MPa; producing case 7.876 / 1.946 kg/s | `tests/test_well_model.py::test_bistable_well_*`, `test_well_with_no_steady_state_is_reported_unstable` |

### 24.3 Tubing temperature (Ramey)

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Incompressible static column vs Ramey closed form | Own derivation; f(t) and relaxation distance "READ IN A SECONDARY SOURCE (OpenGeoSys benchmark documentation 'Wellbore Heat Transport - EUBHE', which reproduces Ramey 1962)" | analytical ("a verification of the integrator against the analytical solution") | Error < 1e-9 K; tolerance 0.001 K (numerical) | `tests/test_well_thermal.py::test_incompressible_static_column_reproduces_ramey` |
| Black-oil tubing first law | ṁ·[h(P_res,T_res) − h(P_wh,T_wh)] = Q_formation + ṁ·g·Δz | consistency | Tolerance ṁ·cp·0.01 K (50–150 W) | `tests/test_well_thermal.py::test_black_oil_tubing_satisfies_the_first_law` |
| Network with a well | Energy balance | consistency | Residual 2e-4 kW against ṁ·cp·0.01 K ≈ 0.16 kW | `tests/test_well_thermal.py::test_well_delivers_its_wellhead_temperature_to_the_network` |
| Flowing temperatures vs field data | Sagar, Doty & Schmidt 1991 "not accessible" | not validated | — | `tests/test_well_thermal.py` (module docstring) |

### 24.4 Compositional wells

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| PVT table vs direct EOS flash march (gas and condensate wells, 10–90 % open flow) | Direct flash at every evaluation | consistency | P_wh within 0.05 %, T_wh within 0.01 K. Tolerances 0.2 % and 0.1 K: "well below the ±0.3 % P_wh scatter already caused by the tubing segmentation" | `tests/test_well_compositional.py::test_pvt_table_matches_direct_flash` |
| Dry-gas tubing vs independent integration | Hydrostatic + Haaland with EOS density/viscosity | consistency | 0.045 %; tolerance 0.5 % (Churchill vs Haaland) | `tests/test_well_compositional.py::test_dry_gas_tubing_matches_single_phase_integration` |
| Mass per stock-tank m³ vs black-oil formula ρ_o + GOR·ρ_g | Formula | consistency | 0.1 %; tolerance 0.5 % (black-oil ideal gas at sc, Z ≈ 0.997) | `tests/test_well_compositional.py::test_standard_volumes_agree_with_black_oil_formula` |
| Heat from table march vs EOS enthalpies | Independent | consistency | 0.07 %; tolerance 0.5 % | `tests/test_well_compositional.py::test_gas_well_in_a_network` |
| IPR-only link reproduces back-pressure law | Formula | consistency | — | `tests/test_well_compositional.py::test_ipr_only_link_with_compositional_fluid` |

### 24.5 Gas lift

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Continuity as lift rate → 0 | Natural well | consistency | 1 Sm³/d changes P_wh by 0.003 % | `tests/test_gas_lift.py::test_lift_rate_to_zero_recovers_the_natural_well` |
| First law of the two-section march | Energy | consistency | Residual 0.01 W of 500 kW | `tests/test_gas_lift.py::test_two_section_march_satisfies_the_first_law` |
| Performance curve shape (rise, peak, fall) | Qualitative | consistency | — | `tests/test_gas_lift.py::test_lift_revives_a_dead_well_and_has_an_optimum` |
| Against published gas-lift data | "Brown's gradient curves, Economides' examples not accessible" | not validated | — | `tests/test_gas_lift.py` (module docstring) |

## 25. Nodal analysis and sensitivity studies

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Nodal operating point = full network solve (black-oil tutorial 12) | Network solve | consistency | Well 1: 6.6933 vs 6.6932 kg/s (1.5e-5), 4.04543 vs 4.04548 MPa; well 2: 2.2189 vs 2.2190 kg/s. Tolerance 1e-3 (bisection 1e-4 of rate; hydraulic criterion 1e-3 per node) | `tests/test_studies.py::test_nodal_operating_point_is_the_network_solution` |
| Compositional tutorial 09 | Network solve | consistency | 6.2268 = 6.2268 kg/s; 6.97104 vs 6.97114 MPa (1.4e-5) | `tests/test_studies.py::test_nodal_compositional_gas_condensate_well` |
| Review cases (PVT correlation, shared reservoir, isothermal solver, no reservoir thermal BC, shared wellhead) | Network solve | consistency | First version errors 2 % to 127 % (shared reservoir +2.7 %, Standing substitution +3.9 %); now ≤ 3.8e-5 or refused | `tests/test_studies.py::test_nodal_with_a_shared_reservoir_node`, `test_nodal_keeps_the_case_pvt_correlation`, `test_nodal_with_an_isothermal_solver`, `test_nodal_refuses_a_shared_wellhead`; `STATUS.md` v1.32.0 |
| Against a published nodal-analysis example | Economides et al., Beggs "not in the repository" | not validated | — | `tests/test_studies.py` (module docstring) |

## 26. Line sizing and per-connection profiles

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Line sizing picks smallest diameter meeting criteria; segmented line sized as a whole | Inverse of solved case | consistency | — | `tests/test_calibration.py::test_line_sizing_selects_the_smallest_diameter_that_meets_the_criteria`, `test_line_sizing_a_segmented_line` |
| Profile report: segments chain, lengths/rises/drops add up | Sums of the same numbers | consistency | Tolerance 1e-9 | `tests/test_profiles.py::*` |

## 27. Units and standard-volume rates

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Conversion factors (barrel, gallon, ft³, ft, in, lb, psi, cP) | NIST SP 811 (2008 edition), Appendix B.8, nist.gov (retrieved 2026-11-27), 7 significant figures | read in source | Every factor within half a unit in the 7th printed digit (the tolerance: the table's rounding) | `tests/test_units.py::test_conversion_factors_against_nist_sp811` |
| Black-oil standard rates by definition (q_o = ṁ/(ρ_o,sc + GOR·ρ_g,sc + WOR·ρ_w,sc)) | Definition, recomputed independently | consistency | First version lost Standing's residual Rs: −0.7 % gas in tutorial 14 (fixed) | `tests/test_units.py::test_black_oil_rates_by_definition` |
| Stream feeding a separator carries the separator's rates (tutorial 14) | Separator rates | consistency | Oil 705.49 = 705.49 Sm³/d (≤ 1e-6); water 102.1589 vs 102.1587 m³/d (1.6e-6); gas through separator within 0.11 % | `tests/test_units.py::test_black_oil_standard_rates_match_the_separator` |
| Compositional water rate | All stream water as liquid at sc | consistency | First version dropped 24 % (fixed) | `tests/test_units.py::test_compositional_water_counts_all_the_water` |

## 28. Isothermal mode

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Every node and well tubing at the given temperature; ignored thermal devices named; out-of-range runs end converged or diagnosed | Itself ("Internal consistency only: an isothermal solve is a modelling choice, with nothing to compare against but itself") | consistency | — | `tests/test_isothermal.py::*`; `tests/test_two_loop_networks.py::test_isothermal_runs_hold_every_node` |

## 29. Tutorial regression (golden outputs)

| What | Reference | Label | Result | Where |
|---|---|---|---|---|
| Every tutorial `run.py` output | Stored reference `tests/tutorial_reference/scripts` | regression | Every printed number within one unit of its last printed digit | `tests/test_tutorials.py::test_tutorial_script` |
| Every tutorial `.gui.json` | Stored reference `tests/tutorial_reference/gui` | regression | Node P, T and mass flows match; mass and energy balance < 0.01 % | `tests/test_tutorials.py::test_tutorial_gui_file` |
| Tutorial-suite numeric checks (e.g. compressible T01 pressures 667 310 Pa ± 500, flows 2.8206 ± 0.01 kg/s; black-oil manifold 3 723 268 Pa ± 50 000; non-isothermal T05 trunk 9.92 ± 0.1 kg/s) | Stored values | regression | — | `tests/test_tutorial_suite.py::*` |

---

## Not validated and known deviations

Everything the repository states as not validated, unverified, or deviating, with its numbers.

### Hydrates and flow assurance
- Natural-gas mixtures (hydrate curve): "NOT validated" (`src/angelica/properties/hydrate.py`, `tests/test_hydrate.py`); roadmap criterion "two natural-gas mixtures within 1 K" open (`ROADMAP.md` 2.1).
- Pure methane: "NOT VALIDATED (no measured data accessible): pure methane between 285 K and 295 K" (`tests/test_hydrate.py`); module: "methane above 289 K are NOT validated".
- Methane hydrate warm-end deviation: model 0.60 to 0.91 K too cold, 274–289 K, 2.9–14.2 MPa (Semenov 2024).
- Methane Q1 +7.3 % in pressure; quadruple points overall "within 7.3 % in pressure at Q1 and 0.8 K at Q2".
- Dissolved-gas term raises methane formation pressure by 0.7 % at 273 K, 2 % at 285 K, 4 % at 290 K, 7 % at 295 K; "not the cause of the model running high at the warm end"; whether Munck et al. fitted with it "could not be established".
- CO₂ above about 40 MPa: predicted hydrate temperature falls; 280.8 K at 100 MPa against 289.7 K (~9 K low). A Poynting factor (v̄ = 32 cm³/mol) would put methane solubility 5–10 % below measurements, so not added.
- High-pressure secondary-source comparisons: methane 40 MPa −1.56 K (295.4 K vs 297 K, "−1.6 K"); CO₂ 10 MPa +0.36 K, 40 MPa −1.35 K — regression guard 2 K, "NOT a validation".
- Katz chart: not validated against measurements ("chart accuracy per literature is a few K", ROADMAP v1.20.0); vs Towler-Mokhatab differences up to −1.81 °F at 500 psia, 2.4–8.5 °F at 2000 psia.
- Hydration number: 5.895 vs 6.04 ± 0.03 (−2.4 %, "unas 5 veces su incertidumbre" per STATUS) — "COMPARED, NOT VALIDATED (a recorded discrepancy)".
- Hydrate flash: one structure only (other flagged); amount varies steeply near structure change (0.300 → 0.160 mol between 6.85 and 6.90 °C in one case); dense gas labelled liquid carries no water vapour, amount jumps (+32 % between 130 and 140 bar in one case; C1/C2/C3 with 0.05 % water at 15 °C: 3.79e-4 mol at 130 bar, 5.0e-4 at 140 bar); with very little gas (≤ 1e-5) 13–27 s; no ice (T ≥ 0 °C); compositional only.
- Water content at pressure: "NOT VALIDATED" (ideal-vapour model neglects non-ideal enhancement; underestimates at high pressure; not quantified). Hydrate–vapour equilibria not validated. Estimated bias: omits Poynting (+5.5 % in a_w at 70 bar and 10 °C, calculated) and water fugacity coefficient (+10–15 %, estimated): undersaturated hydrate temperature ≈ 0.6–2 K high (estimate, safe side) (`STATUS.md` v1.35.0).
- Inhibitors: methanol above 13 MPa "unvalidated and likely more on the unsafe side" (model depression 30.2 K at 3 MPa, 32.7 at 13, 34.4 at 30, while measured is flat, 31.53–31.64 K at 50 %); over-predicted (unsafe) +0.3 to +1.02 K at 9.7–12.8 MPa; absolute hydrate T up to 1.89 K cold (50 %, 12.76 MPa); max depression error 1.39 K at 50 %.
- MEG: consistency only (no MEG data from another source read).
- Inhibitor applied to gases other than methane without validation: at 7 MPa and 30 % methanol 16.4 K (methane), 14.9 K (C1/C2/C3 sII), 20.9 K (CO₂); black-oil shift ignores gas dependence (1.45 K between the two reference gases; smaller taken).
- Ice points implied by methanol activity 1.0–8.9 K too high at 20–50 %; not used as a limit. Inhibitor loss to gas/oil not modelled; no salts; no ice phase (water treated as liquid down to −60 °C).
- Tried and dropped: activity from ice points / UNIFAC over-predicted depression by up to 6.3 K at 50 % methanol.
- Ng & Robinson (1985), the roadmap criterion for inhibitors, "not accessible".

### Multiphase flow and wells
- Beggs-Brill: not compared with a worked example or measurements; H_L = 0.574 agreement with Brill & Mukherjee Example 4.2 is "FROM MEMORY — not counted as validation". Acceleration term omitted ("<1% for most pipelines"); no-slip friction uses Churchill, not the smooth-pipe Beggs-Brill form.
- Payne et al. factors (0.924 uphill, 0.685 downhill): "UNVERIFIED", not read in a source.
- Beggs-Brill flow-pattern gaps: black-oil tutorial 14 oil line has no solution between 6.87 and 7.55 kg/s, two flows for inlet pressures 0.730–0.799 MPa; flow jumps 6.06 → 7.57 kg/s between 0.72 and 0.73 MPa (ROADMAP v1.26.1). No smoothing planned.
- Hagedorn-Brown and Gray: "NOT VALIDATED against measured well data"; no textbook example accessible; pyResToolbox cross-check "denied by the permission system".
- Duns-Ros, Mukherjee-Brill (3.4): not implemented; Jones IPR not included ("no accessible source").
- IPR (Fetkovich, back-pressure): not validated against a multirate-test worked example.
- Ramey tubing temperature: not validated against measured flowing temperatures (Sagar, Doty & Schmidt 1991 not accessible); default formation properties "typical values, not sourced"; heat of solution of released gas neglected. Default `tubing_u_w_per_m2k = inf` is "a bound, not a best estimate".
- Compositional wells: "NOT VALIDATED against a published gas-well worked example (none accessible)".
- Gas lift: "NOT VALIDATED against published gas-lift performance data"; lift gas assumed at formation temperature at the valve, no mass exchange; annulus not modelled (surface injection pressure not computed).
- `Well` model limitations: P_wh varies ±0.3 % non-monotonically with the number of tubing segments; tabulated curve (60 points) errs up to 2.5 % in P_wh at the high-rate end; compositional liquids use untranslated PR densities (< 1 % for tutorial-09 condensate, 12 % for an n-C16-rich stock-tank oil; untranslated PR made it 671 kg/m³ vs translated 760, overstating oil rates by 13 %).
- Black-oil vs compositional wellhead pressures for oil 2: EOS 13–22 % higher (Standing bubble point 8.1 vs EOS 10.7 MPa).
- Laminar–turbulent jump at Re 2300 (0.0278 vs 0.0475, ×1.71): no solution possible near Re 2300; maintainer's choice pending (Churchill figures "from memory, to be read in the source before use").

### PVT and properties
- Beggs-Robinson dead-oil viscosity below 21.1 °C: extrapolated (API 35: 31 cP at 21 °C, 68 at 15 °C, 5 600 at 0 °C, 3·10⁶ at −6.7 °C); one source says it overstates there; frozen at 21.1 °C tutorial results change ≤ 1.4 % except tutorial 07 (cold oil 80 vs 14 mPa·s; "31.8 % more flow hot" becomes 11.2 %). Original paper not read.
- Black-oil gas (∂h/∂P)_T: −3.9 to −6.1 % at 10 MPa (largest deviation of the data read in source); tolerance 10 % set with results known.
- Black-oil cp(P,T) is not exactly ∂h/∂T: up to 8 % at 200 bar (ROADMAP v1.22.0).
- Heat of solution of dissolved gas neglected (black-oil); black-oil choke temperatures not validated against field or lab data; black-oil compression not validated (no entropy function).
- Compositional liquid viscosity for light liquids with dissolved methane: −25 % to −70 % vs CoolProp for C₁/C₃–C₆ liquids (CHANGELOG v1.10.5 known limitation).
- Surface tension: n-alkane check uses CRC values from memory (regression guard); "mixture IFT data not available" (ROADMAP v1.16.0).
- Phase envelope of mixtures not validated against measured data (the only comparison is pure methane with NIST values cited in the test, read status not stated) (`STATUS.md` v1.28.0, `ROADMAP.md` v1.28.0).
- Kinetic energy not in the energy balance (≈ 1 K at 65 m/s in compositional tutorial 01) (ROADMAP v1.18.0).

### Equipment and tools
- Compressor: no published case read ("GPSA's examples are not in the repository"); no map, surge, staging or intercooling; ratio 10 / 30 MW and 1–30 kW in the test station end not converged; parallel compressors with different specifications do not converge.
- Heater: duties in gas, black-oil and compositional not validated against a published heater (rest on each fluid's enthalpy model); no firing efficiency or approach temperature.
- Separator: "Not compared with published separator tests"; one gas gravity for free and dissolved gas; gas outlets carry ~0.1 % phantom oil (GOR 10⁶); back-flow composition rule is a modelling choice.
- Regulator: "Not compared with a published regulator test or a manufacturer's curve"; wide-open K = 5 is an assumption; no droop, lock-up or valve curve.
- Nodal analysis: no published example; gas-lifted wells not supported.
- Calibration factors: not validated against field data (user's to fit).
- Compositional standard rates: not validated against a published flash; single-stage flash to 60 °F, 1 atm (not a plant's multi-stage path).
- Inner film coefficient: not validated against measurements; two-phase film not modelled (homogeneous); entrance effects not modelled (flagged); transition Re 2300–3000 linear interpolation (no correlation); soil moisture, freezing, transient soil not modelled.
- Looped gas and multiphase networks: not validated against a published case ("No encontré una red mallada de gas ni multifásica publicada con resultados").
- Network at rest with a square-law device (heater, fitting, choke): does not settle (water, K = 135: 0.2 kg/s while pipes carry 1e-5); reported not converged.
- Black-oil cooled below about 0 °C: not converged (extrapolated viscosity + laminar/turbulent jump), diagnosed.
- First outer iteration uses the case default fluid (root cause of tutorial 14 spurious state at 9–10 MPa); not fixed.

---

## Corrections made while writing this manual (v1.39.0)

Collecting the evidence turned up statements that did not agree between
files.  They were resolved as follows (the tables above already use the
corrected labels):

1. **Hanoi network** — its expected heads and flows were generated by
   Angelica itself (Hazen-Williams, C = 130), as the tutorial README says:
   a *regression* baseline, not a comparison with EPANET.  The test comment,
   class docstring and name (`test_hanoi_benchmark_matches_regression_baseline`)
   now say so.  Older CHANGELOG/ROADMAP/STATUS entries that call it an
   "EPANET reference" are corrected by this note.
2. **EPANET pump benchmark** — the values are those of the EPANET 2.2 user
   manual, Appendix C, Listing 2 (the report of the tutorial input of
   Listing 1, same 1000 gpm / 200 ft pump): *read in source* (solver against
   solver).  The Quick Start chapter, also cited, uses a 600 gpm / 150 ft
   pump and is a different case.  An intermediate draft of this manual had
   wrongly labelled it "not verified"; the website's old "within 0.15 psi"
   was also wrong (node 2 differs by 0.212 psi).
3. **Hydrate quadruple points** — the Sloan & Koh table itself was not read
   (secondary source: a search summary of it); the current test values
   (+7.3 % methane Q1, …) supersede the v1.19.0 CHANGELOG figures, which
   predate the dissolved-gas term.
4. **n-alkane surface tension** — the CRC values in the test were written
   from memory: a regression guard, not the validation "within 1–5 % of CRC
   values" that ROADMAP v1.16.0 and the CHANGELOG state.
5. **Methanol ice points** — 0.07 K (test, obtained) and "within 0.1 K"
   (module) agree; the obtained value is 0.07 K.
6. **Z-factor** — the test is a cross-check with the Papay formula, not
   Standing-Katz chart readings; the class docstring now says so.
7. **Phase envelope** — the "HYSYS session data" are not in the repository
   and the test checks broad ranges only: *not validated* for mixtures
   (the only comparison is pure methane with NIST values cited in the test,
   whose read status is not stated).
8. **Isothermal tutorial suite** — its README no longer calls the whole
   suite "validated"; most cases are regression only.
9. **Compressor energy balance** — 0.01 % is the test tolerance, 1.5·10⁻³ %
   the obtained value.
