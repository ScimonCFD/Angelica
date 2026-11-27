# Angelica roadmap — towards an open alternative to PIPEPHASE

Goal: cover PIPEPHASE's main use — steady-state oil and gas gathering
networks with wells. Out of scope for now: transients, optimisation
(NETOPT-style), steam networks, reservoir time-stepping.

This file is the working source of truth. Each item is released on its own
(minor version), validated against published data, and logged at the bottom.

## 1. Where Angelica stands (v1.31.0, checked against the code)

| Area | Present today |
|---|---|
| Fluids | Single-phase liquid, thermal liquid, compressible gas (ideal / PR), black-oil (Standing, Hall-Yarborough, Beggs-Robinson, Lee-Gonzalez-Eakin; enthalpy model), compositional (PR/SRK through `thermo`, free water, phase envelope and quality lines) |
| Pressure drop | Colebrook, Hazen-Williams, Weymouth, Panhandle A/B; Beggs-Brill for black-oil and (opt-in) compositional fluids, with a per-pipe flow-pattern, holdup and erosional-velocity report |
| Heat transfer | Thermal or isothermal (every mode; isothermal option for gas, black-oil and compositional since v1.28.0); fixed U or series layers (flat-wall form); Joule-Thomson, friction heating, latent heat and potential energy along pipes; isenthalpic valves; closed global energy balance |
| Flow assurance | Hydrate curve (vdW-Platteeuw, compositional) and gas-gravity chart (black-oil); hydrate margin per pipe |
| Wells | `WellInflow` and `Well`: IPR linear, Vogel, Fetkovich, gas back-pressure; tubing Beggs-Brill, Hagedorn-Brown or Gray; Ramey temperature profile; black-oil and compositional fluids; gas lift; dead / bistable well states |
| Equipment | Pipe, K-fitting, choke (Sachdeva), check valve, pump by head curve, gas compressor (ratio, rise or power; isentropic or polytropic efficiency), heater/cooler (fixed duty or outlet temperature, with a duty limit), separator (two-phase, three-phase, free-water knockout), pressure regulator (reducing and back-pressure) |
| Boundaries | Pressure, mass flow, temperature, per-inlet black-oil fluid or composition |
| GUI | Canvas editor, all five physics modes, wells and separators, profiles, phase envelope with hydrate curve, SI, US and oil-field units (psia, °F, ft, in, lb/s), Excel export; standard-volume rates (Sm³/d, STB/d, MMscf/d) in the exports |
| Quality | 932 tests, 51 tutorials under golden-output regression (scripts and GUI files), energy balance < 0.01 %; a non-converged solve names its stuck node and alternating flow patterns or laminar/turbulent regimes; correlations used outside their data range are reported as warnings |

**Validation debt (2026-11-25).**  Most models since v1.19 are validated
by internal consistency only: of the last ten log summaries, eight say
"validated against data read in a source: none".  Data read in a source
cover only the black-oil enthalpy (∂h/∂P)_T of gas and water and the
methane solubility in water.  Without a source:
- Hagedorn-Brown, Gray, Ramey, Fetkovich and gas back-pressure;
- gas lift and compositional wells;
- the separator.

Beggs-Brill's coefficients and Ramey's closed form were read in secondary
sources, and the code was verified against them by hand calculation and
the analytic solution.  Neither is compared with a worked example or with
measurements.  The worked examples that
would close this (Brill & Mukherjee, Economides et al., GPSA) are not in
the repository's reference folder (`refs/` does not exist).  Recalling
them from memory does not count.  When the books are available, a
validation pass of these correlations against the worked examples comes
before further stage-7 work, together with item 3.4.

## 2. Gap table against PIPEPHASE

Priority: A = needed to be credible for gathering networks, B = expected by
regular users, C = nice to have. Effort: S ≤ 1 day, M = a few days, L = a week
or more (human-equivalent).

| # | Feature (PIPEPHASE has it) | Angelica today | Pri | Effort | Validation |
|---|---|---|---|---|---|
| 1 | Two-phase slip correlations for compositional fluids | Homogeneous only | A | M | Same case in black-oil vs compositional Beggs-Brill with matched phase properties; Brill & Mukherjee, *Multiphase Flow in Wells* (1999) worked Beggs-Brill example |
| 2 | Holdup, flow pattern, superficial velocities, liquid inventory reported per segment | Computed inside Beggs-Brill, not reported | A | S | Beggs & Brill (1973) flow-pattern boundaries; hand calculation |
| 3 | Hydrate formation curve and margin | None | A | L | Sloan & Koh, *Clathrate Hydrates of Natural Gases* (3rd ed.) data tables: CH4, C2H6, C3H8, CO2 and natural-gas mixtures (Deaton & Frost 1946; Wilcox et al. 1941) |
| 4 | Hagedorn-Brown (wells) | None | A | M | Brill & Mukherjee worked example; Economides et al., *Petroleum Production Systems* examples |
| 5 | Gray (gas-condensate wells) | None | A | M | Gray (1974) / API 14B manual formulae; published gas-well examples (Guo, Lyons & Ghalambor, *Petroleum Production Engineering*) |
| 6 | Duns-Ros, Mukherjee-Brill, Orkiszewski, Aziz, Dukler-Eaton-Flanigan | None | B | M each | Brill & Mukherjee worked examples |
| 7 | Mechanistic models and flow-regime map (Taitel-Dukler-Barnea, Ansari, Xiao; OLGAS is proprietary) | None | B | L | Taitel & Dukler (1976) air-water transition map; Ansari et al. (1994) examples |
| 8 | Well tubing temperature profile (geothermal gradient, Ramey) | Isothermal tubing | A | M | Ramey (1962) analytical solution; Sagar, Doty & Schmidt (1991) field cases |
| 9 | IPR models: Fetkovich, Jones, gas back-pressure C(Pr²−Pwf²)ⁿ | Linear, Vogel | B | S | Fetkovich (1973); Beggs, *Production Optimization Using Nodal Analysis* examples |
| 10 | Wells with compositional fluids (gas and condensate wells) | `Well` is black-oil only | A | M | Consistency with item 5; gas-well deliverability examples |
| 11 | Gas lift (injection point in tubing) | None | B | M | Brown, *The Technology of Artificial Lift Methods*; Economides examples (gradient curves vs GLR) |
| 12 | ESP (stage curves, viscosity and gas derating) | None | C | M | Takacs, *Electrical Submersible Pumps Manual* examples; affinity laws |
| 13 | Separator (phase split node, gas / oil / water outlets) | Done in v1.26.0 (`SeparatorOutlet`) | A | M | Mass and component balance; outlet streams equal an independent `thermo` flash |
| 14 | Compressor (polytropic / isentropic, efficiency, power, discharge T; curves) | Done in v1.30.0 (`Compressor`; no map / curves) | A | M | GPSA *Engineering Data Book* §13 examples; isentropic limit against EOS entropy (done) |
| 15 | Pressure regulator / back-pressure controller | Done in v1.27.0 (`PressureRegulator`) | B | M | Analytical: set-point held when feasible, open or closed otherwise |
| 16 | Heater/cooler with outlet-temperature spec | Done in v1.29.0 (`HeatSource(outlet_temperature_c=…)`) | B | S | Energy balance; water duty against NIST enthalpies |
| 17 | Buried pipe, cylindrical layers, computed inner film coefficient, water/air surroundings | Flat-wall layers, user film coefficients | B | M | Conduction shape factor (Incropera); GPSA / Gnielinski; cylindrical-wall analytic solution |
| 18 | Sensitivity / case-study runner (nodal analysis plots) | Scripts only | A | M | Reproduces tutorial 12 part C and tutorial sweeps from the GUI |
| 19 | Line sizing | None | C | S | Inverse of a solved case |
| 20 | Calibration factors (friction and holdup multipliers) | None | B | S | Factor 1 reproduces current results; monotonic response |
| 21 | Erosional velocity check (API RP 14E) | None | A | S | API RP 14E formula, hand calculation |
| 22 | Oil-field units (bbl/d, MMscf/d, psig, °F, in, ft, cP) and standard-volume rates in results | Done in v1.31.0 (`angelica.units`, per-link standard rates, GUI "Oil field") | A | M | Round-trip conversion tests; factors against NIST SP 811 |
| 23 | Liquid viscosity with dissolved gas (compositional) | Low with dissolved methane | B | M | Pedersen et al. (1984, 1987) corresponding-states method; Lee, Gonzalez & Eakin (1966) methane–n-decane data |
| 24 | Slug and pigging liquid volumes | None | C | M | Brill et al. (1981) Prudhoe Bay slug correlation |
| 25 | Completion pressure drops (perforations, gravel pack) | None | C | S | McLeod (1983); Jones et al. (1976) |
| 26 | Reports: per-link profile tables, standard-volume rates, CSV/Excel | Partial | B | S | Regression tests |
| 27 | Validation manual (benchmarks collected in the docs) | Scattered across tutorials | A | M | — |

PIPEPHASE features deliberately left out: steam, optimisation, time-stepping
with reservoir decline, proprietary models (OLGAS, TACITE), wax deposition.

## 3. Plan

Each item is "done" only when: validated as stated, new tests added, full
suite green, tutorial added where it teaches something, available in the
GUI, released (version, CHANGELOG, tag, website) and CI green.

### Stage 0 — small debts
- 0.1 Round tick values in the phase-envelope plot. *Done when:* axes show
  1-2-5 multiples for any range.
- 0.2 Website guide: document the compositional mode, choke, check valve,
  wells and temperature-dependent properties in the GUI sections.

### Stage 1 — multiphase flow where it is missing (items 1, 2, 21)
- 1.0 Beggs-Brill inclination correction per the published correlation
  (added after review: the code used the Payne multipliers as C).  **Done
  in v1.14.0.**
- 1.1 **Done in v1.15.0.** Beggs-Brill for compositional fluids, using per-phase flash properties
  (three-phase: oil and free water combined as the liquid).
  *Done when:* black-oil and compositional runs with matched phase properties
  agree within 1 %; homogeneous remains the default so existing tutorials do
  not change.
- 1.2 **Done in v1.16.0.** Per-segment report: holdup, flow pattern, superficial velocities,
  mixture velocity, liquid inventory; shown in the GUI link results.
- 1.3 **Done in v1.17.0.** Erosional velocity ratio (API RP 14E) per link, with a GUI warning.
- 1.4 **Done in v1.18.0.** Joule-Thomson cooling along pipes (compressible and compositional).
  Found while doing 1.1: the pipe energy equation has no JT term (only
  valves/fittings are isenthalpic), so a gas line expanding 25 bar stays at
  constant temperature apart from wall heat loss — PIPEPHASE includes it.
  *Done when:* an adiabatic gas pipe follows h(T, P) = const (plus friction
  work consistent with the energy balance), checked against the EOS
  isenthalpic temperature, and the global energy balance still closes.

### Stage 2 — flow assurance (item 3)
- 2.1 **Done in v1.19.0 (validation partly open — see log).** Hydrate curve by van der Waals–Platteeuw (sI and sII, published
  Parrish-Prausnitz / Munck parameters). *Done when:* pure CH4, C2H6, C3H8,
  CO2 and at least two natural-gas mixtures fall within about 1 K of the
  experimental data over 273–295 K; larger errors are reported, not tuned away.
- 2.2 **Done in v1.19.0.** Hydrate curve on the phase-envelope plot; hydrate margin per link in
  the results and on the canvas.
- 2.3 Inhibitors (methanol, MEG). *Done when:* matches published depression
  data (Ng & Robinson 1985) within the accuracy of the chosen method.

### Stage 2b — advanced flow assurance (added 2026-11-23, after stage 4)
- 2b.1 Water content of the gas as an input (ppm, mg/Sm³ or saturation
  fraction); water dew point; hydrate check for undersaturated gas (no free
  water: hydrate–vapour equilibrium with the water in the gas phase).
  *Done when:* a gas at its water saturation reproduces the current curve,
  and a dehydrated gas shifts the hydrate temperature down consistently with
  its water dew point.
- 2b.2 Phase calculation with hydrate (and ice) as additional phases
  (Gibbs minimisation): amount of hydrate formed, water-limited cases,
  composition change. *Done when:* with excess water the incipient boundary
  equals the vdW-P curve; mass balance of water and formers closes.
- 2b.3 Inhibitors (methanol, MEG) on top of 2b.2 — replaces item 2.3.

### Stage 3 — well correlations (items 4, 5, 6)
- 3.1 **Done in v1.21.0.** Hagedorn-Brown (with Griffith bubble-flow correction, no-slip floor).
- 3.2 **Done in v1.21.0.** Gray.
- 3.3 **Done in v1.21.0 (Well part).** Selectable tubing correlation in `Well`; flowline correlation per case.
- 3.4 Duns-Ros and Mukherjee-Brill.
  *Done when:* each reproduces its textbook worked example within 2 %.

### Stage 4 — wells (items 8, 9, 10, 11)
- 4.1 **Done in v1.22.0.** Tubing temperature profile (geothermal gradient + Ramey); the well
  delivers its wellhead temperature to the network.
- 4.2 **Done in v1.23.0 (Jones not included — no accessible source).** IPR: Fetkovich, Jones, gas back-pressure.
- 4.3 **Done in v1.24.0.** `Well` with compositional fluids (gas-condensate wells).
- 4.4 **Done in v1.25.0.** Gas lift: injection depth and rate in `Well`; deliverability vs
  injection rate tutorial.
- Known limitations of the `Well` model (recorded in v1.23.1; fix when one
  gets in the way):
  - Bistable wells (corrected in v1.23.2; v1.23.1 wrongly said every well
    below its minimum stable rate had no steady state).  Between the
    static wellhead pressure and the curve's maximum a well is dead or
    producing.  Rule: producing when consistent; else dead when consistent
    (shut in, its wellhead stays above the static pressure);
    `result.well_states` names the state.  Only when neither is consistent
    (dead it would start, producing it would die — a weak well on a thin
    flowline) is there no steady state: state "unstable",
    converged=False.  Rule chosen as proposed by the maintainer; with
    several bistable wells the retry forces all of them to produce
    together, not every combination.
  - P_wh varies ±0.3 % non-monotonically with the number of tubing
    segments (Beggs-Brill flow-pattern discontinuities; also isothermal).
  - The tabulated curve (60 points, linear interpolation) errs up to 2.5 %
    in P_wh at the high-rate end.
  - Compositional liquids flow with untranslated Peng-Robinson densities
    (pre-existing in the compositional flow model; Peneloux is applied to
    standard volumes and surface tension only).  Below 1 % for the
    tutorial-09 condensate in the tubing; 12 % for an n-C16-rich
    stock-tank oil.  A volume-translation option for the flow model would
    change existing compositional results, so it would be opt-in: not done,
    maintainer to decide.

### Stage 5 — equipment (items 13, 14, 15, 16)
- 5.1 **Done in v1.26.0.** Separator node with gas / liquid (/ water) outlets (`SeparatorOutlet`).
- 5.2 **Done in v1.30.0.** Compressor: polytropic or isentropic with efficiency; power and
  discharge temperature; specified ΔP, pressure ratio or power.  (No
  compressor map, surge or staging.)
- 5.3 **Done in v1.29.0.** Heater/cooler with outlet-temperature spec (and maximum duty).
- 5.4 **Done in v1.27.0.** Pressure regulator (reducing and back-pressure).

- Known limitations (recorded 2026-11-25, from the v1.26.1 review):
  - **Beggs-Brill flow-pattern gaps.**  The pattern boundaries are
    discontinuous, so a pipe can have no solution at a flow imposed on it
    and two solutions at one pressure drop (black-oil tutorial 14's oil line:
    no solution between 6.87 and 7.55 kg/s; two flows for inlet pressures
    between 0.730 and 0.799 MPa).  Which of two solutions a solve finds
    depends on its starting point.  The diagnosis names such a pipe; no
    smoothing of the boundaries is planned (it would no longer be the
    published correlation).
  - **First outer iteration on the default fluid.**  Compositions are
    propagated before the first hydraulic pass, when every flow is zero, so
    the first separator split and regulator check see the case's default
    fluid.  This was the root cause of tutorial 14's spurious state at 9–10
    MPa export (v1.26.1).  The back-flow rule keeps the solver from staying
    there, and regulators skip their first check, but the initialisation
    itself is unchanged.  Not fixed: a fix would change the iteration path
    of every per-inlet case.
  - **Laminar-turbulent jump (recorded 2026-11-25, v1.27.1).**  The
    homogeneous model takes f = 64/Re below Re = 2300 and Colebrook above:
    at ε/D = 3e-4, 0.0278 against 0.0475 (×1.71).  A viscous line whose Re
    sits near 2300 can have no solution; the temperature then moves it
    across and the solve cycles (20 km black-oil line, ambient 6, 5 and
    2 °C).  The diagnosis names the pipe and both regimes.  **Maintainer's
    choice** (not decided here):
    (a) leave it — the published laminar and turbulent laws, a stated gap;
    (b) interpolate f linearly between Re 2300 and 4000 — continuous, no
    published basis for the values in between;
    (c) Churchill (1977), one equation for all Re — 0.0309 at Re 2300,
    within 1.7 % of Colebrook at Re 4000 and 0.6 % at 10⁴ (formula and
    figures from memory, to be read in the source before use).  (b) and (c)
    change results of every pipe between Re 2300 and 4000.
  - **Black-oil below the dead-oil viscosity's data range (recorded
    2026-11-25, v1.27.1).**  Beggs-Robinson's data cover 21.1–146.1 °C and
    16–58 °API (secondary sources; the original not read).  Below 21.1 °C
    the formula is extrapolated (API 35: 31 cP at 21 °C, 68 at 15 °C,
    5 600 at 0 °C, 3·10⁶ at −6.7 °C), and one source says it overstates
    the viscosity there.  The solver now warns.  8 of 15 black-oil tutorials
    are partly below 21.1 °C; frozen at 21.1 °C their results change by
    ≤ 1.4 %, except tutorial 07 (cold oil at 10 °C: 80 against 14 mPa·s;
    "31.8 % more flow hot" becomes 11.2 %).  **Maintainer's choice** (not
    decided here): (a) extrapolate and warn (now); (b) freeze at 21.1 °C
    (no source either; tutorial 07 changes as above); (c) another
    correlation with low-temperature data (to be read in a source first);
    (d) refuse a solve below the range.  The numerical floor −6.7 °C (no
    overflow) stays whatever the choice.

  - **Network at rest with a square-law device (recorded 2026-11-27,
    v1.29.0).**  With every boundary at the same pressure, a heater,
    fitting or choke (ΔP ∝ v²: infinite conductance at v = 0) turns the
    pressure differences left at the stopping tolerance into flow (water,
    K = 135: 0.2 kg/s while the pipes carry 10⁻⁵), and the mass criterion
    is never met.  Reported as not converged with the criteria named, no
    exception (an exception in the gas solver was fixed).  A remedy would
    be a linear (viscous) law below a small velocity — a model change, not
    done.
  - **Black-oil cooled below about 0 °C (recorded 2026-11-27, v1.29.0).**  A
    cooler set to −5 … −20 °C in the two-loop test network: the
    extrapolated viscosity turns the export line laminar and the outer loop
    cycles; not converged, diagnosed.  Follows the two decisions above.

### Stage 6 — heat transfer (item 17)
- 6.1 Cylindrical layers, buried pipe, computed inner film coefficient.

### Stage 7 — engineering tools (items 18, 19, 20, 22, 26)
- 7.1 **Done in v1.31.0.** Oil-field unit system in the GUI; standard-volume rates in results.
- 7.2 Sensitivity runner (API and GUI): vary one or two parameters, plot.
- 7.3 Calibration factors; line sizing.
- 7.4 Per-link profile report export.

### Stage 8 — second wave (items 7, 12, 23, 24, 25, 27)
- 8.1 Pedersen corresponding-states viscosity.
- 8.2 Taitel-Dukler flow-regime map; one mechanistic model (Ansari, vertical).
- 8.3 ESP; completion pressure drops; slug volumes.
- 8.4 Validation manual in the docs.

### Order of work (maintainer review, 2026-11-22)

The objective is PIPEPHASE's main use — oil and gas gathering networks with
wells — so the remaining stages are worked in this order:

1. **Done in v1.20.0. Black-oil brought level** (new item B): Joule-Thomson and friction
   heating in black-oil pipes; hydrate check from the gas gravity (or an
   equivalent composition) with the same outputs as compositional.
2. Stage 3 — well correlations: Hagedorn-Brown and Gray selectable in
   `Well`; then Duns-Ros and Mukherjee-Brill.
3. Stage 4 — wells: tubing temperature profile, Fetkovich and gas IPR,
   `Well` with compositional fluids, gas lift.
4. Stage 5 — equipment: separator, compressor, pressure regulator, heater
   by outlet temperature.
5. Stage 7 — tools: oil-field units, sensitivity studies in the GUI,
   calibration factors.
6. Stage 2b — advanced flow assurance (water content, hydrate phase flash,
   inhibitors).
7. The rest: detailed heat transfer (stage 6), stage 8.

## 4. Standing rules

- No parameter is tuned to make a validation pass. A mismatch is explained
  or recorded as a limitation.
- A change in an existing tutorial result needs a physical reason, written
  in the CHANGELOG.
- New models are opt-in when they would change existing results.
- A vertical two-phase riser as a plain network pipe does not converge
  (multivalued J-curve); wells go through the `Well` component.
- PyPI publication is the last step and is decided by the maintainer.
- (Review 2026-11-22) No tag before the full check passes in both
  environments (ruff, mypy, pytest with thermo 0.5 and 0.6).
- One release per finished roadmap item; bugs found on the way go into
  the same release under "Fixed".
- Every new feature works in black-oil as well as compositional unless it
  is physically impossible — then the reason is recorded.
- Validation data are labelled in the tests: read in the source, read in a
  secondary source, or from memory. Memory does not count as validation.
- Each validation test states the error obtained and where its tolerance
  comes from (data uncertainty, published accuracy of the method); a
  tolerance is never relaxed to make a test pass.
- At the end of each stage the log summarises: validated against data read
  in a source / validated only by internal consistency / not validated.
- (Review 2026-11-25) Every new component is also run outside its
  tutorial's range before tagging: zero flow, reversed flow, dead wells,
  set-points that cannot be met.  Each case must end converged or with
  converged=False and a stated reason, never with an exception.

## 5. Progress log

### 2026-11-22 — Stage 0 and first Stage 1 fixes

- **0.1 Round ticks (v1.13.3).** `angelica.gui.plot_axes` (1-2-5 steps);
  phase-envelope axes use it. Checked visually on a lean gas and a rich
  condensate envelope.
- **0.2 Website guide (with v1.13.3).** GUI sections rewritten from the code:
  case types, material dialogs per mode, all eight devices with their
  availability, link context menu (profiles, phase envelope), unit system
  (inputs follow the selected units — the old guide said SI only), results
  on the canvas.
- **Found while starting 1.1 — three-phase density bug (v1.13.4).** With
  free water the bulk density used the feed molar mass for both HC phases:
  15–18 % low in gas-condensate flow. Fixed (mass over volume, per-phase
  molar mass); trace-water test reproduces the dry EOS density exactly.
  Compositional tutorials 04/05/06 changed (CHANGELOG).
- **1.0 Beggs-Brill inclination (v1.14.0).** Replaced the constant C by the
  published C(λ_L, N_LV, N_Fr) with uphill/downhill coefficient sets; added
  Baker-Swerdloff / gas-water surface tension; Payne multipliers as an
  option (GUI too).
  - Validation: the coefficient tables were cross-checked against the
    whitson wiki Beggs-Brill page. No published worked example could be
    retrieved (Brill & Mukherjee not accessible); the tests use a
    step-by-step hand calculation of their 6-in. vertical data set
    (H_L = 0.574, which matches the value I recall from their Example 4.2 —
    **unverified**).
  - Tutorial 12 changed: +14 % oil (lighter tubing column).
  - Decision: original Beggs-Brill stays the default; Payne is opt-in.
  - Pending: acceleration term (E_k) still omitted; no-slip friction uses
    Churchill (rough Moody) rather than the smooth-pipe Beggs-Brill form.
- **1.1 Beggs-Brill for compositional fluids (v1.15.0).** Flash supplies
  gas/liquid split (condensate + free water as liquid, thermo condensate
  surface tension). Validated by consistency: no-slip density = bulk density,
  C_L = flash phase volumes, network ΔP = Beggs-Brill gradient of the flash
  phases. The roadmap's "black-oil vs compositional with matched phase
  properties" check was not done as such: a black-oil PVT cannot be matched
  to an EOS fluid without tuning; the plumbing tests cover the same risk
  (both paths call the same `beggs_brill_gradient`).
  - Tutorial compositional 07 (hilly gas-condensate line).
  - Also fixed: `free_water_split` put all water in the vapour when the HC
    flash had no vapour (no tutorial affected).
  - Found: no Joule-Thomson term in the pipe energy equation → new item 1.4.
  - CI failed on v1.15.0: thermo's condensate surface tension differs
    between versions (local 0.5 vs CI 0.6). Lesson: anything taken from
    thermo beyond the EOS core (σ, transport models) must be checked in the
    py311 scratch env (`scratchpad/py311`, thermo 0.6.1) before release.
- v1.16.0 lint job failed on mypy (not installed locally). Pre-release
  check now: `scratchpad/precheck.sh` = ruff + mypy (py311) + pytest in
  both environments. Fixed in v1.16.1.
- **1.2 Multiphase report (v1.16.0).** `MultiphaseFlowResult` per pipe,
  CSV/Excel section, GUI window (checked visually by capturing the real
  Toplevel with ImageMagick `import -window`; PIL grabs come out black on
  this desktop). Also replaced thermo σ by Macleod-Sugden + Peneloux
  (pure n-alkanes within 1–5 % of CRC values; mixture IFT data not
  available to me — limitation noted).
- **1.3 Erosional velocity (v1.17.0).** API RP 14E, C = 100, on the no-slip
  density; in the multiphase report, exports, GUI table and a red canvas
  label (the ⚠ glyph is missing from Tk's default font — plain text used).
  - Decision: C fixed at 100 in the report (the function takes any C);
    a per-case C setting can come with the oil-field units work (7.1).
  - Observation: black-oil tutorial 06 runs at ~4.9× V_e (≈ 490 kg/s through
    its lines). Left unchanged — tutorials show physics, not sizing — but
    worth revisiting when the validation manual (8.4) is written.
- **1.4 Joule-Thomson / latent heat / potential energy (v1.18.0).** Lagged
  enthalpy-difference source per pipe (exact at convergence), inner energy
  iteration (relaxation 0.5), enthalpy-conserving junction offsets, real-
  enthalpy global balance. Validated: adiabatic pipes vs exact first law
  (PR methane) and vs thermo PH flash (compositional, incl. condensation and
  elevation) within 0.001 K.
  - Found and fixed on the way: (a) ṁ·cp·T junction mixing made T jump up to
    1.9 K between consecutive pipes in compositional runs; (b) compositions
    mixed by mass instead of molar flow (methane 77.8 vs 79.5 % in tutorial
    02); (c) valves/chokes/pumps did not carry or propagate compositions;
    (d) compositional valves were pass-through (no JT).
  - First attempt (real-gas cp in `specific_heat_for_link`) was abandoned:
    with a P-dependent cp the ṁ·cp·T junctions jump; the enthalpy offsets are
    the right fix.
  - Decision: the new terms are on by default (they complete the first law,
    not an optional model) with a settings switch to turn them off. Black-oil
    and incompressible liquids are unchanged (no enthalpy model): liquid
    friction heating and black-oil JT remain unmodelled — candidates for a
    later item.
  - Pending: kinetic energy (≈1 K at 65 m/s in compositional tutorial 01) is
    not in the balance.
- **2.1/2.2 Hydrates (v1.19.0).** vdW-P with Munck Langmuir constants
  (verified verbatim from Antunes et al. 2018, Table 2, via pdftotext of the
  open PDF) and Holder reference properties (883/808/3.4 sII from Dongre et
  al. 2022; sI 1264/1151/3.0). Hydrate curve on the envelope, per-pipe
  subcooling, canvas label, CSV columns, tutorial compositional 08.
  - Validation: quadruple points (Sloan & Koh, via a ScienceDirect table):
    all within ~0.75 K (or 7 % in pressure).
  - **Not validated (data not accessible to me):** pure methane above 288 K
    and natural-gas mixtures (Nature and ScienceDirect sources refused
    automated access). The roadmap's "two natural-gas mixtures within 1 K"
    criterion is therefore open. Qualitative check only: 3 % propane switches
    to sII and cuts the formation pressure by ~4×.
  - Decision: sI Δμ0 = 1264 J/mol (Munck/Holder-consistent) rather than
    1297 (Handa & Tse) — documented in the CHANGELOG with both errors.
  - Pending: ice branch (T < 0 °C); inhibitors (2.3); gas solubility in
    water (x_w < 1).

### 2026-11-22 — Review corrections (v1.19.1)

Maintainer review of v1.13.3–v1.19.0. Done:
- **CO2 hydrates.** Cause confirmed: dissolved gas ignored (x_w = 1).
  Added Henry's-law water activity with NIST WebBook constants (Sander
  compilation, READ IN SOURCE; no Poynting factor — no sourced partial molar
  volumes). CO2 Q1 now −1.2 %. Before the fix the model was 14–29 % low
  against CO2 points at 277–283 K, but those points are from memory and do
  not count as validation. Methane Q1 +7.3 %.
- `_hydrate_report` returns None when no hydrate is stable above 0 °C
  (it used to report −T as "subcooling").
- Tests: hydrate, surface-tension and Beggs-Brill validation tests now state
  the standing of their data and the errors obtained. The hydrate tolerance
  (11 % in pressure) is the published accuracy of vdW-P + Munck for single
  guests (Antunes et al. 2018, read in source), replacing the 8 % chosen
  after seeing a 6.7 % result. The surface-tension test is labelled a
  regression guard (CRC values from memory; 6 % chosen after the results).
- Energy balance: docstring explains that real-basis and sensible-basis
  percentages are not comparable; first-law tests also bound the absolute
  residual by ṁ·cp·0.01 K.
- **Still open — data needed from the maintainer:** measured Lw-H-V points
  for pure methane 273–295 K, CO2 above 278 K, and two natural-gas mixtures.
  Tried: Nature Sci. Rep. 2020 and ScienceDirect Results in Engineering 2022
  (login/403), MDPI Energies 2025 and Crystals 2021 (403), Deaton & Frost
  1946 (PRCI members only). A full copy of Sloan & Koh found on a course web
  page was not used (unclear licence).

Summary for stage 2 so far:
- Validated against data read in a source: none for hydrates (the
  quadruple points came through a secondary-source summary).
- Validated by consistency: Henry solubility ordering, structure switch
  with propane, curve/inverse consistency.
- Not validated: methane above 288 K, CO2 above 273 K, all mixtures.

### 2026-11-23 — Item B: black-oil brought level (v1.20.0)

- Black-oil enthalpy h(P,T) from (∂h/∂P)_T = Σ w_i(v_i − T ∂v_i/∂T); all of
  the v1.18.0 energy machinery now applies to black-oil (pipes, junctions,
  chokes/valves, real-enthalpy balance). Wells stay isothermal until 4.1;
  their implied heat is a new balance term.
- Hydrate check for black-oil via the Katz gas-gravity chart (no
  composition in black-oil — a composition-based check is impossible).
- Fixed: per-inlet black-oil fluids lost the PVT correlation / gas props.

Summary:
- **Validated against data read in a source:** gas and water (∂h/∂P)_T vs
  NIST WebBook (Setzmann-Wagner, IAPWS-95): 2–6 % and 2–3 %.
- **Validated by consistency:** adiabatic pipes h_in = h_out + g·Δz (0.01 K);
  liquid friction heating between 0 and ΔP/(ρcp) (T·β ≈ 0.31 for the test
  oil); choke isenthalpic; Katz chart reproduces its digitised nodes.
- **Not validated:** the Katz chart itself against measurements (chart
  values read in a secondary source; chart accuracy per literature is a few
  K); the neglected heat of solution of dissolved gas; black-oil choke
  temperatures against field or lab data.

### 2026-11-23 — Stage 3, items 3.1–3.3 (v1.21.0)

- Hagedorn-Brown (modified) and Gray in `closures/well_correlations.py`,
  selectable per `Well`; GUI combobox; tutorial 12 Part D.
- Sources: Gray equations and HB groups/structure from the whitson wiki
  (secondary, read); HB chart-fit coefficients and the Gray roughness-floor
  unit from pyResToolbox source (GPL-3.0, secondary, read — coefficients
  only, own implementation).
- pyResToolbox also confirmed our Beggs-Brill inclination coefficients, and
  notes that the 0.924 uphill "Payne" factor is credited by secondary
  sources but absent from Payne et al. (1979), and that no downhill factor
  is applied. Our optional `payne_correction` applies 0.924/0.685 — to be
  reviewed (pending).
- Attempted cross-check against pyResToolbox installed in the scratch env:
  **denied** by the permission system; not worked around.

Summary:
- Validated against data read in a source: none.
- Validated by consistency: hand calculations of both correlations; dry-gas
  limits; HB and BB share the static column at zero rate.
- Not validated: either correlation against measured well data or a
  textbook example (maintainer could supply Brill & Mukherjee examples).
- Pending in stage 3: Duns-Ros and Mukherjee-Brill (3.4); network-pipe use
  of HB/Gray (3.3, flowline part).

### 2026-11-23 — Review follow-ups (v1.21.3)

No numerical change; four open questions from the review closed or
recorded.
- **Methane hydrates and the Henry term.** The term adds 0.7 % (273 K) to
  7 % (295 K) to the formation pressure, so it does not explain the model
  running high at the warm end. Munck et al. (1988) not accessible —
  whether their fit included solubility is unknown; Antunes et al. (2018)
  include it for all gases with the same constants. Kept for all gases.
- **Hagedorn-Brown coefficients.** Same numbers in pyResToolbox and rNodal;
  rNodal credits Guo, Lyons & Ghalambor (2007). Book not read.
- **Katz chart.** The dissertation never states that its "Experiment"
  column is a chart reading (it says "experimental data retrieved from the
  literature"); recorded as an inference.
- **Payne factors.** Not confirmed in any source; option kept, off by
  default, labelled "unverified" in the GUI and the docstrings.
- Useful find: the NIST ThermoML archive (trc.nist.gov/ThermoML) is open
  and machine-readable — a source of measured data for later validations.

Summary:
- Validated against data read in a source: methane solubility in water
  (Böttger et al. 2016, 13 points inside their expanded uncertainty).
- Compared with secondary-source values: hydrate temperatures at 10–40 MPa
  (methane −1.6 K; CO2 +0.4 / −1.3 K); Katz chart vs Towler-Mokhatab at
  500 psia (within 1.8 °F).
- Not validated: methane hydrate 285–295 K, mixtures, the Katz chart
  against measurements, Hagedorn-Brown / Gray / Beggs-Brill against data.
- Known limitation: CO2 hydrate at 100 MPa is 9 K low (Henry's law without
  Poynting factor).

### 2026-11-23 — Item 4.1: tubing temperature profile (v1.22.0)

- `Well`: isenthalpic inflow, then a coupled pressure–temperature march up
  the tubing with heat loss to the formation (Ramey: wellbore resistance +
  transient rock conduction, linear geothermal profile). The well hands its
  wellhead temperature to the network; `WellInflow` is isenthalpic.
- `isothermal_device_heat_kw` removed; formation heat goes into the wall
  loss and the lift into the potential-energy term.
- Numerics: explicit temperature march diverged at near-zero rates (stiff:
  relaxation distance ≪ segment) → exponential integrator per segment. The
  enthalpy is tracked exactly; the temperature is re-synchronised with it
  every 5 segments (black-oil cp(P,T) is not exactly ∂h/∂T: up to 8 % at
  200 bar). Cost: deliverability curve 0.1 → 0.4 s per well, one 25 ms
  march per well and energy solve.
- Black-oil only for now: there is no compositional `Well` until 4.3; the
  integrator takes generic callables so 4.3 can reuse it.
- Tutorials 11 and 12 changed (CHANGELOG).

Summary:
- Validated against data read in a source: none.
- Verified against an analytical solution: Ramey's profile (own derivation
  of the closed form; f(t) and relaxation distance read in a secondary
  source, the OpenGeoSys benchmark documentation).
- Validated by consistency: first law on black-oil wells and on a network
  with a well; trends with rate, time, insulation.
- Not validated: flowing wellhead temperatures against field data (Sagar,
  Doty & Schmidt 1991 not accessible); default formation properties are
  typical values, not sourced; heat of solution of released gas neglected.

### 2026-11-23 — Item 4.2: Fetkovich and back-pressure IPR (v1.23.0)

- Four IPR models behind one `IPR` object; the reference phase (oil, or gas
  for the back-pressure law) is converted to mass by the fluid. In 4.3 the
  same object will give the gas rate of a compositional well directly.
- Jones left out: neither the original paper nor a worked example was
  accessible, and the roadmap rule asks for a source.
- Done on the way: the pending 4.1 GUI check. Both well dialogs captured;
  labels shortened and defaults stated.
- Note for later validations: ipt.ntnu.no times out from here and PetroWiki
  redirects to onepetro.org (login).

Summary:
- Validated against data read in a source: none.
- Validated by consistency: inverses, derivatives, n = 1 limit, GOR
  conversion, solver vs formula, back-pressure law inside a `Well`.
- Not validated: against a multirate-test worked example (Fetkovich 1973,
  Beggs's examples) — the maintainer could supply one.

### 2026-11-23 — Fix: well dead by back-pressure (v1.23.1)

- Found by an independent review of v1.21.3–v1.23.0: a `Well` whose
  wellhead is held above the maximum of its curve (high separator or
  manifold pressure) gave the right answer (no flow) with converged=False.
- Cause: `ColebrookPipeCorrelation._assign_pipe_state` stored the 1e-6 m/s
  floor of the iteration guesses as the pipe velocity.  The stagnant
  flowline carried ±ρ·A·1e-6 kg/s, flipping sign each iteration; with dense
  oil (~4e-6 kg/s) that exceeded the solver's 'network at rest' threshold
  (1e-6 kg/s).  The other dead-well case (column heavier than P_r) converged
  only because the gassy fluid at 1.5 MPa was light enough to stay below it.
  Fix: store the solved velocity (Hagen-Poiseuille as ΔP → 0); the
  coupling stays the laminar one.  No tutorial result changed.
- Found on the way: in a two-well network with the weak well at the edge of
  its stable rate the energy solve raised OverflowError (1e10 °C at a node
  receiving ~1 kg/s and passing on 1e-9 kg/s on a non-converged flow field).
  Fix: such a node (outflow < 1/1000 of inflow) is mixed as a sink.  A first
  attempt (denominator max(in, out)) changed every tutorial by up to 3.7 K —
  in and out are Σ ṁ·cp with each link's own cp and differ in a converged
  network — and was discarded.  The case itself has no steady state (see
  the Stage 4 limitations).
- Tutorial 12: Part F (separator pressure; Well 2 dies against the
  manifold); README states that the default `tubing_u_w_per_m2k = inf` is
  an upper bound on the heat loss.

Summary:
- Validated against data read in a source: none (numerical fix).
- Verified: stagnant pipe velocity equals Hagen-Poiseuille to round-off;
  dead-well cases converge with Colebrook and Beggs-Brill flowlines; new
  tests fail on v1.23.0 and pass now.

### 2026-11-23 — Fix: bistable wells (v1.23.2)

- Review of v1.23.1: the case it left unconverged has a physical answer.
  Weak well: static wellhead 3.04 MPa, curve maximum 5.49 MPa at
  0.94 kg/s; shut in, the manifold sits at 5.42 MPa, above its static
  pressure, so it cannot start.  Dead is a valid steady state; the link law
  forced the well to produce whenever P_wh < P_max.
- Fix: hysteresis as an active set.  A producing well dies when its
  wellhead reaches the maximum (only in the turbulent phase: the laminar
  warm-start pressures swing 1.5–11 MPa and are not physical — a first
  version that let them switch the well cycled 200 times).  A dead well
  starts again between outer iterations when its wellhead falls below the
  static pressure (at most 3 restarts per solve).  A well left dead inside
  the band is retried forced to produce; producing wins if that solve
  converges.
- `SolveResult.well_states` (`WellOperatingState`): producing /
  dead_column / dead_back_pressure / dead_minimum_rate / unstable, with
  the pressures, the minimum stable rate and a message.
- Item 5 of the review — a real case with no steady state: weak reservoir
  205 bar on its own 1-in., 3 km flowline, separator 5 bar.  Dead: wellhead
  5.25 MPa < static 5.87 MPa (it would start).  Forced to produce: the
  solve does not converge (its minimum stable rate of 0.96 kg/s through the
  thin line exceeds the 7.34 MPa maximum).  Reported "unstable",
  converged=False.  A first version reported it converged-dead because the
  outer loop ran out with the well cycling; fixed (a dead well below its
  static pressure makes the solve unconverged).
- The 1e-6 kg/s absolute 'network at rest' threshold did not decide
  convergence in any of the cases run for this fix.
- No tutorial number changed; tutorial 12 text corrected.

Summary:
- Validated against data read in a source: none (logic of the steady
  states).
- Verified: the five two-well cases of the review and the thin-flowline
  cases give the state the rule predicts; the bistable case where producing
  is consistent gives v1.23.1's rates unchanged.

### 2026-11-23 — Item 4.3: compositional wells (v1.24.0)

- One fluid adapter for the well models (`angelica.closures.well_fluid`):
  black-oil wraps the existing correlations (all black-oil results
  bit-identical, checked by regenerating every tutorial reference);
  compositional uses a PVT table per well — a flash costs 10–15 ms and a
  curve needs tens of thousands of evaluations (direct flashing: ~3 min per
  curve; table: 5–10 s).
- Table error, measured against a march flashing at every step: P_wh
  0.05 %, T_wh 0.01 K.  A first, geometric-only pressure grid had 55 bar
  between nodes at 250 bar and enthalpy errors up to 2.9 kJ/kg in dense gas
  (P_wh 0.45 %, well heat 2.8 % off); now at most 10 bar between nodes.
- Standard volumes from a flash to 15.56 °C, 1 atm with Peneloux-translated
  phase volumes: untranslated PR made an n-C16-rich stock-tank oil 671
  kg/m³ (translated 760; n-hexadecane is 773), overstating oil rates by
  13 %.
- Energy: a table fluid's well heat is taken from the EOS enthalpies at
  both ends, so the network balance closes exactly; the independent check
  is that it equals the table march's own heat (0.07 %).
- The bistable-well logic moved to a mixin shared by both solvers
  (`angelica.solvers.wells`); `WellInflow` also works with compositional
  streams.
- Black-oil vs compositional, same volatile oil (black-oil GOR, API and gas
  gravity taken from the EOS stream's standard flash), linear IPR, 2 500 m
  tubing: the mass per stock-tank m³ agrees to 0.1 %.  Wellhead pressures:
  oil 1 (GOR 203, API 55) within 2–5 %; oil 2 (GOR 86, API 52) EOS 13–22 %
  higher.  Causes: the EOS column is lighter (untranslated liquid density,
  see the limitation above) and Standing's bubble point differs from the
  EOS one (8.1 vs 10.7 MPa for oil 2).  Nothing adjusted.
- Tutorial compositional 09 (gas and condensate wells).  GUI: the Well and
  Well IPR buttons now appear in compositional mode; the dialog is the same.

Summary:
- Validated against data read in a source: none.
- Validated by consistency: table vs direct flash; dry gas vs an
  independent single-phase integration (0.045 %); standard-volume mass
  against the black-oil formula (0.1 %); first law; IPR-only link
  reproduces the back-pressure law.
- Not validated: against a published gas-well example (none accessible).

### 2026-11-24 — Item 4.4: gas lift (v1.25.0)

- The lift gas enters the network at the wellhead node as an internal
  source.  Rejected alternatives: (a) the well link carrying production +
  lift gas — the reservoir node would report the lift gas as production
  and the inlet enthalpy would count it at reservoir conditions; (b) a
  third terminal for a lift-gas supply node — not expressible with
  two-node links.
- Tubing: two sections, lift gas mixed without mass transfer, arriving at
  the formation temperature of the valve depth (no annulus model).  For
  compositional fluids this avoids a PVT table per rate (the mixture
  composition changes along the deliverability curve); the lift gas has its
  own table.
- Found on the way: black-oil junction mixing averaged API, gas gravity,
  GOR and WOR by mass — not conserving standard volumes, and the mixed
  fluid's enthalpy 2–4 K·cp off; a gas-only stream could not be added.
  Fixed (standard-volume mixing); five black-oil tutorials changed (largest:
  tutorial 05 mixing temperature +3 K, now between the two inlet
  temperatures).
- Checked: the performance curve rises, peaks and falls; deeper valve =
  more oil and higher valve pressure (tutorial 13, numbers checked before
  writing the text).

Summary:
- Validated against data read in a source: none.
- Validated by consistency: continuity, first law of the march, network
  mass/energy balances, standard-volume mixing (exact) and its enthalpy
  consistency, compositional molar mixing.
- Not validated: gas-lift performance against published data; lift-gas
  temperature at the valve (assumed equal to the formation's).

### 2026-11-24 — Item 5.1: separator (v1.26.0)

- A separator is a node with `SeparatorOutlet` links leaving it
  ("liquid", "oil" + "water", or "water" alone).  The node's feed is
  split at its P and T; an outlet's flow is its phases' mass fraction of
  the current feed (the level-control valve: flow set by the split, ΔP
  reported), the gas leaves through the node's other links or its
  pressure boundary by mass balance.  Rejected alternative: a separator
  as a boundary on both outlets — overdetermined (P_sep and the split
  cannot both follow from two downstream pressures); a real separator's
  liquid valve absorbs that difference.
- The pressure correction includes the outlet's dependence on the feed
  (δṁ_out = f·Σ δṁ_feed): with the outlet flow lagged, an LP separator in
  series (f ≈ 1) diverged.  Compositions are propagated again after the
  hydraulic pass when there are separators, so the first energy solve does
  not see default compositions downstream of the outlets.
- Black-oil split = the PVT's own in-situ split (oil with Rs, free gas,
  water); compositional = PR flash + immiscible water.  The outlet streams
  propagate as their own compositions.
- Found on the way: the energy equation of a boundary node with
  throughflow (a pressure boundary that also passes flow on) divided the
  inflowing enthalpy by the outflow only: −256 °C in the new separator
  train tutorial, 3.6 K too hot in a black-oil test.  Fixed (mixing
  temperature of the arriving streams).  No existing tutorial had such a
  node: references unchanged.
- Checked with numbers before writing: stock-tank liquid has a maximum in
  HP separator pressure (≈18 bar one stage, ≈30 bar two stages) and a
  second stage gains 2–5 % (compositional tutorial 10); in black-oil the
  stock-tank oil does not depend on the separator pressure (tutorial 14).

Summary:
- Validated against data read in a source: none.
- Validated by consistency: black-oil gas fraction = in-situ free-gas mass
  fraction (round-off); compositional split = independent phase-split
  routine (1e-9) with species conserved; enthalpy additivity of the split
  (≤ 0.0015 K·cp black-oil, 6e-6 K·cp compositional); network mass and
  energy balances; GUI file = script.
- Not validated: against a published separator test or plant data.
  Black-oil: one gas gravity for free and dissolved gas; gas outlets carry
  ~0.1 % phantom oil (GOR 10⁶ convention).

### 2026-11-25 — Fix: separator outside the tutorials' range (v1.26.1)

An external review swept black-oil tutorial 14's export pressure beyond
the tutorial's 0.5–5 MPa: an exception at 7–8 MPa, converged=False at
9–10 MPa.  Compositional tutorial 10 was swept too (export 30–150 bar,
flare 5–40 bar): an exception with the export above the manifold.  Each
cause was found before changing anything:

- 7–8 MPa: the per-link black-oil fluid was looked up by object identity,
  and the energy assembly passes links wrapped with a temperature, so every
  link had the default fluid's cp.  This is pre-existing and affects every
  black-oil network with per-inlet fluids.  The enthalpy offsets kept
  results right at convergence to within 0.09 K and 190 Pa (the tutorial
  reference changes), but the water outlet line, cooled with the mixture's
  cp, fell below the dead-oil viscosity correlation's range.
- 9–10 MPa: the first outer iteration splits the default fluid (flows are
  still zero), the gas line reverses, and the export boundary then
  injected the default fluid — water and low-GOR oil — which kept the gas
  fraction near zero.  The state was spurious and self-sustaining.  Rule
  adopted: back-flow from an outlet boundary with no declared fluid is the
  separator's own gas.  With a gas phase present, the back-flow cannot be
  steady (feed + B = f·(feed + B) needs f = 1), so the iteration leaves
  that state.  Without a gas phase it is undetermined: converged=False with
  a stated reason (shown in black-oil and compositional).
- Export above the manifold (compositional): the hydraulics converged
  (gas flows back through the HP separator to the inlet), but the stagnant
  liquid branch kept ~1e-11 kg/s round-off flows, which the energy
  equation divided by: −10⁸ °C.  Flows below 1e-9 of the largest are now
  zero for the energy solve.
- Found on the way: tutorial 14's 0.5 MPa row (v1.26.0) was printed as a
  solution but had not converged.  The oil LCV imposes 7.12 kg/s, and
  Beggs-Brill's intermittent/distributed boundary makes the oil line's
  flow jump from 6.06 to 7.57 kg/s between 0.72 and 0.73 MPa inlet
  pressure: no solution exists in the model.  Measured by solving the
  line alone at 0.70–0.94 MPa.  Not a solver defect; the row now says
  so.
- After the fixes: tutorial 14 is converged from 0.6 to 24 MPa (wells dying
  at ~5 and ~12 MPa); tutorial 10 is converged at export 30–150 bar and at
  flare 5–25 bar.  Flare at 32–40 bar (above the HP separator) ends with
  converged=False and the no-steady-state reason.

Summary:
- Validated against data read in a source: none (numerical fixes).
- Validated by consistency: the separator's water fraction at 7 and 9 MPa
  equals Well 1's split (only Well 1 is alive); balances close.
  Tutorial changes stay within 0.09 K and 190 Pa, and the energy residuals
  fell.
- Not validated: the separator against plant data or a published test
  (unchanged); the back-flow composition rule is a modelling choice, not
  a measured behaviour.

### 2026-11-25 — Item 5.4: pressure regulator (v1.27.0)

- `PressureRegulator`, reducing (holds its end node) or back-pressure
  (holds its start node), in black-oil and compositional.  Three states
  (controlling / open / closed) chosen between outer iterations.  While
  controlling, the held node is a fixed-pressure row and the regulator's
  flow enters the pressure correction implicitly, as for separator outlets.
  Rejected alternative: a lagged flow, which is the separator mistake of
  v1.26.0.
- Problems met while sweeping beyond the tutorials, each fixed at its
  cause:
  - Next to a separator node the regulator's flow is (1 − F) of a change in
    the feed: without that factor in the Jacobian the hydraulics diverged.
  - State checks after an unconverged pass saw transients (−2 000 kg/s
    through a regulator carrying 0.7) and closed it: states are checked only
    after a converged pass, or after two unconverged ones in a row (the
    state itself may be what cannot converge).
  - Imposing a demand the supply cannot meet drove the upstream pressure to
    zero within a pass (compositional flash raised): the imposed flow is
    capped at the wide-open valve's flow.  A capped regulator in an
    unconverged pass opens, unless it reached controlling from open (then
    control is feasible by continuity).
  - Energy solves on unconverged flow fields gave absurd temperatures that
    set the next pass's densities: in networks with regulators the
    temperatures are held after an unconverged pass.  Networks without
    regulators are unchanged.
  - Zero flow with the held node exactly at the set point: closed (flows
    below 1e-5 kg/s are zero).
- Found on the way:
  - The black-oil enthalpy inversion returned the bracket midpoint (2.5 K
    off, alternating): a check-valve case cycled.  Fixed; tutorial results
    move by at most 1e-4 K.
  - Black-oil properties raised or overflowed outside −6.7…500 °C during
    unsettled iterations: they are now evaluated within it, and a solution
    outside it is reported as not converged.  (Corrected in v1.27.1: −6.7 °C
    is a numerical floor, not the correlation's limit, and out of range is
    now a warning.)
- Review items of v1.26.1: the tutorial 14 explanation was corrected (a
  no-solution band in flow, with two solutions in pressure); the
  non-convergence diagnosis was added; the initialisation root cause was
  recorded under stage 5.
- Checked with numbers before writing:
  - Separator held at exactly 3 MPa below 3 MPa export; open above, equal
    to tutorial 14 (459.5 m³/d at 5 MPa in both).
  - City gate: 0.62 K/bar let-down cooling, about 260 kW of preheat for
    0 °C at the outlet, control from about 36 bar of supply.
  - The flow still rises slightly with supply pressure while controlling,
    through the outlet temperature.  A first draft of the tutorial said it
    did not; corrected before release.

Summary:
- Validated against data read in a source: none.
- Validated by consistency:
  - controlling flow = the downstream network solved alone at the set
    point (1.2e-7);
  - open = a check valve of the same K (1e-7);
  - isenthalpic outlet against the fluid's h(P, T) (cp·0.01 K);
  - fine sweeps converge and are monotonic in both modes and both fluid
    models.
- Not validated: against a published regulator or station case.  The
  wide-open K (5) is an assumption; there is no droop, lock-up or valve
  curve.

### 2026-11-25 — Review fixes of v1.27.0 (v1.27.1)

- Regulators in series: a BPR ahead of a reducer left the line between
  them with no pressure of its own (singular, 10¹⁰ bar, exception in the
  compositional solver).  The giving side of such a line now opens.
  Sweeping set points found three more failures (two reducers switching
  back to control together; a reducer closing on a transient zero flow; a
  network at rest at −337 °C), all fixed.  Off-suite sweep: 5 layouts ×
  (black-oil 81 + compositional 64) = 725 solves, all converged and
  consistent.
- Black-oil enthalpy below −6.7 °C: (∂h/∂P)_T jumped (h fell 35 kJ/kg in
  1 K): +16.6 K let-down heating and a node at −22.7 °C.  Frozen with the
  properties; h continuous and increasing.
- Beggs-Robinson range: 21.1–146.1 °C, 16–58 °API (secondary sources);
  −6.7 °C was unsourced.  Out of range is now `SolveResult.warnings`, not
  converged=False.  Effect of the extrapolation on the tutorials measured
  (see limitations).
- A non-converging 20 km black-oil line at 5 °C ambient: the
  laminar-turbulent jump, not the regulator.  Diagnosed now; options in
  the limitations.
- Found on the way: tutorial black-oil 10 prints 1.289 MPa where its
  reference has 1.290 — the same at v1.27.0 run directly, so the
  environment the reference was generated in; within the golden test's
  last-digit tolerance.

Summary:
- Validated against data read in a source: the data range of the dead-oil
  viscosity (secondary sources, two agreeing); nothing else.
- Validated by consistency: regulator pair states against their pressures
  (725 solves); BPR + reducer flow = the downstream pipe solved alone
  (1.3e-7); h continuous and increasing below the floor, let-down there
  equal to the one at the floor (0.02 K).
- Not validated: the regulator against a published case; the
  low-temperature viscosity (no data below 21.1 °C).

### 2026-11-27 — Goal check, part 1 (v1.28.0)

Autonomous run towards the goal "PIPEPHASE's basic capabilities for oil,
gas and water gathering networks" (summaries in Spanish in `STATUS.md`).
Part 1 checks what exists, with numbers, before new items.

- Looped networks: a network with two independent loops and two inlets in
  every fluid mode (water, methane PR, black-oil with gas and water,
  compositional with free water and condensate), thermal and isothermal.
  Junction imbalance ≤ 1.07e-5; each pipe solved alone with the network's
  flow gives the network's pressure within 3.02e-5 of its ΔP (the loop
  law itself holds by construction with one pressure per node — a loop-sum
  test was circular and was removed in review); energy balance ≤ 2.77e-4 %;
  hydrate margin per pipe equal to the stand-alone functions.
- Found missing: no isothermal option in the gas, black-oil and
  compositional solvers.  Added (`isothermal_temperature_c`, GUI
  Numerics), off by default.
- Found on the way: the laminar-turbulent jump also stops an isothermal
  black-oil solve (5–6 °C in the two-loop network) — the discontinuity
  itself, not the thermal coupling.  Recorded for the maintainer's choice.
- GUI Numerics dialog: an entry hidden under another (fixed).
- Review before the tag: GUI still required a fixed-temperature boundary in
  black-oil and compositional; NaN / −300 °C accepted (exceptions in gas,
  a "converged" black-oil at −300 °C) — temperatures are now checked;
  devices whose temperature change the isothermal mode drops are named in
  the warnings.  Not fixed (pre-existing): the gas solver returns no
  diagnosis when it does not converge.

Summary:
- Validated against data read in a source: nothing new (the published
  looped network is still Hanoi, water, EPANET results).
- Validated by consistency: everything above.
- Not validated: looped gas and multiphase networks against a published
  case; mixture phase envelopes against measured data.

### 2026-11-27 — Item 5.3: heater/cooler by outlet temperature (v1.29.0)

- `HeatSource(outlet_temperature_c, max_duty_w)` in every thermal mode; duty
  = ṁ·Δh with the fluid's enthalpy; "duty limited" beyond the maximum;
  `SolveResult.heaters`; GUI fields and canvas label; CSV/workbook.
- Black-oil tutorial 16 (heater station), compositional tutorial 11 Part D.
- Found on the way, fixed: fixed-duty heaters in gases did not conserve
  energy (0.39 % at 500 kW, constant cp); a transparent heater stalled a
  water network (coupling cap); a gas network at rest with a fitting raised
  (starved junction, −310 °C); gas and thermal-liquid solvers gave no
  diagnosis.
- Recorded: network at rest with a square-law device; black-oil cooled
  below ~0 °C.
- Review before the tag: my first coupling fix (0.01 m/s floor) broke a
  heater after a dead well (phantom 0.025 kg/s) — floor set to 1e-4 m/s,
  test added; the gas fixed-duty error was 7.1 % (methane) and 5.9 %
  (condensate) in the test network, not only 0.39 %; node-after-heater
  accuracy is ≈ 1e-3 × T(°C), not the 0.004 K of one network.

Summary:
- Validated against data read in a source: water duty against NIST
  WebBook enthalpies (+0.048 %, −0.044 %; tolerance 0.51 % from the NIST cp
  spread).
- Validated by consistency: outlet = set point; duty = ṁ·Δh (7.4e-6);
  energy balance ≤ 1.9e-4 %; duty limit; reversed flow; set points from
  −60 to 300 °C.
- Not validated: duties in gas, black-oil and compositional fluids against
  a published heater (they rest on each fluid's enthalpy model).

### 2026-11-27 — Item 5.2: compressor (v1.30.0)

- `Compressor` (ratio, rise or power; isentropic or polytropic) in gas,
  black-oil and compositional; ∫v·dP marched in ln P on each fluid's h and
  v; `SolveResult.compressors`; GUI; compressible tutorial 06 Part B.
- Development found: a ratio taken from the inlet pressure inside the
  pressure correction ran away (5 000 kg/s) and a power at zero flow asked
  for the ratio limit — both now set once per outer iteration by a
  bracketed false-position search on the rise; a node at 0 Pa made the EOS
  divide by zero (floor 1 kPa in absolute-pressure solvers).

- Review before the tag found the first version judged convergence on the
  imposed rise, not the node pressures (a station with no steady state
  reported converged; a small nozzle missed its ratio), the 1 kPa floor
  overriding a 0.5 kPa boundary and a flow-boundary exception — fixed.
- Limitations recorded: ratio 10 / 30 MW and 1–30 kW in the test station;
  parallel compressors with different specifications.

Summary:
- Validated against data read in a source: none.
- Validated analytically / against the EOS: ideal-gas isentrope (8.8e-5);
  isentrope against EOS entropy, PR (1.3e-4) and compositional (4.7e-5).
- Validated by consistency: specification met (1e-4), power = ṁ·Δh,
  energy balance ≤ 1.5e-3 %.
- Not validated: against a published compressor case (GPSA examples not
  available); black-oil compression (no entropy function).

### 2026-11-27 — Item 7.1: oil-field units and standard rates (v1.31.0)

- `angelica.units` (exact definitions); per-link standard-volume rates in
  black-oil and compositional results and exports (Sm³/d, STB/d, MMscf/d,
  bbl/d); GUI "Oil field" unit system with °F and psia.

Summary:
- Validated against data read in a source: conversion factors against
  NIST SP 811 Appendix B.8 (7 significant figures).
- Validated by consistency: black-oil standard rates equal the separator's
  (oil ≤ 1e-6, water 1.6e-6; gas 0.11 %, the gas-only representation).
- Not validated: compositional standard rates against a published flash
  (single-stage flash to 60 °F, 1 atm; not a plant's multi-stage path).

