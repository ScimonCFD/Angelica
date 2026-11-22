# Angelica roadmap — towards an open alternative to PIPEPHASE

Goal: cover PIPEPHASE's main use — steady-state oil and gas gathering
networks with wells. Out of scope for now: transients, optimisation
(NETOPT-style), steam networks, reservoir time-stepping.

This file is the working source of truth. Each item is released on its own
(minor version), validated against published data, and logged at the bottom.

## 1. Where Angelica stands (v1.13.2, checked against the code)

| Area | Present today |
|---|---|
| Fluids | Single-phase liquid, thermal liquid, compressible gas (ideal / PR), black-oil (Standing, Hall-Yarborough, Beggs-Robinson, Lee-Gonzalez-Eakin), compositional (PR/SRK through `thermo`, free water, phase envelope and quality lines) |
| Pressure drop | Colebrook, Hazen-Williams, Weymouth, Panhandle A/B; Beggs-Brill (black-oil only). Compositional two-phase flow is homogeneous no-slip |
| Heat transfer | Fixed U or series layers (film + wall + insulation + film, flat-wall form); Joule-Thomson; shaft work; closed global energy balance |
| Wells | `WellInflow` (linear, composite Vogel); `Well` = IPR + isothermal Beggs-Brill tubing, black-oil only |
| Equipment | Pipe, K-fitting, choke (Sachdeva), check valve, pump by head curve (also used as a "compressor"), fixed-duty heater/cooler |
| Boundaries | Pressure, mass flow, temperature, per-inlet black-oil fluid or composition |
| GUI | Canvas editor, all five physics modes, profiles, phase envelope, SI and US units (ft, psi, lb/s), Excel export |
| Quality | 511 tests, 43 tutorials under golden-output regression, energy balance < 0.01 % |

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
| 13 | Separator (phase split node, gas / oil / water outlets) | None | A | M | Mass and component balance; outlet streams equal an independent `thermo` flash |
| 14 | Compressor (polytropic / isentropic, efficiency, power, discharge T; curves) | Pump head curve only | A | M | GPSA *Engineering Data Book* §13 examples; isentropic limit against EOS entropy |
| 15 | Pressure regulator / back-pressure controller | None | B | M | Analytical: set-point held when feasible, open or closed otherwise |
| 16 | Heater/cooler with outlet-temperature spec | Fixed duty only | B | S | Energy balance |
| 17 | Buried pipe, cylindrical layers, computed inner film coefficient, water/air surroundings | Flat-wall layers, user film coefficients | B | M | Conduction shape factor (Incropera); GPSA / Gnielinski; cylindrical-wall analytic solution |
| 18 | Sensitivity / case-study runner (nodal analysis plots) | Scripts only | A | M | Reproduces tutorial 12 part C and tutorial sweeps from the GUI |
| 19 | Line sizing | None | C | S | Inverse of a solved case |
| 20 | Calibration factors (friction and holdup multipliers) | None | B | S | Factor 1 reproduces current results; monotonic response |
| 21 | Erosional velocity check (API RP 14E) | None | A | S | API RP 14E formula, hand calculation |
| 22 | Oil-field units (bbl/d, MMscf/d, psig, °F, in, ft, cP) and standard-volume rates in results | SI and partial US | A | M | Round-trip conversion tests |
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
- 1.1 Beggs-Brill for compositional fluids, using per-phase flash properties
  (three-phase: oil and free water combined as the liquid).
  *Done when:* black-oil and compositional runs with matched phase properties
  agree within 1 %; homogeneous remains the default so existing tutorials do
  not change.
- 1.2 Per-segment report: holdup, flow pattern, superficial velocities,
  mixture velocity, liquid inventory; shown in the GUI link results.
- 1.3 Erosional velocity ratio (API RP 14E) per link, with a GUI warning.

### Stage 2 — flow assurance (item 3)
- 2.1 Hydrate curve by van der Waals–Platteeuw (sI and sII, published
  Parrish-Prausnitz / Munck parameters). *Done when:* pure CH4, C2H6, C3H8,
  CO2 and at least two natural-gas mixtures fall within about 1 K of the
  experimental data over 273–295 K; larger errors are reported, not tuned away.
- 2.2 Hydrate curve on the phase-envelope plot; hydrate margin per link in
  the results and on the canvas.
- 2.3 Inhibitors (methanol, MEG). *Done when:* matches published depression
  data (Ng & Robinson 1985) within the accuracy of the chosen method.

### Stage 3 — well correlations (items 4, 5, 6)
- 3.1 Hagedorn-Brown (with Griffith bubble-flow correction, no-slip floor).
- 3.2 Gray.
- 3.3 Selectable tubing correlation in `Well`; flowline correlation per case.
- 3.4 Duns-Ros and Mukherjee-Brill.
  *Done when:* each reproduces its textbook worked example within 2 %.

### Stage 4 — wells (items 8, 9, 10, 11)
- 4.1 Tubing temperature profile (geothermal gradient + Ramey); the well
  delivers its wellhead temperature to the network.
- 4.2 IPR: Fetkovich, Jones, gas back-pressure.
- 4.3 `Well` with compositional fluids (gas-condensate wells).
- 4.4 Gas lift: injection depth and rate in `Well`; deliverability vs
  injection rate tutorial.

### Stage 5 — equipment (items 13, 14, 15, 16)
- 5.1 Separator node with gas / liquid (/ water) outlets.
- 5.2 Compressor: polytropic or isentropic with efficiency; power and
  discharge temperature; specified ΔP, pressure ratio or power.
- 5.3 Heater/cooler with outlet-temperature spec.
- 5.4 Pressure regulator.

### Stage 6 — heat transfer (item 17)
- 6.1 Cylindrical layers, buried pipe, computed inner film coefficient.

### Stage 7 — engineering tools (items 18, 19, 20, 22, 26)
- 7.1 Oil-field unit system in the GUI; standard-volume rates in results.
- 7.2 Sensitivity runner (API and GUI): vary one or two parameters, plot.
- 7.3 Calibration factors; line sizing.
- 7.4 Per-link profile report export.

### Stage 8 — second wave (items 7, 12, 23, 24, 25, 27)
- 8.1 Pedersen corresponding-states viscosity.
- 8.2 Taitel-Dukler flow-regime map; one mechanistic model (Ansari, vertical).
- 8.3 ESP; completion pressure drops; slug volumes.
- 8.4 Validation manual in the docs.

## 4. Standing rules

- No parameter is tuned to make a validation pass. A mismatch is explained
  or recorded as a limitation.
- A change in an existing tutorial result needs a physical reason, written
  in the CHANGELOG.
- New models are opt-in when they would change existing results.
- A vertical two-phase riser as a plain network pipe does not converge
  (multivalued J-curve); wells go through the `Well` component.
- PyPI publication is the last step and is decided by the maintainer.

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
