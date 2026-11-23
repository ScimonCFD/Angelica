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
- 4.3 `Well` with compositional fluids (gas-condensate wells).
- 4.4 Gas lift: injection depth and rate in `Well`; deliverability vs
  injection rate tutorial.
- Known limitations of the `Well` model (recorded in v1.23.1; fix when one
  gets in the way):
  - No steady state below the minimum stable rate.  A well whose
    deliverability curve peaks at q > 0 can flow only at or above that
    rate; if the network needs less, the solver reports converged=False
    (the real well heads).  Options, not decided: report it as a named
    state ("unstable"), or a time-averaged closure.
  - P_wh varies ±0.3 % non-monotonically with the number of tubing
    segments (Beggs-Brill flow-pattern discontinuities; also isothermal).
  - The tabulated curve (60 points, linear interpolation) errs up to 2.5 %
    in P_wh at the high-rate end.

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
