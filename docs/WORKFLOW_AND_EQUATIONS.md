# Angelica v1: High-Level Workflow and Governing Equations

This note summarises what the current `Angelica` implementation is doing at a
high level and which equations it is solving.

The scope of the current project is:

- steady
- incompressible
- pressure-driven pipe networks with pipes, fittings, pumps, and inline heat sources
- isothermal mode
- non-isothermal mode with an outer temperature loop

This note focuses primarily on the hydraulic core implemented in
`src/angelica/solvers/steady_isothermal_incompressible.py`, because that same
pressure-correction machinery is reused inside the non-isothermal solver.

## 1. High-Level Workflow

The current isothermal workflow in
`src/angelica/solvers/steady_isothermal_incompressible.py` is:

1. Build the network state from the case definition.
2. Initialise a pressure field for all internal nodes.
3. Run a laminar initialisation stage.
4. Use that state as the initial condition for the turbulent stage.
5. Return nodal pressures and component flow rates.

The non-isothermal workflow in
`src/angelica/solvers/steady_non_isothermal_incompressible.py` wraps that
hydraulic solve in an outer temperature loop:

1. Build the network state and initialise thermal boundary conditions.
2. Solve the hydraulic problem with the current temperature-dependent
   properties.
3. Solve the steady energy equation on that flow field.
4. Update nodal and component temperatures.
5. Repeat until the temperature change becomes sufficiently small.
6. Run one final synchronised hydraulic + energy pass for the reported result.

In more detail:

### Step 1. Build the network

The case defines:

- nodes
- pressure boundaries
- components (`Pipe`, `Fitting`, `Pump`)
- fluid properties

This is converted into a mutable internal state with:

- nodal pressures
- component velocities
- Reynolds numbers
- mass flow rates

### Step 2. Initialise internal pressures

For nodes without prescribed pressure, the code generates a simple
interpolated pressure field between the available pressure boundaries.

This is only a numerical starting point. It is not yet the hydraulic
solution.

### Step 3. Laminar initialisation

The solver computes a provisional flow field using:

- laminar pipe physics for `Pipe`
- the fitting local-loss relation for `Fitting`
- the pump curve relation for `Pump`

Then it assembles and solves a pressure-correction system and updates the
nodal pressures.

This is repeated:

- once for networks with pipes only
- several times for networks with fittings

The idea is:

- pure laminar pipes are effectively linear
- fittings are already nonlinear
- therefore a network with fittings benefits from a few laminar-style
  initialisation passes before the turbulent stage

### Step 4. Turbulent solve

The solver then switches the pipes to the turbulent Darcy-Weisbach plus
Colebrook relation.

At each iteration it:

1. computes a velocity field from the current pressures
2. computes mass flows
3. linearises each component into a coupling coefficient
4. assembles the pressure-correction system
5. solves for pressure correction
6. updates the nodal pressures
7. checks convergence

### Step 5. Convergence

The current turbulent stage stops when the pressure correction and nodal
mass imbalance both become sufficiently small.

Colebrook itself has its own internal residual tolerance.

## 2. Equations Implemented

The current code follows the structure from the original draft reasonably
closely.

Notation:

- `P_i`, `P_j`: pressures at the start and end nodes of a component
- `DeltaP = P_i - P_j`
- `z_i`, `z_j`: elevations
- `Delta z = z_j - z_i`
- `rho`: density
- `mu`: dynamic viscosity
- `g`: gravitational acceleration
- `D`: pipe diameter
- `L`: pipe length
- `A`: cross-sectional area
- `f`: Darcy friction factor
- `K`: fitting loss coefficient
- `V`: average velocity in the component
- `m_dot = rho A V`: mass flow rate

## 2.1 Pipe, laminar

Implemented in `src/angelica/closures/friction.py`.

The current laminar pipe velocity is:

```text
V = (D^2 / (32 mu L)) * (DeltaP - rho g Delta z)
```

This matches the intended form in the draft:

```text
V ~ P_in - P_out - gamma (z_out - z_in)
```

where `gamma = rho g`.

Important:

- the elevation term enters the estimated velocity
- it does not enter the laminar pressure-correction coefficient directly

The laminar coupling used in the pressure-correction system is:

```text
C_lam = -(rho / (32 mu)) * (A D^2 / L)
```

## 2.2 Pipe, turbulent

Also implemented in `src/angelica/closures/friction.py`.

The turbulent pipe update is based on Darcy-Weisbach:

```text
DeltaP - rho g Delta z = (rho f L / (2 D)) * V * abs(V)
```

The present code rearranges this using the previous velocity estimate `V*`:

```text
V_new = 2 D (DeltaP - rho g Delta z) / (rho f L abs(V*))
```

This is consistent with the usual linearisation of:

```text
V * abs(V)
```

instead of:

```text
V^2
```

which is important because `V * abs(V)` preserves flow direction.

The turbulent coupling used in the pressure-correction equation is:

```text
C_turb = -2 A D / (f V L)
```

In practice the implementation uses the current signed velocity. This is one
place where sign sensitivity matters.

## 2.3 Pipe friction factor

For Re < 2300 (laminar regime) the friction factor is taken directly from the
Hagen-Poiseuille result:

```text
f = 64 / Re
```

For Re ≥ 2300 the friction factor is found from the Colebrook equation:

```text
1 / sqrt(f) = -2 log10( eps/(3.7 D) + 2.51/(Re sqrt(f)) )
```

with:

```text
Re = rho abs(V) D / mu
```

The current implementation solves Colebrook iteratively by Newton-style
finite-difference updates.

## 2.4 Fittings

Implemented in `src/angelica/closures/minor_losses.py`.

The fitting velocity is currently computed from:

```text
V = sign(DeltaP) * sqrt( 2 abs(DeltaP) / (K rho) )
```

So fittings are not treated as a strictly linear laminar law. They retain
their local-loss relation even during the laminar initialisation stage.

This is important:

- the current "laminar initialisation" is laminar in the pipes
- but not fully linear across the whole network if fittings are present

The fitting coupling used in the pressure-correction system is:

```text
C_fit = -2 A / (K abs(V))
```

## 2.5 Pumps

Implemented in `src/angelica/closures/pump.py`.

Pumps are treated as pressure-changing devices driven by a Q-Head curve.

The solver first converts the current pressure difference across the pump into
the head the pump must deliver:

```text
H_req = -DeltaP / (rho g)
```

Then it evaluates the pump curve locally to obtain the volumetric flow rate:

```text
Q = G(H_req)
```

Two pump modes are supported:

- a one-point EPANET-style design point, expanded internally to a full curve
- a piecewise-linear multi-point Q-Head table

The velocity follows from:

```text
V = Q / A
```

For pressure correction, the pump contributes the local slope of the curve:

```text
C_pump = dQ / d(DeltaP) = (dQ / dH) * 1/(rho g)
```

This is the pump analogue of the hydraulic conductance used for pipes and
fittings. It lets the pump enter the same global continuity system as any
other component.

## 2.6 Mass conservation and pressure correction

The global system is assembled in `src/angelica/numerics/assembly.py`.

For each component:

```text
m_dot = rho A V
```

and each component contributes:

- a mass imbalance term to the right-hand side
- a coupling coefficient to the matrix

In matrix form the solver is building:

```text
M p' = b
```

where:

- `p'` is the pressure correction
- `b` is built from the current mass-flow imbalance
- `M` comes from the component coupling coefficients

Boundary-pressure nodes are pinned by imposing:

```text
p' = 0
```

at those nodes.

## 3. Non-Isothermal Thermal Model

### 3.1 Governing assumptions

The non-isothermal solver retains the **incompressible** assumption from the
hydraulic core.  Concretely:

- Fluid density is **not** a function of pressure.  ρ, μ, cₚ, and k may all
  vary with temperature, but they are evaluated at the local pipe temperature
  only — not at the local pressure.
- This is the correct model for **liquids** (water, thermal oils, crude oil)
  under pressures typical of building and pipeline networks (up to a few tens
  of bar), where pressure-induced property changes are negligible compared to
  temperature-induced ones.
- **Compressible gases** — where ρ = ρ(P, T) — are outside the scope of this
  solver.  Using it for a gas network would give incorrect density and mass
  flow results.
- Only **steady-state** conditions are modelled.  Thermal capacitance, pipe
  wall heat storage, and hydraulic inertia are not included.

### 3.2 Energy equation (pipe FV discretisation)

For each pipe segment the solver assembles the steady convection–diffusion
energy equation:

```text
ṁ cₚ (T_E - T_W) = k A / Δx (T_E - T_P) - k A / Δx (T_P - T_W)
                  + U π D Δx (T_amb - T_P)
```

where the last term is the Moukalled-style source linearisation for heat loss
to the surroundings:

```text
S = Sc + Sp · T_P
Sc = U π D Δx T_amb
Sp = -U π D Δx   (≤ 0, strengthens the diagonal)
```

For inline heat sources the wall heat-loss term is replaced by the fixed power
source `Q / n_segments`.

### 3.3 Convection scheme

Three face-interpolation schemes are available (selectable via `convection_scheme`):

| Scheme | Face coefficients |
|--------|-------------------|
| Upwind (default) | `a_W = max(ṁ cₚ, 0)`, `a_E = max(-ṁ cₚ, 0)` |
| Hybrid | blend of central and upwind depending on Peclet number |
| Power Law | Patankar power-law blend |

Upwind is unconditionally stable and first-order accurate in space.  For
networks where diffusion is negligible compared to convection (Pe ≫ 1, which
is typical of pipe flow), all three schemes give essentially the same answer.

### 3.4 Analytical reference solution (NTU method)

For a single pipe with uniform U, ṁ, and cₚ, the exact steady-state solution is:

```text
T_out = T_amb + (T_in - T_amb) · exp(-NTU)
NTU = U π D L / (ṁ cₚ)
```

The automated test suite verifies that the FV solver converges to this
solution as the number of segments increases (`test_mesh_convergence_approaches_ntu_analytical`).

An external textbook benchmark (Cengel & Ghajar, *Heat and Mass Transfer*,
5th ed., McGraw-Hill, Example 8-3 — oil through an icy lake) provides an
independent reference with published result T_out = 19.74 °C, reproduced by
Angelica within 0.05 K.
