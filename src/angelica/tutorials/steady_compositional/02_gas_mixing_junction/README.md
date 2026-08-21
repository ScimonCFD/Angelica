# Tutorial 02 — Compositional Gas Mixing Junction

Two gas sources with different three-component compositions feed a common
junction.  The solver propagates each inlet's mole-fraction vector
downstream and computes a molar-flow-weighted blend at the junction node.

## Network

```
[Node 1] ──PipeA (5 km, D=0.12 m)──┐
 P=120 bar, T=70 °C                  ├── [Node 3] ──PipeC (8 km, D=0.18 m)── [Node 4]
 Rich gas: CH₄ 90%, C₂H₆ 8%        │    junction                              P=25 bar
                                     │
[Node 2] ──PipeB (5 km, D=0.10 m)──┘
 P=110 bar, T=60 °C
 Lean gas: CH₄ 60%, C₂H₆ 30%, C₃H₈ 10%
```

## Fluid

| Property | Value |
|---|---|
| Components | methane, ethane, propane |
| Source A (rich)  | CH₄ 90%, C₂H₆ 8%, C₃H₈ 2% |
| Source B (lean)  | CH₄ 60%, C₂H₆ 30%, C₃H₈ 10% |

## Composition propagation

In each outer iteration the solver:
1. Seeds node compositions from the inlet boundary conditions.
2. Follows the flow direction through each pipe.
3. Mixes incoming streams at junctions by molar-flow-weighted average.
4. Writes the result to `PipeState.zs` before the hydraulic solve.

## Run

```bash
python tutorials/steady_compositional/02_gas_mixing_junction/run.py
```

## Expected output (abridged)

```
Converged:         True
Outer iterations:  4

Node   P (bar)    T (°C)
   3    69.23      50.61   ← junction pressure from SIMPLE
   4    25.00      35.82

Blended composition in pipe_C (mass-weighted):
  methane  : 78.43 mol%
  ethane   : 16.49 mol%
  propane  :  5.09 mol%
```

The blend is weighted by the mass flows of pipes A and B, which differ
because of the pressure and diameter difference between the two inlet pipes.
