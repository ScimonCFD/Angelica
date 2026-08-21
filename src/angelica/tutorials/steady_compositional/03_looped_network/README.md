# Tutorial 03 — Compositional Looped Gas Network

Two gas sources at different pressures and compositions feed a five-pipe
looped gathering network.  The loop pipe (PipeD) creates a hydraulic
cycle — the SIMPLE solver must find the pressure distribution that
satisfies both nodal mass balance and Kirchhoff's pressure law around
the loop.  The composition in the loop pipe, and its flow direction,
emerge from the hydraulic solution.

## Network

```
[Node 1] ──PipeA (5 km, D=0.12 m)──→ [Node 3] ──PipeC (8 km, D=0.15 m)──→ [Node 5]
 P=120 bar, T=70 °C                      │                                    P=20 bar
 Rich gas: CH₄ 90%                      PipeD (3 km, D=0.10 m)  ← loop pipe     ↑
                                          │                                       │
[Node 2] ──PipeB (5 km, D=0.12 m)──→ [Node 4] ──PipeE (8 km, D=0.15 m)────────┘
 P=110 bar, T=60 °C
 Lean gas: CH₄ 60%, C₂H₆ 30%, C₃H₈ 10%
```

One independent loop: Node 3 → Node 4 → Node 5 → Node 3 (via PipeD, PipeE, PipeC).

## Fluid

| Property | Value |
|---|---|
| Components | methane, ethane, propane |
| Source A (rich)  | CH₄ 90%, C₂H₆ 8%, C₃H₈ 2% |
| Source B (lean)  | CH₄ 60%, C₂H₆ 30%, C₃H₈ 10% |

## Key things to observe

**PipeD flow direction** — PipeD is defined start→end as 3→4.  Since
Node 1 feeds Node 3 at 120 bar and Node 2 feeds Node 4 at 110 bar,
Node 3 ends up at a slightly higher pressure than Node 4, so the loop
pipe carries gas in the 3→4 direction.  This is not prescribed — it
emerges from the pressure solution.

**Composition at Node 4** — receives lean gas from PipeB *and* a small
stream of rich gas from PipeD, producing an intermediate blend.

**Delivered composition at Node 5** — a weighted mix of Node 3 (pure
rich gas via PipeC) and Node 4 (blend via PipeE).

## Run

```bash
python tutorials/steady_compositional/03_looped_network/run.py
```

## Expected output (abridged)

```
Converged:         True
Outer iterations:  4

Node  Label           P (bar)    T (°C)
   3  Junction 3       67.42      55.33
   4  Junction 4       67.07      46.15   ← slightly lower than Node 3

Pipe       ṁ (kg/s)   Direction
pipe_D      0.63      ** 3 → 4 **   ← loop flow direction confirmed

Compositions (mol%):
            Node 3    Node 4   Node 5 (outlet)
methane      90.00     61.61       74.89
ethane        8.00     28.82       19.08
propane       2.00      9.57        6.03

Global mass balance error: 0.0002 %
```

The small loop flow (0.63 kg/s vs ~11 kg/s from each source) is
determined entirely by the pressure imbalance around the loop.
