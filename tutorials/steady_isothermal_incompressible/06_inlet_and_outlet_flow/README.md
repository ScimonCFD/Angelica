# Tutorial 06 — Inlet and Outlet Flow Boundaries

Same network as tutorial 02 (fittings, no elevation), with **both** the inlet
at node 1 and the outlet at node 6 replaced by flow boundaries. Only node 5
keeps a pressure boundary, which anchors the absolute pressure level.

## Network layout

Identical topology to tutorial 02.

| Node | Type   | Condition                    |
|------|--------|------------------------------|
| 1    | Source | 1.7053 kg/s (flow inlet)     |
| 5    | Source | 201 300 Pa                   |
| 6    | Sink   | 2.5530 kg/s (flow outlet)    |

Mass balance: Q_in_node1 + Q_in_node5 = Q_out_node6
(the flow at node 5 is implicit; pressure boundary sets the head reference).

## What this tutorial demonstrates

- Mixed flow + pressure boundary configuration with two flow boundaries
- Only a single pressure boundary is required to fix the absolute pressure level
- Most general boundary-condition setup short of an all-flow-boundary system

## Run

```bash
python3 tutorials/steady_isothermal_incompressible/06_inlet_and_outlet_flow/run.py
```
