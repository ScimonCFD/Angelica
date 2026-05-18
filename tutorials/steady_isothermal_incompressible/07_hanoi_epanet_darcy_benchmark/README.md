# Hanoi EPANET/Darcy Benchmark

This folder contains a `Angelica` benchmark case for the classical Hanoi water
distribution network, configured to match a published EPANET/Darcy-style
reference as closely as possible.

The goal of this case is not only to provide a large `Open -> Run` GUI example,
but also to serve as a reproducible validation case for the current
steady, incompressible, isothermal, single-component solver.

## What the Hanoi Case Is

The Hanoi network is a classical benchmark from the water-distribution-network
literature. It is a looped municipal-style network with:

- one fixed-head reservoir
- many demand nodes distributed across the network
- several interconnected loops
- a mix of long trunk pipes and smaller branch pipes

It is widely used because it is large enough to be nontrivial, but still small
enough to be studied, reproduced, and benchmarked across different solvers.
For a network solver, Hanoi is useful as a test of:

- mass conservation over a multi-loop network
- pressure/head propagation over long paths
- flow redistribution through loops
- robustness on a case much larger than toy examples

In short, this is not just a tutorial network. It is a recognized benchmark for
checking whether a hydraulic network solver behaves correctly on a realistic
looped system.

## Reference

The benchmark data used here come from:

- Vuta, L.I. and Piraianu, V. (2008), *Infoworks WS and EPANET v2 - Modeling
  the water distribution networks*

In that reference:

- the Hanoi network is treated as a flat network
- the source is a single fixed-head reservoir at `100 m`
- the hydraulic comparison is given in Darcy/EPANET style
- published tables include expected pipe flows and nodal heads

## Benchmark Definition Used in `Angelica`

This `Angelica` version of the benchmark uses:

- a single `100 m` fixed-head source
- Darcy-Weisbach / Colebrook-White hydraulic closure
- equivalent absolute roughness `0.2 mm`
- trunk pipe lengths from the published table
- trunk pipe diameters from the published table
- nodal demands from the published table

## Adaptation to the `Angelica` Node Model

The classical Hanoi benchmark places demands directly on network nodes.

The current `Angelica` core uses a stricter node model:

- `junction` nodes are conservative by construction
- net inflow/outflow is imposed only through `source` and `sink` nodes

Because of that, each published nodal demand is represented here as:

- one conservative trunk `junction`
- one short side branch
- one dedicated `sink`

To keep this adaptation hydraulically close to the reference network:

- the side branches are intentionally made almost lossless
- terminal demand nodes `13` and `22` are collapsed into their sink nodes
  to keep the scene cleaner without changing the effective terminal behavior

## What Is Compared

The benchmark comparison is performed on the **trunk network only**:

- trunk-node heads, in `m`
- trunk-link volumetric flows, in `m^3/h`

This is important:

- the lateral `sink` nodes are part of the `Angelica` adaptation
- they are not benchmark nodes from the published Hanoi reference
- their displayed pressures are therefore not the validation target

## Expected vs Calculated Tables

The script
[benchmark_tables.py](/home/simon/Documents/Angelica/tesisIca/ClasesTesis/Angelica/Angelica/tutorials/steady_isothermal_incompressible/07_hanoi_epanet_darcy_benchmark/benchmark_tables.py)
prints:

- expected trunk-node heads from the published reference
- calculated trunk-node heads from `Angelica`
- the difference `delta = Angelica - expected`
- expected trunk-link flows from the published reference
- calculated trunk-link flows from `Angelica`
- the difference `delta = Angelica - expected`

One subtle point:

- some published pipe flows are reported with a direction convention that may
  differ from the internal link orientation used in `Angelica`
- therefore a few links can appear with opposite sign while still matching in
  magnitude
- in those cases, the hydraulic interpretation is still consistent; only the
  orientation convention differs

## Practical Interpretation

For this benchmark, a good result means:

- the case converges
- trunk-node heads are close to the published values
- trunk-link flows are close to the published values

In the current validated setup, `Angelica` reproduces the reference very closely:

- trunk heads typically differ by only a few tenths of a meter
- trunk flows typically differ by only tiny fractions of `m^3/h`

That makes this case a strong validation benchmark for the current solver on a
nontrivial looped network.

## Useful Commands

Run the tutorial case:

```bash
python3 tutorials/steady_isothermal_incompressible/07_hanoi_epanet_darcy_benchmark/run.py
```

Print the benchmark comparison tables:

```bash
python3 tutorials/steady_isothermal_incompressible/07_hanoi_epanet_darcy_benchmark/benchmark_tables.py
```

Open the case in the GUI:

```bash
PYTHONPATH=src python3 -m angelica.gui.app
```

then use:

- `File -> Open`
- open
  [hanoi_epanet_darcy_benchmark.gui.json](/home/simon/Documents/Angelica/tesisIca/ClasesTesis/Angelica/Angelica/tutorials/steady_isothermal_incompressible/07_hanoi_epanet_darcy_benchmark/hanoi_epanet_darcy_benchmark.gui.json)
- press `Run`
