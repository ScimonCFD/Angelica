# Network Setup

A `NetworkCase` collects all simulation inputs: components, boundary conditions, fluid model, and solver settings hints.

## Components

Components are immutable dataclasses (frozen). Each connects two node IDs.

```python
from angelica import Pipe, Fitting, Pump, HeatSource

# Pipe: the main flow element
pipe = Pipe(
    start_node=1, end_node=2,
    diameter_m=0.2,
    length_m=500.0,
    roughness_m=4.6e-5,          # Darcy friction (steel)
    height_change_m=10.0,        # elevation rise (positive = uphill)
    heat_transfer_coefficient_w_per_m2k=5.0,  # wall U-value (non-isothermal)
    ambient_temperature_c=15.0,
)

# Fitting: minor loss (K coefficient)
fitting = Fitting(start_node=2, end_node=3, diameter_m=0.2, loss_coefficient=0.5)

# Pump: defined by a Q-H curve (pairs of (flow_m3/s, head_m))
pump = Pump(
    start_node=3, end_node=4,
    diameter_m=0.15,
    curve_points_q_head=((0.0, 80.0), (0.05, 60.0), (0.10, 20.0)),
)

# HeatSource: adds fixed power (W) to the fluid
heater = HeatSource(
    start_node=4, end_node=5,
    diameter_m=0.1, length_m=2.0,
    power_w=50_000.0,
)
```

## Boundary Conditions

```python
from angelica import PressureBoundary, FlowBoundary, ThermalBoundary, InletFluidBC

# Pressure BCs
inlet  = PressureBoundary(node_id=1, pressure_pa=6e5)
outlet = PressureBoundary(node_id=5, pressure_pa=1e5)

# Flow BCs (alternative to pressure at the inlet)
flow_in = FlowBoundary(node_id=1, mass_flow_kg_per_s=5.0)

# Thermal BC (non-isothermal solvers only)
thermal_in = ThermalBoundary(node_id=1, bc_type="fixed_temperature", temperature_c=80.0)
```

## Assembling the Case

```python
from angelica import NetworkCase

case = NetworkCase(
    name="my_network",
    components=(pipe, fitting, pump),
    pressure_inlets=(inlet,),
    pressure_outlets=(outlet,),
    thermal_inlets=(thermal_in,),   # only for thermal/compressible modes
    fluid=fluid,
)
```

All sequences are converted to tuples internally. Node IDs are arbitrary positive integers; they must form a fully connected graph.

## Accessing Results

```python
result = solver.solve(case)

print(result.converged)                           # bool
print(result.node_pressures_pa[2])                # Pa at node 2
print(result.node_temperatures_c[2])              # °C at node 2 (thermal only)
for flow in result.component_flows:
    print(flow.label, flow.mass_flow_kg_per_s)

# Mass and energy balance checks
print(result.global_balance.mass_error_pct)
print(result.global_energy_balance.energy_error_pct)  # thermal only
```
