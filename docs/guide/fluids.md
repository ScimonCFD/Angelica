# Fluid Models

Each solver requires a fluid model that provides density, viscosity, and (for thermal/compressible modes) specific heat and thermal conductivity.

## SingleComponentFluid

For isothermal incompressible simulation with constant properties.

```python
from angelica import SingleComponentFluid

fluid = SingleComponentFluid(
    density_kg_per_m3=998.2,
    viscosity_pa_s=1.002e-3,
)
```

## ThermalFluid

For non-isothermal simulation. Properties can be constants or callables of temperature (°C).

```python
from angelica import ThermalFluid

fluid = ThermalFluid.from_constants(
    density_kg_per_m3=870.0,
    viscosity_pa_s=5e-3,
    specific_heat_j_per_kg_k=2100.0,
    thermal_conductivity_w_per_m_k=0.13,
)

# Or from T-dependent callables:
fluid = ThermalFluid.from_functions(
    density_fn=lambda T: 900.0 - 0.6 * T,
    viscosity_fn=lambda T: 0.01 * (1 - 0.02 * T),
    specific_heat_fn=lambda T: 2000.0 + 2.0 * T,
    thermal_conductivity_fn=lambda T: 0.14,
)
```

Built-in helpers: `build_water_thermal_fluid()`, `build_thermal_dead_oil(api_gravity)`.

## CompressibleFluid

For compressible gas simulation. Density comes from an equation of state.

```python
from angelica import CompressibleFluid
from angelica import IdealGasEOS, PengRobinsonEOS

eos = IdealGasEOS(molecular_weight_kg_per_mol=0.016043)  # methane
fluid = CompressibleFluid.from_constants(
    eos=eos,
    viscosity_pa_s=1.1e-5,
    specific_heat_j_per_kg_k=2200.0,
    thermal_conductivity_w_per_m_k=0.033,
)
```

## BlackOilFluid

For three-phase (gas/oil/water) black-oil simulation.

```python
from angelica import BlackOilFluid

fluid = BlackOilFluid(
    api_gravity=35.0,
    gas_specific_gravity=0.75,
    water_cut=0.0,
    gor_scf_per_stb=600.0,
)
```

## CompositionalFluid

For multi-component EOS-based simulation using the `thermo` library.

```python
from angelica import CompositionalFluid

fluid = CompositionalFluid.from_names(
    component_names=["methane", "ethane", "propane"],
    mole_fractions=[0.85, 0.10, 0.05],
    eos="PR",
)
```
