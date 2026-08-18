# 01 — Natural Gas Pipeline

Branched methane network demonstrating the compressible solver.

**Fluid:** methane (CH₄), ideal gas, T = 15 °C  
**Network:** one inlet, interior junction, two outlets  
**Inlet:** 800 kPa &nbsp;·&nbsp; **Outlets:** 500 kPa each

## Expected output

```
Case: Natural Gas Pipeline (methane, 800→500 kPa)
Converged: True

Node pressures:
  Node 1 (inlet):     800.00 kPa   rho = 5.263 kg/m³
  Node 2 (junction):  ~640  kPa    rho = ~4.2  kg/m³
  Node 3 (outlet A):  500.00 kPa   rho = 3.289 kg/m³
  Node 4 (outlet B):  500.00 kPa   rho = 3.289 kg/m³

Compressibility effect:
  Ratio (rho_in / rho_out): 1.600  (= P_in / P_out for ideal gas)
```

## Running

```bash
python tutorials/steady_compressible/01_natural_gas_pipeline/run.py
```

Or open `natural_gas_pipeline.gui.json` in Angelica via **File → Open**.
