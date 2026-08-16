# 02 — Flow BC Cross-Validation

Same methane network as tutorial 01, but with **mass-flow boundary conditions** at the two outlets instead of prescribed pressures.

**Fluid:** methane (CH₄), ideal gas, T = 20 °C  
**Network:** one pressure inlet, interior junction, two flow outlets  
**Inlet:** 800 kPa &nbsp;·&nbsp; **Outlets:** 1.2555 kg/s and 1.5407 kg/s

The outlet mass flows match the solution from the equivalent all-pressure-BC case (800 kPa → 500 kPa). Running this case with those flow values and checking that the solver recovers ≈ 500 kPa at both outlet nodes validates the mixed pressure/flow BC solver path.

## Expected output

```
Case: Natural Gas Pipeline — Flow BC Cross-Validation
Converged:            True
Density iterations:   2

Node pressures:
  Node 1:   800.00 kPa   rho = 5.266 kg/m³
  Node 2:   667.33 kPa   rho = 4.393 kg/m³
  Node 3:   500.06 kPa   rho = 3.292 kg/m³
  Node 4:   500.06 kPa   rho = 3.292 kg/m³

Cross-validation: outlet pressures with flow BCs
  Node 3: 500.06 kPa  (should match pressure-BC case ≈ 500 kPa)
  Node 4: 500.06 kPa  (should match pressure-BC case ≈ 500 kPa)
```

## Running

```bash
python tutorials/steady_compressible/02_flow_bc_cross_validation/run.py
```

Or open `flow_bc_cross_validation.gui.json` in Angelica via **File → Open**.
