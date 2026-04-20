---
name: pk-validation-cmax-vpc
description: >
  PK validation procedures for the TOURMALINE synthetic datasets: correct
  computation of population median Cmax, VPC construction, NCA targets, and
  documentation of the known Ixazomib Cmax overprediction issue inherent to
  the Gupta 2017 3-compartment model parameters.
applies_when:
  - Validating Ixazomib, Lenalidomide, or Dexamethasone PK against NCA targets
  - Cmax validation is giving an inflated or deflated value
  - Building or interpreting a VPC plot for the PK data
  - Deciding whether to attempt fixing the Cmax overprediction
  - Extending the PK model with new covariates (BSA, CrCL, CYP3A4)
keywords:
  - Cmax, AUC, NCA, VPC, PK, Ixazomib, Lenalidomide, Dexamethasone
  - Gupta 2017, 3-compartment, PCTESTCD, PCSTRESN, BLQ, sdtm_pc, adam_adpc
  - CYP3A4, BSA, CrCL, population PK, IIV
load_cost: medium
---

## Level 2 — Instructions

### Known Limitation: Ixazomib Cmax Overprediction
Population median Cmax from the Gupta 2017 3-cmt model: ~53 ng/mL (MM2) / ~57 ng/mL (MM1).
Published target: **41 ng/mL** (+30%/+38% overprediction).

**Root cause**: Ka=0.5 h⁻¹ (fast absorption) combined with V1 (central volume)
produces a sharp, high peak that the sparse-sampling-based published estimates
smooth out. This is a pre-existing model parameter issue.

**Do not attempt to fix** without refactoring the ODE solver and re-fitting Ka, V1
against dense PK profiles. No simple parameter adjustment will resolve it while
preserving AUC/t½ accuracy.

### Correct Cmax Calculation
Cmax = **population median of per-subject maximum concentrations** (not global max).

```python
# CORRECT
sub_pk = pc[(pc["PCTESTCD"]=="IXAZOMIB") & (pc["BLQ"]=="N")].copy()
sub_pk["CONC"] = pd.to_numeric(sub_pk["PCSTRESN"], errors="coerce")
cmax_median = sub_pk.dropna(subset=["CONC"]).groupby("USUBJID")["CONC"].max().median()

# WRONG — inflated ~4× because global max picks the single highest observation
cmax_wrong = sub_pk["CONC"].max()
```

### VPC Construction
1. Compute 5th/50th/95th percentile of simulated concentrations at each nominal time
2. Overlay individual observed concentrations as scatter points
3. BLQ observations excluded from percentile computation (use `BLQ=="N"` filter)
4. Cycle 1 (dense) and Cycle 3 (sparse) plotted separately
5. 80% prediction interval (10th–90th) is the pharmacometric standard

### Covariate Effects (implemented in generate_pk_v2.py)
| Covariate | Effect | Magnitude |
|-----------|--------|-----------|
| BSA on V4 | V4 = TV_V4 × (BSA/1.73)^0.70 | Power function |
| CrCL on Lenalidomide CL | CL = TV_CL × (CrCL/80)^0.60 | Power function |
| CYP3A4 strong inhibitor | CL × 0.55 | 45% reduction |
| CYP3A4 strong inducer | CL × 2.0 | 100% increase |

---

## Level 3 — Resources

### NCA Targets (Ixazomib, Gupta 2017)
| Parameter | Published | Tolerance | Current sim |
|-----------|-----------|-----------|-------------|
| Cmax median | 41 ng/mL | ±15% | 44–45 ng/mL (+8–10%) ✓ — read from sdtm_pp.csv Cycle 1 NCA |
| AUC0-24 | ~300 ng·h/mL | ±20% | calibrated ≈ ✓ |
| t½ apparent | 228h (9.5 days) | ±20% | ~228h ✓ |
| CL (3-cmt) | 1.86 L/h | ±20% | calibrated ✓ |
| IIV CV% | ~36% | ±10% | 36% ✓ |

### PK Model Parameters (Gupta 2017 Table 3)
```
Ixazomib:    CL=1.86 L/h, V1=?, V2=?, V3=?, Q2=?, Q3=?, Ka=0.5/h, F=58%
             BSA on V4: exponent=0.70, reference=1.73 m²
Lenalidomide: CL/F=8.94 L/h, CrCL covariate, 1-compartment
Dexamethasone: CL/F=16 L/h, 1-compartment
Residual error: proportional + additive
LLOQ: Ixazomib=0.5 ng/mL, Lenalidomide=2.0 ng/mL, Dex=0.2 ng/mL
```

### File Locations
| File | Contents |
|------|----------|
| `scripts/generate_pk_v2.py` | 3-cmt ODE, IIV, NCA, sparse/dense sampling schedule |
| `MM2/sdtm_pc.csv` | PK concentration records |
| `MM2/sdtm_pp.csv` | NCA-derived parameters (Cmax, AUC, t½) |
| `MM2/adam_adpc.csv` | PK analysis dataset with covariates |
| `scripts/validate_data.py` → `plot_pk_vpc()` | VPC figure generation |
| `scripts/validate_data.py` → `write_validation_table()` | Cmax validation row |
