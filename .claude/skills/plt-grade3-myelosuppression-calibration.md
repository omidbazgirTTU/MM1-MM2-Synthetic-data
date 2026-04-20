---
name: plt-grade3-myelosuppression-calibration
description: >
  Calibration of Grade 3 thrombocytopenia (PLT < 50 × 10⁹/L) rates in synthetic
  trial data using a within-cycle sinusoidal dip model. Covers the correct
  per-patient-worst definition, log-linear dip amplitude scaling, and arm/study-
  specific baseline differences between NDMM and RRMM populations.
applies_when:
  - Grade 3 PLT rate is off target for any arm × study combination
  - Changing dip_amp doesn't produce the expected Grade 3 rate change
  - PLT rates drifted after any unrelated change to the LB generator
  - Adding a new study arm that needs independent PLT calibration
keywords:
  - PLT, thrombocytopenia, Grade 3, dip_amp, myelosuppression
  - WEEKNUM, CTCAE, per-patient worst, worst_plt, NEUT, WBC
  - IRd, Rd, NDMM, RRMM, lenalidomide, ixazomib
load_cost: medium   # includes calibration table and dip model code
---

## Level 2 — Instructions

### Critical Definition: Per-Patient Worst, Not Per-Cycle Snapshot
Published Grade 3 rates are **cumulative worst grade across all cycles**.
Validating against per-cycle Day-15 snapshots gives ~6× lower rates.

```python
# CORRECT — used in validate_data.py
worst_plt = plt_lb.groupby('USUBJID')['AVAL'].min()
grade3_rate = (worst_plt < 50).mean() * 100

# WRONG — gives per-observation rate
grade3_rate = (plt_lb[plt_lb['WEEKNUM']==3]['AVAL'] < 50).mean() * 100
```

### Within-Cycle Dip Model
PLT nadir at Day 15 (WEEKNUM=3, lenalidomide days 1-21):
```
dip = dip_amp × |Day1_val| × sin(π × wk_off / 28)
```
- wk_off=0 (Day 1): sin=0 → no dip
- wk_off=7 (Day 8): sin≈0.71 → partial dip  
- wk_off=14 (Day 15): sin=1.0 → full nadir

### Grade 3 Rate vs dip_amp: Log-Linear Relationship
Grade 3 rate grows **super-linearly** — use log-linear interpolation, not quadratic:

```
slope = (ln(G3_b) - ln(G3_a)) / (dip_b - dip_a)
target_dip = dip_a + (ln(target_G3) - ln(G3_a)) / slope
```

Example (MM2 Rd, seed=42):
- dip=0.44 → 10.2%, dip=0.52 → 24.9%
- For target 14%: slope=(ln24.9−ln10.2)/(0.52−0.44)=11.2 → dip=0.47 ✓

### Biological Calibration Principles
| Factor | NDMM vs RRMM | IRd vs Rd |
|--------|-------------|-----------|
| PLT baseline | NDMM~200, RRMM~195 | Same |
| Dip cause | Lenalidomide days 1-21 | IRd adds Ixazomib thrombocytopenia |
| Grade 3 target | Lower in NDMM | Always higher for IRd |
| dip_amp needed | Higher for NDMM (higher baseline means bigger absolute dip needed) | Higher for IRd |

### After Any Change to the LB Generator
PLT is generated downstream of SPEP in the patient loop. Any change to the disease
marker section (more RNG calls, different probabilities) shifts the RNG state before
PLT, changing Grade 3 rates even if PLT code is unchanged.
→ Always re-validate PLT Grade 3 after ANY change to `_sim_trajectory()` or `make_lb()`.

---

## Level 3 — Resources

### Calibrated Dip Amplitudes (TOURMALINE, seeds 42/43)
| Arm | Study | dip_amp | Grade 3 sim | Target |
|-----|-------|---------|-------------|--------|
| IRd | MM2 (NDMM) | 0.45 | 25.4% | 25% |
| IRd | MM1 (RRMM) | 0.48 | 32.5% | 31% |
| Rd  | MM2 (NDMM) | 0.47 | 15.5% | 14% |
| Rd  | MM1 (RRMM) | 0.46 | 16.0% | 16% |

### Code Block (in `make_lb()`, within-cycle dip section)
```python
# Grade 3 empirical reference points for log-linear calibration:
#   IRd NDMM: 0.45→26.8%  |  IRd RRMM: 0.45→25.3%, 0.50→37.5% → 0.48@31%
#   Rd  NDMM: 0.44→10.2%, 0.52→24.9% → 0.47@14%
#   Rd  RRMM: 0.46→16.0%
if test in _MYELO_SUPP:
    if arm == "IRd":
        dip_amp = 0.45 if is_ndmm else 0.48
    else:
        dip_amp = 0.47 if is_ndmm else 0.46
    dip = dip_amp * abs(float(val)) * np.sin(np.pi * wk_off / 28.0)
    val_wc -= dip
```

### PLT Baseline Values
```python
# In _sim_trajectory(), PLT section:
if test == "PLT":
    base_mean = 200.0 if is_ndmm else 195.0
    # RRMM: prior therapy reduces bone marrow reserve → lower baseline
```

### File location
`scripts/generate_v2.py` → `make_lb()` → within-cycle dip block (~line 910)
`scripts/validate_data.py` → `write_validation_table()` → per-patient worst PLT check
