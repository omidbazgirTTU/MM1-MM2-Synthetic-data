---
name: synthetic-clinical-trial-response-calibration
description: >
  Calibration techniques for IMWG response rates (ORR, VGPR, CR) in simulated
  oncology trial data using a bimodal phenotype mixture model. Covers the three
  root causes of response rate underestimation and their fixes, including the
  nadir-override technique and CR/VGPR boundary noise control.
applies_when:
  - ORR, VGPR+, or CR+ rates are below published targets
  - Changing phenotype probabilities does not produce the expected rate change
  - CR patients are appearing as VGPR (or VGPR as non-responders)
  - Best-response waterfall counts differ from expected mixture fractions
  - Recalibrating after any change to the disease-marker trajectory section
keywords:
  - VGPR, CR, ORR, response rate, waterfall, PCHG, SPEP_MPROT
  - nadir override, bimodal phenotype, resp_rate, _sample_resp_phenotype
  - best response, IMWG, -90%, -99%, miss_rate
load_cost: medium   # includes code snippets; load when actively calibrating
---

## Level 2 — Instructions

### Root Causes of Response Rate Underestimation

**1. Organic trajectory doesn't reach nadir in time**
With exponential decay `k=0.30`, a patient needs ~10 cycles to reach 90% of nadir.
Patients with n_cycles < 10 are misclassified even if assigned VGPR/CR phenotype.

→ **Fix**: Insert a guaranteed nadir observation at `min(n_cycles-1, 9)` for ALL
responders (resp_rate ≥ 0.50), including PR. Threshold: `n_cycles ≥ 2`.

**2. MAR miss_rate removes the nadir observation**
At miss_rate=0.08, 8% of VGPR/CR patients lose their override observation per cycle,
reducing VGPR+ by ~8 × (VGPR+ fraction) percentage points.

→ **Fix**: Set `miss_rate = 0.03` for `_DISEASE_BM` and `B2MG`. Disease-biomarker
visits (SPEP, Ig, UPEP) are the primary efficacy endpoint — near-perfect attendance.

**3. CR/VGPR boundary noise contamination**
CR tier range `[0.990, 1.000)` places patients with resp_rate ≈ 0.990 right at the
-99% PCHG threshold. Noise pushes ~50% of these borderline patients above -99% → VGPR.

→ **Fix**: Start CR tier at **0.993** (not 0.990). All CR patients then have nadir
≤ 0.007×base; 3σ upward noise still lands below -99% → correctly classified as CR.

### Decision Guide
| Symptom | Check first |
|---------|-------------|
| Both VGPR and CR low | miss_rate and n_cycles threshold |
| CR low but VGPR OK | CR tier boundary (0.990 → 0.993) |
| ORR low but VGPR OK | PR patients not reaching -50% → extend override to resp_rate ≥ 0.50 |
| Rate lower than mixture prob implies | Post-nadir floor pushing best min back up |
| Rates inconsistent between runs | RNG state contamination (see rng-management skill) |

### Validation Method (how best-response is computed)
```python
mp = adlb[(adlb['PARAMCD']=='SPEP_MPROT') & adlb['AVAL'].notna()]
bl   = mp[mp['EPOCH']=='BASELINE'].groupby('USUBJID')['AVAL'].first()
best = mp[mp['EPOCH']!='BASELINE'].groupby('USUBJID')['AVAL'].min()
pchg = ((best[common] - bl[common]) / bl[common] * 100)

orr_rate  = (pchg <= -50).mean() * 100
vgpr_rate = (pchg <= -90).mean() * 100   # VGPR or better
cr_rate   = (pchg <= -99).mean() * 100   # CR or better
```
Denominator = N patients with any non-baseline SPEP observation (excludes n_cycles=1).

### Phenotype Probability Calibration Rules
1. Published target → inflate by 1/0.97 (miss_rate loss)
2. Inflate further by N_total / N_waterfall (patients lost to n_cycles=1)
3. CR tier: set p_cr so that p_cr × 0.97 ≈ target_CR%; p_vgpr fills remainder to VGPR+ target

### Calibrated Values (TOURMALINE, seeds 42/43)
| Arm | Study | p_nonresp | p_pr | p_vgpr | p_cr* |
|-----|-------|-----------|------|--------|-------|
| IRd | MM2   | 0.15 | 0.22 | 0.35 | 0.28 |
| Rd  | MM2   | 0.22 | 0.18 | 0.43 | 0.17 |
| IRd | MM1   | 0.20 | 0.30 | 0.38 | 0.12 |
| Rd  | MM1   | 0.28 | 0.33 | 0.32 | 0.07 |

*p_cr = 1 − p_nonresp − p_pr − p_vgpr

---

## Level 3 — Resources

### Nadir Override Code (in `_sim_trajectory()`, after main for-t loop)
```python
if resp_rate >= 0.50 and n_cycles >= 2:
    nadir_t = min(n_cycles - 1, 9)   # cycle 10 or last available

    if resp_rate >= 0.90:
        # VGPR/CR: tight noise preserves -90%/-99% boundaries
        # CR patients (nadir≈0.005×base): std≈0.00085×base → deterministic
        # VGPR patients (nadir≈0.055×base): std≈0.0033×base → stays above -99%
        nadir_noise_std = nadir * 0.05 + base * 0.0005
        nadir_val = max(1e-6, nadir + RNG.normal(0, nadir_noise_std))
        vals[nadir_t] = nadir_val
        # Floor post-nadir: prevent noise creating false minima below nadir_val
        for t2 in range(nadir_t + 1, n_cycles):
            if vals[t2] < nadir_val:
                vals[t2] = nadir_val * RNG.uniform(1.001, 1.02)
    else:
        # PR: moderate noise matching organic trajectory (base×4%)
        nadir_val = max(1e-6, nadir + RNG.normal(0, base * 0.04))
        vals[nadir_t] = nadir_val
        # Do NOT floor post-nadir for PR — relapse is clinically expected
```

### CR Tier Definition (in `_sample_resp_phenotype()`)
```python
# CR/sCR tier — starts at 0.993 NOT 0.990
# Margin: resp_rate=0.993 → nadir=0.007×base
#         noise_std=0.007×0.05+0.0005×base=0.00085×base
#         3σ upward → 0.00955×base → PCHG=-99.04% → CR ✓
return float(RNG.uniform(0.993, 1.000))
```

### miss_rate Block (in `make_lb()` missingness section)
```python
if is_baseline:
    miss_rate = 0.03
elif test in _DISEASE_BM or test == "B2MG":
    miss_rate = 0.03   # key efficacy endpoint — near-perfect attendance
elif wk_off == 0:
    miss_rate = 0.05   # Day 1 visit
else:
    miss_rate = 0.10   # Day 8/15 intermediate visits
```

### File locations
- Phenotype sampler: `scripts/generate_v2.py` → `_sample_resp_phenotype()` (~line 570)
- Trajectory + nadir override: `scripts/generate_v2.py` → `_sim_trajectory()` (~line 614)
- Missingness block: `scripts/generate_v2.py` → `make_lb()` (~line 927)
- Validator: `scripts/validate_data.py` → `plot_mprotein()` + `write_validation_table()`
