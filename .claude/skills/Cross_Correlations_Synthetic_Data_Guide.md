---
name: cross-correlations-synthetic-data-guide
description: >
  Implementation patterns for correlated baseline covariates (7×7 MVN Cholesky),
  OMEGA matrix for PK etas, PK→PD mechanistic equations, and validation strategies
  for cross-correlations in the TOURMALINE synthetic dataset.
applies_when:
  - Implementing or debugging baseline covariate correlations (Age↔CrCl, Weight↔BSA, M-protein↔HGB, PLT↔HGB)
  - Adding OMEGA matrix Cholesky decomposition for correlated PK etas (CL, V2, V4)
  - Implementing mechanistic PK→PD equations (Emax M-protein kill, linear PLT precursor model)
  - Verifying correlation structure after any baseline generation change
  - Adding IC50_i / EC50_plt_i correlated PD heterogeneity
keywords:
  - MVN, multivariate normal, Cholesky, OMEGA, correlation matrix, covariance
  - age, CrCl, BSA, weight, M-protein, HGB, PLT, baseline covariates
  - Emax, IC50, k_kill, EC50_plt, indirect response, copula
  - IXAZ_CL_I, auc_rel, patient-level PK-PD link
load_cost: medium
---

## Level 2 — Instructions

### Why Cross-Correlations Are Required
Generating variables independently produces unrealistic data that:
1. Fails to capture physiological relationships (old people have lower CrCl)
2. Breaks covariate effects in the ML/PK framework (BSA → V4 requires realistic BSA range)
3. Misrepresents disease biology (high M-protein → anemia, thrombocytopenia)

### Critical Correlations — MUST PRESERVE

#### Baseline Covariate Correlations (7×7 MVN)

| Pair | r | Priority | Justification |
|------|---|----------|---------------|
| Age ↔ CrCl | **−0.45** | MUST | Renal function declines with age (Cockcroft-Gault mechanism) |
| Weight ↔ BSA | **+0.85** | MUST | Formula-based: BSA = √(ht × wt / 3600) |
| M-protein ↔ HGB | **−0.30** | MUST | Tumor burden displaces red cell production |
| M-protein ↔ BM plasma cells | +0.70 | MUST | Direct production relationship |
| PLT ↔ HGB | **+0.25** | Should | Shared bone marrow suppression |
| Age ↔ Weight | +0.10 | Optional | Weak association |
| CrCl ↔ Weight | +0.20 | Optional | Body size effect |

**Full 7×7 correlation matrix (indices: Age=0, CrCl=1, Weight=2, BSA=3, M-prot=4, PLT=5, HGB=6):**
```
           Age   CrCl  Weight  BSA  M-prot  Plt    Hgb
Age     [1.00 -0.45  0.10   0.05 -0.05  -0.10 -0.15]
CrCl   [-0.45  1.00  0.20   0.15  0.05   0.15  0.20]
Weight  [0.10  0.20  1.00   0.85 -0.10   0.10  0.15]
BSA     [0.05  0.15  0.85   1.00 -0.08   0.08  0.12]
M-prot [-0.05  0.05 -0.10  -0.08  1.00  -0.20 -0.30]
Plt    [-0.10  0.15  0.10   0.08 -0.20   1.00  0.25]
Hgb    [-0.15  0.20  0.15   0.12 -0.30   0.25  1.00]
```

Study-specific means and SDs:
```
MM2 (NDMM): age=73/7, crcl=65/30, weight=72/17, bsa=1.82/0.24, mprot=3.8/2.3, plt=203/84, hgb=9.85/1.85
MM1 (RRMM): age=66/10, crcl=75/32, weight=75/18, bsa=1.85/0.25, mprot=3.2/2.1, plt=222/81, hgb=10.15/1.75
```

#### PK Parameter Correlations (OMEGA matrix — assumed, not published for Ixazomib)

Off-diagonals represent covariances between log-normal etas:
- Cov(CL, V2) = 0.35 × 0.30 × 0.40 = 0.042 → r_CL,V2 ≈ 0.40
- Cov(CL, V4) = 0.35 × 0.45 × 0.29 = 0.046 → r_CL,V4 ≈ 0.29
- Cov(V2, V4) = 0.30 × 0.45 × 0.20 = 0.027 → r_V2,V4 ≈ 0.20

### Implementation Strategy

**Preferred method: MVN Cholesky (Method 2)** — more numerically stable than
`scipy.stats.multivariate_normal.rvs()`, and integrates with the existing
`RNG = np.random.default_rng(seed)` pattern in the generators.

```python
cov = np.outer(sds, sds) * CORR_MATRIX
L = np.linalg.cholesky(cov)
z = RNG.standard_normal((n, 7))
samples = (L @ z.T).T + means   # shape (n, 7)
```

**Sequential conditional sampling (Method 3)** is an alternative for causal chains
but harder to control exact correlations. Prefer MVN for the 7-variable baseline block.

### PD Heterogeneity — Correlated IC50/EC50

Resistant patients tend to have **both** high IC50 (M-protein insensitive to ixazomib)
AND high EC50_plt (platelets also insensitive). Implemented via 2×2 OMEGA_PD:

```
OMEGA_PD = [[0.40², 0.40×0.40×0.35],   r_IC50,EC50 = 0.40
            [0.40×0.40×0.35, 0.35²]]
```

After Cholesky draw: `IC50_i = 200 × exp(eta_pd[:,0])`, `EC50_plt_i = 60 × exp(eta_pd[:,1])`

### Patient-Level PK-PD Link (IXAZ_CL_I)

The current generators have a disconnect: `generate_v2.py` draws `auc_rel` independently,
while `generate_pk_v2.py` samples its own `CL_i`. Fix:

1. `make_dm()` computes `IXAZ_CL_I` (per-patient CL) from the OMEGA Cholesky draw and
   stores it as a column in DM/ADSL.
2. `make_lb()` reads `IXAZ_CL_I` from the subject dict to compute `auc_i = (F×dose×1000)/cl_i`.
3. `sample_ixaz()` in `generate_pk_v2.py` reads `IXAZ_CL_I` from `adam_adsl.csv` instead
   of re-sampling CL, ensuring both generators use the same per-patient CL.

---

## Level 3 — Resources

### Full MVN Baseline Implementation (make_dm() pattern)

```python
import numpy as np

BASELINE_CORR = np.array([
    # Age    CrCl   Weight  BSA    M-prot  Plt    Hgb
    [ 1.00, -0.45,  0.10,  0.05, -0.05, -0.10, -0.15],
    [-0.45,  1.00,  0.20,  0.15,  0.05,  0.15,  0.20],
    [ 0.10,  0.20,  1.00,  0.85, -0.10,  0.10,  0.15],
    [ 0.05,  0.15,  0.85,  1.00, -0.08,  0.08,  0.12],
    [-0.05,  0.05, -0.10, -0.08,  1.00, -0.20, -0.30],
    [-0.10,  0.15,  0.10,  0.08, -0.20,  1.00,  0.25],
    [-0.15,  0.20,  0.15,  0.12, -0.30,  0.25,  1.00],
])

STUDY_PARAMS = {
    "MM2": dict(
        means=np.array([73, 65, 72, 1.82, 3.8, 203, 9.85]),
        sds  =np.array([ 7, 30, 17, 0.24, 2.3,  84, 1.85]),
        bounds=[(23,93),(20,220),(40,140),(1.3,2.6),(0.5,10),(50,500),(6,15)]
    ),
    "MM1": dict(
        means=np.array([66, 75, 75, 1.85, 3.2, 222, 10.15]),
        sds  =np.array([10, 32, 18, 0.25, 2.1,  81,  1.75]),
        bounds=[(23,93),(20,220),(40,140),(1.3,2.6),(0.5,10),(50,500),(6,15)]
    ),
}

def draw_correlated_baseline(n, study, RNG):
    p = STUDY_PARAMS[study]
    cov = np.outer(p["sds"], p["sds"]) * BASELINE_CORR
    L = np.linalg.cholesky(cov)
    z = RNG.standard_normal((n, 7))
    raw = (L @ z.T).T + p["means"]
    for i, (lo, hi) in enumerate(p["bounds"]):
        raw[:, i] = np.clip(raw[:, i], lo, hi)
    # Returns columns: age, crcl, weight, bsa, mprot, plt, hgb
    return raw
```

### OMEGA Cholesky for PK Etas (sample_ixaz() pattern)

```python
OMEGA = np.array([
    [0.1225, 0.042,  0.046],   # CL: ω²=0.35², off-diag = covariances
    [0.042,  0.09,   0.027],   # V2: ω²=0.30²
    [0.046,  0.027,  0.2025],  # V4: ω²=0.45²
])
L_omega = np.linalg.cholesky(OMEGA)

def sample_ixaz_etas(n, RNG):
    eta_raw = RNG.standard_normal((n, 3))
    etas = (L_omega @ eta_raw.T).T   # shape (n, 3): [η_CL, η_V2, η_V4]
    return etas

# Usage:
etas = sample_ixaz_etas(n, RNG)
cl_i = CL_pop * np.exp(etas[:, 0]) * ddi_mult          # CYP3A4 modifier applied here
v2_i = V2_pop * np.exp(etas[:, 1])
v4_i = v4_bsa  * np.exp(etas[:, 2])   # v4_bsa already has BSA covariate applied
```

### Correlated PD Heterogeneity (IC50_i, EC50_plt_i)

```python
OMEGA_PD = np.array([
    [0.40**2,          0.4 * 0.40 * 0.35],
    [0.4 * 0.40 * 0.35, 0.35**2],
])
L_pd = np.linalg.cholesky(OMEGA_PD)

eta_pd = (L_pd @ RNG.standard_normal((2, n))).T   # shape (n, 2)
IC50_i     = 200.0 * np.exp(eta_pd[:, 0])   # ng·h/mL, population IC50 for M-protein
EC50_plt_i =  60.0 * np.exp(eta_pd[:, 1])   # ng/mL, population EC50 for platelet
```

### Correlation Verification (post-generation)

```python
import pandas as pd

def verify_baseline_corr(adsl, study):
    targets = [
        ("AGE", "BASE_CREACL",    -0.45, 0.10, "Age ↔ CrCl"),
        ("WEIGHT", "BSA",          0.85, 0.05, "Weight ↔ BSA"),
        ("BASE_SPEP_MPROT","BASE_HGB", -0.30, 0.10, "M-prot ↔ HGB"),
        ("BASE_PLT","BASE_HGB",    0.25, 0.10, "PLT ↔ HGB"),
    ]
    print(f"\n--- Baseline Correlation Verification [{study}] ---")
    for c1, c2, tgt, tol, label in targets:
        r = adsl[[c1, c2]].corr().iloc[0,1]
        ok = abs(r - tgt) <= tol
        print(f"{'PASS' if ok else 'FAIL'}  {label}: r={r:.3f} (target {tgt:+.2f} ±{tol})")
```

### Physiological Bounds Reference

| Variable | Min | Max | Unit |
|----------|-----|-----|------|
| Age | 23 | 93 | years |
| CrCl | 20 | 220 | mL/min |
| Weight | 40 | 140 | kg |
| BSA | 1.3 | 2.6 | m² |
| M-protein | 0.5 | 10.0 | g/dL |
| PLT | 50 | 500 | ×10⁹/L |
| HGB | 6.0 | 15.0 | g/dL |

### File Locations

| File | Role |
|------|------|
| `scripts/generate_v2.py` → `make_dm()` ~L290 | MVN draw + IXAZ_CL_I computation |
| `scripts/generate_v2.py` → `make_lb()` ~L869 | Reads IXAZ_CL_I → auc_rel |
| `scripts/generate_v2.py` → `_sim_trajectory()` ~L740 | HGB/PLT baseline adjusted by M-protein |
| `scripts/generate_pk_v2.py` → `sample_ixaz()` ~L170 | Reads IXAZ_CL_I from ADSL |
| `scripts/validate_data.py` → `check_cross_correlations()` | Criteria 49–68 |
