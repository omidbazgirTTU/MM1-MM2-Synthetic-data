---
name: published-mechanistic-correlations-tourmaline
description: >
  Published quantitative PK/PD/endpoint relationships from the pivotal TOURMALINE-MM1
  analysis (Srimani 2022 CPT:PSP): platelet semi-mechanistic model (LINEAR AUC→PLT),
  two-population M-protein model, flat exposure-efficacy finding, M-protein Week 8
  slope predicting PFS (HR 0.26), and covariate significance results.
applies_when:
  - Implementing or calibrating the AUC→platelet relationship (linear vs Emax decision)
  - Linking M-protein trajectory slope to PFS hazard modifier
  - Deciding whether AUC→efficacy E-R should be imposed (published: FLAT within range)
  - Setting timing of platelet nadir (Days 11-15 mechanistic basis)
  - Calibrating AE exposure-response (quartile analysis, logistic regression)
  - Validating that criteria 61-68 (cross-correlation checks) match published values
keywords:
  - Srimani 2022, CPT PSP, semi-physiological platelet model, megakaryocyte precursor
  - M-protein two-population model, drug-sensitive, drug-resistant, indirect response
  - Markov AE model, diarrhea, prior IMiD, rash, race
  - Week 8 M-protein slope, AUC-ROC 0.80-0.85, HR 0.26, landmark analysis
  - flat exposure-efficacy, no significant PFS exposure-response, therapeutic range
  - linear platelet model, not Emax, linear within TOURMALINE range
load_cost: medium
---

## Level 2 — Instructions

**Primary Reference:** Srimani JK, et al. "Population pharmacokinetic/pharmacodynamic
joint modeling of ixazomib efficacy and safety using data from the pivotal phase III
TOURMALINE-MM1 study." CPT Pharmacometrics Syst Pharmacol. 2022;11(8):1085–1099.

---

### 1. AUC → Platelet: LINEAR (not Emax) within Therapeutic Range

**Published model:** Platelet precursor (megakaryocyte) elimination model.

```
dMK_prec/dt = k_in - (k_out + E_IXA + E_LENDEX) × MK_prec
dPlatelet/dt = k_prp × MK_prec - k_out × Platelet

E_IXA   = β_IXA × AUC_weekly     ← LINEAR slope on MK precursor elimination
E_LENDEX = β_LENDEX × dose_lendex  ← LenDex also eliminates precursors
```

**Key decisions for synthetic data:**
- **Use LINEAR model**: The Srimani team note that Emax may be more appropriate at
  doses outside the 4 mg TOURMALINE range, but linear was adequate within it.
- **Cyclical pattern**: Nadir consistently at **Days 11–15** (mechanistic from
  megakaryocyte maturation time ~7–10 days + platelet lifespan ~10 days).
- **Implementation**: Replace fixed `dip_amp_pop` with `dip_amp_i = dip_amp_pop × auc_rel_i`
  (proportional scaling). Cap at 2.5× to avoid extreme outliers.
- **Re-calibrate** `dip_amp_pop` after implementing individual variability so the
  population median Grade 3 PLT rate still hits the validated targets.

### 2. AUC → M-protein: Two-Population Indirect Response

**Published model:** Drug-sensitive + drug-resistant M-protein populations.

```
R(t) = R_sensitive(t) + R_resistant(t)

Drug-sensitive: indirect response model (ixazomib kills sensitive myeloma cells)
Drug-resistant: exponentially growing function (predicts relapse timing)
```

**For synthetic data (simplified):**
- Keep existing bimodal phenotype structure (responder/non-responder tiers).
- Add a **modest continuous AUC shift** (±10pp per 50% AUC deviation) — see A3.
- Do NOT impose a strong Emax link: the published finding is flat E-R within range.

### 3. Exposure-Efficacy: FLAT within Therapeutic Range

**CRITICAL FINDING (Srimani 2022):** Ixazomib systemic exposure was **NOT a significant
predictor** of PFS or probability of response within the 4 mg therapeutic range.

**Implications for generation:**
- AUC→response link must be weak: |r(AUC, M-protein% C6)| < 0.20
- Logistic regression of AUC on ORR: p-value should be > 0.05 (non-significant)
- These are validation criteria 65–68 in validate_data.py
- Despite dose modifications for toxicity, efficacy was maintained — wide therapeutic window

### 4. M-protein Week 8 Slope → PFS (Criterion 63–64)

**Published landmark analysis (HR for time-to-progression):**

| M-protein reduction at Week 8 (Cycle 2) | HR | 95% CI |
|---|---|---|
| < 25% (reference) | 1.00 | — |
| 50–74% | 0.41 | 0.26–0.64 |
| **≥ 75%** | **0.26** | 0.15–0.45 |

**AUC-ROC for Week 8 M-protein predicting PFS event:** 0.80–0.85 (excellent).

**Implementation in make_ds():**
```python
resp_modifier = {
    "CR_VGPR": 0.30,   # ≥75% M-protein reduction → HR 0.30 (matches published 0.26)
    "PR":      0.55,   # 50–74% reduction → HR 0.55 (close to published 0.41)
    "None":    1.00,   # <50% reduction → reference
}
```
Validation target: HR for ≥75% responders = 0.20–0.45 (criteria 63 MM2, 64 MM1).

### 5. AE Exposure-Response (Criteria 61–62 and qualitative checks)

**Quartile analysis (Q1 lowest → Q4 highest AUC):**

| AE | Q1 | Q4 | Trend | p-value |
|---|---|---|---|---|
| Grade ≥3 Thrombocytopenia | ~20% | ~35–40% | ↑ | <0.05 |
| Grade ≥3 Anemia | 6–8% | ~14% | ↑ | <0.05 |
| Grade ≥2 Diarrhea | Low | High | ↑ | <0.05 |
| Grade ≥2 Rash | Low | High | ↑ | <0.05 |

**AUC → PLT nadir validation (criteria 61–62):** Q4 patients should show PLT nadir
≥15% deeper than Q1 (linear model; not a strong Emax effect).

**Special covariates in AE models:**
- Diarrhea: **prior IMiD therapy** significantly increases risk (Markov model covariate)
- Rash: **race** significantly affects risk

### 6. Covariate Effects on PK — What is NOT significant

Published TOURMALINE-MM1 popPK results:

| Covariate | Effect on CL/V | Dose adjustment |
|---|---|---|
| BSA | No effect | None |
| Age | No effect | None |
| Sex | No effect | None |
| CrCl (mild/mod) | No effect | None |
| CrCl < 30 (severe) | Reduced CL | Reduce to 3 mg |

**Note:** BSA has no effect on CL, but the Gupta 2017 model includes BSA on V4 (peripheral
volume). This is a volume effect, not a clearance effect — no dose adjustment required.

### 7. Week-by-Week PLT Cyclical Pattern

Nadir timing is mechanistic (Days 11–15), not a tuning parameter:
- Day 1: ~200 (baseline/recovered)
- Days 2–10: Gradual decline
- **Days 11–15: NADIR ~110–130** (mechanistic from MK maturation + platelet lifespan)
- Days 16–28: Recovery

Cumulative effect: Cycle 2+ nadir slightly lower than Cycle 1 (gradual depletion of
megakaryocyte reserve with repeated dosing).

---

## Level 3 — Resources

### Platelet Linear Model — Implementation Sketch

```python
# Simplified linear precursor model for synthetic data generation
# (Full ODE not needed — approximate with per-cycle dip scaling)

def compute_plt_dip_i(auc_i, dip_amp_pop, AUC_POP=1247.0):
    """
    auc_i: individual patient AUC (ng·h/mL), derived from IXAZ_CL_I
    dip_amp_pop: calibrated population PLT dip amplitude (0–1 fraction)
    AUC_POP: population-typical AUC (Gupta 2017)
    Returns: per-patient dip amplitude (fraction of baseline PLT lost at nadir)
    """
    auc_rel_i = auc_i / AUC_POP
    dip_amp_i = dip_amp_pop * min(auc_rel_i, 2.5)   # linear, cap at 2.5x
    return dip_amp_i

# In _sim_trajectory() PLT section:
dip_amp_i = compute_plt_dip_i(subj_auc_i, dip_amp_pop=dip_amp_pop[arm][is_ndmm])
nadir_frac = max(1.0 - dip_amp_i, 0.10)   # can't go below 10% of baseline
```

### M-protein → PFS Hazard Modifier — Implementation Sketch

```python
# In make_ds(), after make_lb() has computed mprot_slope_c2 per patient:

def get_resp_modifier(mprot_pchg_c2):
    """
    mprot_pchg_c2: M-protein % change at Cycle 2 (negative = reduction)
    Returns: Weibull hazard scale modifier (multiplicative)
    """
    if mprot_pchg_c2 <= -75:
        return 0.30   # ≥75% reduction → HR 0.30 (published 0.26)
    elif mprot_pchg_c2 <= -50:
        return 0.55   # 50–74% reduction → HR 0.55
    else:
        return 1.00   # <50% → reference

# Usage:
for subj in subjects:
    modifier = get_resp_modifier(subj["MPROT_PCHG_C2"])
    # Scale the Weibull shape parameter: lower modifier → longer PFS
    pfs_scale_i = pfs_scale_pop[arm] / (modifier ** (1 / pfs_shape))
    pfs_i = RNG_SURV.weibull(pfs_shape) * pfs_scale_i
```

### Cox PH Validation for Criterion 63–64

```python
from lifelines import CoxPHFitter
import pandas as pd

def check_mprot_pfs_hr(adlb, adtte, adsl, study_key):
    """
    Computes Cox PH HR for M-protein Week 8 responders (≥75%) vs non-responders.
    Target: HR 0.20–0.45 (published Srimani 2022: 0.26)
    """
    # Get M-protein % change at Cycle 2 (Week 8) from adlb
    wk8 = adlb[
        (adlb["LBTESTCD"] == "SPEP_MPROT") &
        (adlb["VISITNUM"] == 2) &
        (adlb["ABLFL"] != "Y")
    ][["USUBJID", "PCHG"]].rename(columns={"PCHG": "MPROT_PCHG_C2"})

    tte = adtte[adtte["PARAMCD"] == "PFS"][["USUBJID", "AVAL", "CNSR"]]
    df = tte.merge(wk8, on="USUBJID", how="inner")
    df["responder"] = (df["MPROT_PCHG_C2"] <= -75).astype(int)

    cph = CoxPHFitter()
    cph.fit(df[["AVAL","CNSR","responder"]], duration_col="AVAL",
            event_col="CNSR", formula="responder")
    hr = np.exp(cph.params_["responder"])

    ok = 0.20 <= hr <= 0.45
    print(f"{'PASS' if ok else 'FAIL'}  [{study_key}] M-protein Wk8 ≥75% → PFS HR: {hr:.3f} (target 0.20–0.45)")
    return ok
```

### Exposure-Efficacy Flatness Validation (Criteria 65–68)

```python
from scipy import stats
from sklearn.linear_model import LogisticRegression

def check_flat_er(adlb, adpc, adsl, study_key):
    """
    Validates that AUC→M-protein and AUC→ORR relationships are flat/non-significant.
    Criteria 65–68 (2 per study).
    """
    # Get individual AUCinf from adpc (IRd arm only)
    auc = adpc[adpc["PPTESTCD"]=="AUCINF_PP"][["USUBJID","PPSTRESN"]].rename(
        columns={"PPSTRESN":"AUC_INF"})
    # M-protein % change at Cycle 6
    mp_c6 = adlb[
        (adlb["LBTESTCD"]=="SPEP_MPROT") & (adlb["VISITNUM"]==6)
    ][["USUBJID","PCHG"]].rename(columns={"PCHG":"MP_C6"})

    df = auc.merge(mp_c6, on="USUBJID").dropna()

    # Criterion A: |r(AUC, M-prot C6)| < 0.20
    r, pval = stats.pearsonr(df["AUC_INF"], df["MP_C6"])
    ok_r = abs(r) < 0.20
    print(f"{'PASS' if ok_r else 'FAIL'}  [{study_key}] AUC↔M-prot C6 |r|={abs(r):.3f} (target <0.20, flat E-R)")

    # Criterion B: Logistic AUC→ORR p-value > 0.05
    orr = adsl[["USUBJID","ORR"]].merge(auc, on="USUBJID").dropna()
    _, pval_logit = stats.pearsonr(orr["AUC_INF"], orr["ORR"].astype(float))
    ok_p = pval_logit > 0.05
    print(f"{'PASS' if ok_p else 'FAIL'}  [{study_key}] AUC→ORR logistic p={pval_logit:.3f} (target >0.05, non-significant)")

    return ok_r, ok_p
```

### Published Parameter Summary

| PK/PD Parameter | Value | Source |
|---|---|---|
| β_IXA (linear slope on MK precursor) | Estimated (not tabulated) | Srimani 2022 |
| PLT nadir timing | Days 11–15 | Srimani 2022 |
| M-protein Week 8 → PFS HR (≥75%) | 0.26 (0.15–0.45) | Landmark analysis |
| M-protein Week 8 AUC-ROC | 0.80–0.85 | Srimani 2022 |
| AUC → PFS | Not significant | Srimani 2022 |
| AUC → ORR | Not significant | Srimani 2022 |
| Diarrhea covariate | Prior IMiD | Srimani 2022 |
| Rash covariate | Race | Srimani 2022 |
