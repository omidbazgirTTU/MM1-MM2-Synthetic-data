---
name: validation-framework
description: Defines how validation criteria are structured, how tolerances are set, what the automated validator checks, and what additional checks the agent sign-off layer adds. Used by the QSP Scientist and Literature Search agents to agree on end conditions.
metadata:
  type: rules
---

# Validation Framework

## Two-Layer Validation Architecture

```
Layer 1 — Automated Quantitative Criteria
  engine/validate_trial.py reads config.yaml → computes all metrics → emits pass/fail table
  Machine-readable output: trials/<name>/outputs/validation_summary.csv
  Criterion count: config-driven (typically 60–80 per trial)

Layer 2 — Agent Sign-Off (qualitative + mechanistic)
  All 4 agents independently review validation output and generated data
  Output: trials/<name>/outputs/agent_signoff.md
  End condition: all 4 agents APPROVED
```

---

## Criterion Categories and Tolerances

### Category 1 — Survival Endpoints
**Tolerance**: ±5% relative (tightest — primary endpoint)

| Metric | How computed | Pass condition |
|--------|-------------|----------------|
| PFS median per arm | KM estimator on adam_adtte | abs(sim − target) / target ≤ 0.05 |
| OS median per arm | KM estimator on adam_adtte | abs(sim − target) / target ≤ 0.05 |
| PFS HR (active vs control) | lifelines CoxPHFitter | HR within published 95% CI |
| Censoring rate | CNSR=1 fraction | within ±10 pp of target |

### Category 2 — Enrollment
**Tolerance**: Exact match (0% tolerance)

| Metric | How computed | Pass condition |
|--------|-------------|----------------|
| N per arm | sdtm_dm row count by ARMCD | sim == target exactly |
| Total N | sum | sim == target exactly |

### Category 3 — Demographics
**Tolerance**: ±10% relative

| Metric | Source column | Pass condition |
|--------|--------------|----------------|
| Median age | adam_adsl AGE | abs(sim − target) / target ≤ 0.10 |
| % female | adam_adsl SEX='F' | abs(sim − target) ≤ 0.10 (absolute) |
| ECOG distribution | adam_adsl ECOGBL | each level within ±10 pp |

### Category 4 — Disease-Specific Staging
**Tolerance**: ±10% relative

| Example metrics by indication | |
|-------------------------------|--|
| MM: ISS Stage I/II/III fractions | adam_adsl ISS |
| CLL: Rai Stage 0/I–IV fractions | adam_adsl RAI |
| Solid tumor: ECOG 0 vs 1 | |

### Category 5 — Biomarker Distributions
**Tolerance**: ±15% relative on prevalence rates

| Metric | Pass condition |
|--------|----------------|
| Cytogenetic abnormality rates (del17p, t414, etc.) | ±15% |
| FLT3-ITD prevalence (RATIFY) | ±15% |
| BRCA mutation prevalence (SOLO-1) | ±15% |
| PD-L1 TPS distribution (KEYNOTE-189) | ±15% |
| CrCL ≤60 mL/min fraction | ±15% |

### Category 6 — Efficacy / Response Rates
**Tolerance**: ±15% relative

| Metric | How computed | Pass condition |
|--------|-------------|----------------|
| ORR (≥PR) | adam_adrs best response ≤−50% (or CR flag for AML) | ±15% |
| VGPR+ | adam_adrs best PCHG ≤−90% | ±15% |
| CR+ | adam_adrs best PCHG ≤−99% | ±15% |
| CR rate (AML) | blast% < 5% AND ANC ≥1.0 AND PLT ≥100 | ±15% |
| MRD negativity rate | if published | ±20% |

### Category 7 — PK NCA
**Tolerance**: ±15% relative

| Metric | Source | Pass condition |
|--------|--------|----------------|
| Cmax median (Cycle 1, primary drug) | sdtm_pp PPSTRESN | ±15% |
| AUCinf median | sdtm_pp | ±15% |
| IIV CV% on CL | sdtm_pp CL/F population SD / mean | ±10 pp absolute |

### Category 8 — Safety / Toxicity
**Tolerance**: ±20% relative

| Metric | How computed | Pass condition |
|--------|-------------|----------------|
| Grade ≥3 PLT rate | sdtm_lb PLT worst < 50 per patient | ±20% |
| Grade ≥3 ANC rate | sdtm_lb ANC worst < 0.5 per patient | ±20% |
| Grade ≥3 AE rate (drug-specific) | sdtm_ae max AETOXGR ≥3 per patient | ±20% |
| Dose modification rate | sdtm_ex EXDOSMOD='Y' fraction | ±20% |

### Category 9 — Baseline Covariate Cross-Correlations
**Tolerance**: ±0.10 absolute on Pearson r

| Criterion | Target | Tolerance |
|-----------|--------|-----------|
| r(Age, CrCL) | −0.45 | ±0.10 |
| r(Weight, BSA) | 0.85 | ±0.13 (wider — analytical BSA) |
| r(M-prot/ALC/SLD, HGB) | indication-specific | ±0.10 |
| r(PLT, HGB) | 0.25 (MM) | ±0.10 |

### Category 10 — PK NCA Cross-Correlations
**Tolerance**: r-based thresholds

| Criterion | Pass condition |
|-----------|----------------|
| r(Cmax, AUCinf) | > 0.35 (3-cmt model; > 0.70 for 1-cmt) |
| r(BSA, Cmax) | within target ± 0.13 |
| r(CrCL, AUC_Lenalidomide) | < −0.30 (renal covariate confirmed) |

### Category 11 — AUC → Toxicity Link
**Tolerance**: Directional test

| Criterion | Pass condition |
|-----------|----------------|
| Q4 vs Q1 AUC: PLT/ANC nadir depth | Q4 ≥15% deeper than Q1 |

### Category 12 — Biomarker → Survival
**Tolerance**: Published 95% CI

| Criterion | Pass condition |
|-----------|----------------|
| Cox HR: deep Cycle 6 responders vs non-responders | within published 95% CI (or 0.20–0.45 if CI not available) |

### Category 13 — Exposure-Efficacy Flatness
**Pass condition**: Both must be true within therapeutic range

| Criterion | Pass condition |
|-----------|----------------|
| \|r(AUC, response% at Cycle 6)\| | < 0.20 |
| AUC → ORR logistic regression p-value | > 0.05 (not significant) |

*Note: Only applies where flat E-R has been published (ixazomib, many checkpoint inhibitors).
If the drug has a known steep E-R (e.g., midostaurin at sub-therapeutic doses), this criterion
is omitted from the config.*

---

## VPC / GOF Checks (PK-Specific)

Generated by `engine/pk_vpc_gof.py`. These are not in the automated criterion table but
are reviewed by QSP Scientist during sign-off.

| Check | Pass criterion |
|-------|----------------|
| % observations within 5th–95th PI | ≥80% (all drugs, both arms) |
| GOF: r(PRED, OBS) | > 0.80 |
| CWRES distribution | N(0,1) by Shapiro-Wilk p > 0.05, mean ∈ [−0.5, 0.5] |
| CYP3A4 inhibitor patients: Cmax shift | Visibly elevated vs non-DDI patients |
| CrCL covariate: AUC vs CrCL slope | Negative slope with p < 0.05 |

---

## Recalibration Rules

When a criterion FAILs:

```
Step 1: Identify the failing metric and its current simulated value
Step 2: Identify which config parameter controls this metric
Step 3: Binary search within ±30% of published value
        - If metric is proportional to param: adjust linearly
        - If metric is non-linear (Grade 3 rate): use binary search
Step 4: Re-run generate + validate with updated param
Step 5: Confirm the change did NOT break previously passing criteria
        (run full validation suite, not just the failing criterion)
Step 6: Document in agent_signoff.md: [param], [old value], [new value], [rationale]
```

**One parameter at a time rule**: Never adjust more than one parameter per recalibration
iteration. Simultaneous multi-parameter changes make root cause analysis impossible.

**Seed invariance rule**: Never change a seed. Seeds are frozen after first validated run.
If a seed produces an unlucky random draw (e.g., t(4;14) p=0.090 → 10.1%), adjust the
corresponding probability parameter, not the seed.

**Cross-criterion contamination check**: When any parameter changes, re-run ALL criteria
to verify no regressions. The global RNG means some parameter changes can indirectly
shift criteria that were previously passing.

---

## Agent Sign-Off Log Format

`trials/<name>/outputs/agent_signoff.md` must contain:

```markdown
# Agent Sign-Off Log — TOURMALINE-MM1

## Literature Search Agent
**Status**: APPROVED
**Date**: 2026-06-08
**Evidence sources confirmed**: 8/8 validation targets grounded in published source
**Gaps**: None
**Notes**: —

## QSP Scientist Agent
**Status**: APPROVED
**Date**: 2026-06-08
**Mechanistic checks**:
- [x] AUC–PLT slope positive and linear: confirmed (Q4 27% deeper than Q1)
- [x] IC50 saturation pattern reproduced: |r(AUC, M-prot)| = 0.10 < 0.20 ✓
- [x] AR(1) on M-protein: ρ=0.60 confirmed in residual autocorrelation plot
- [x] VPC: 89.8% within 5th–95th PI ✓
**Calibration changes made**:
- T(4;14) probability: 0.090 → 0.075 (MM1 only; RNG drift from pk_series addition)

## Oncologist Agent
**Status**: APPROVED
**Date**: 2026-06-08
**Clinical plausibility**:
- [x] PLT nadir Day 11–15: confirmed
- [x] M-protein U-shape trajectory in progressors: confirmed
- [x] ALC redistribution (ibrutinib): N/A for this trial
**Concerns raised and resolved**: None

## Drug Developer Agent
**Status**: APPROVED
**Date**: 2026-06-08
**Protocol fidelity**:
- [x] Ixazomib Day 1/8/15 schedule: confirmed
- [x] 4mg→3mg→2.3mg reduction steps: confirmed
- [x] CYP3A4 inhibitor rate 8%: confirmed
**Concerns raised and resolved**: None

## FINAL DECISION: ALL AGENTS APPROVED — TRIAL COMPLETE
```

---

## Total Expected Criterion Count Per Trial

| Category | Count per study | Applies to |
|----------|----------------|------------|
| Survival | 4–6 | All trials |
| Enrollment | 2 | All trials |
| Demographics | 4 | All trials |
| Staging | 3–6 | Indication-specific |
| Biomarkers | 3–8 | Indication-specific |
| Efficacy | 6–12 | All trials |
| PK NCA | 2–4 | All trials |
| Safety | 2–6 | All trials |
| Covariate cross-correlations | 4–6 | All trials |
| PK NCA cross-correlations | 2–4 | All trials |
| AUC → toxicity | 1–2 | All trials |
| Biomarker → survival | 1–2 | All trials |
| Exposure-efficacy flatness | 0–4 | Where flat E-R published |
| **Total** | **~40–70** | Per trial |
