---
name: trial-config-schema
description: Canonical YAML schema for a synthetic clinical trial config. Every trial in the multi-agent framework must have a config.yaml that conforms to this schema. Used by the engine to generate SDTM/ADaM data reproducibly.
metadata:
  type: reference
---

# Trial Config Schema (YAML)

Every trial config lives at `trials/<trial_name>/config.yaml`. The engine reads this file
to generate all SDTM and ADaM domains. All parameters with a `[source]` annotation must
be populated by the Literature Search agent before generation.

---

## Full Annotated Schema

```yaml
# ─────────────────────────────────────────────────────────────────────────────
# TRIAL IDENTITY
# ─────────────────────────────────────────────────────────────────────────────
trial:
  name: TOURMALINE-MM1              # Canonical trial name (matches branch name)
  nct_id: NCT01564537               # ClinicalTrials.gov identifier
  phase: 3                          # Must be 3 for all trials in this project
  sponsor: Takeda                   # Sponsor company
  indication: RRMM                  # Short indication label
  icd10: C90.0                      # ICD-10 code for primary indication
  response_criteria: IMWG           # IMWG | RECIST1.1 | iwCLL | Lugano | AML-CR | RECIST+CNS
  primary_endpoint: PFS             # PFS | OS | EFS | CR_rate | DFS
  studyid: TML-MM1                  # Used in USUBJID construction: STUDYID-SITEID-SUBJECTID

# ─────────────────────────────────────────────────────────────────────────────
# REPRODUCIBILITY
# ─────────────────────────────────────────────────────────────────────────────
seeds:
  main: 43                          # Global RNG seed (integer, never change after validation)
  survival: 77                      # Survival-specific RNG seed (independent)
  mvn: 1043                         # MVN baseline covariate RNG (main + 1000)
  pk: 5043                          # PK generator RNG (main + 5000)

# ─────────────────────────────────────────────────────────────────────────────
# STUDY DESIGN
# ─────────────────────────────────────────────────────────────────────────────
study:
  n_cycles: 26                      # Maximum cycles per patient
  cycle_duration_days: 28           # Days per cycle (21 or 28 typically)
  censoring_cycle: 26               # Administrative censoring at this cycle
  pk_substudy_fraction: 0.20        # Fraction of subjects with dense PK sampling
  dense_cycles: [1, 3]              # Cycle numbers with dense PK sampling
  sparse_timepoints_h: [0, 4, 24]   # Sparse PK sampling timepoints (hours post-dose)
  dense_timepoints_h: [0, 0.5, 1, 2, 4, 8, 24, 48, 72, 168]  # Dense timepoints

# ─────────────────────────────────────────────────────────────────────────────
# ARMS
# ─────────────────────────────────────────────────────────────────────────────
arms:
  - name: IRd                       # Arm name (used throughout SDTM/ADaM)
    armcd: IRD                      # CDISC ARMCD (≤8 chars, uppercase)
    n: 360                          # [source: published Table 1]
    is_active: true                 # true = contains investigational drug
    treatment:
      - drug: IXAZOMIB
        dose: 4.0
        dose_unit: mg
        route: ORAL
        cycle_days: [1, 8, 15]      # Days within cycle (1-indexed)
      - drug: LENALIDOMIDE
        dose: 25.0
        dose_unit: mg
        route: ORAL
        cycle_days_range: [1, 21]   # Continuous dosing days 1 through 21
      - drug: DEXAMETHASONE
        dose: 40.0
        dose_unit: mg
        route: ORAL
        cycle_days: [1, 8, 15, 22]

  - name: Rd
    armcd: RD
    n: 362                          # [source: published Table 1]
    is_active: false
    treatment:
      - drug: PLACEBO
        dose: 0.0
        dose_unit: mg
        route: ORAL
        cycle_days: [1, 8, 15]
      - drug: LENALIDOMIDE
        dose: 25.0
        dose_unit: mg
        route: ORAL
        cycle_days_range: [1, 21]
      - drug: DEXAMETHASONE
        dose: 40.0
        dose_unit: mg
        route: ORAL
        cycle_days: [1, 8, 15, 22]

# Step-up dosing example (venetoclax):
# step_up_dosing:
#   drug: VENETOCLAX
#   schedule:
#     - week: 1
#       dose: 20.0
#     - week: 2
#       dose: 50.0
#     - week: 3
#       dose: 100.0
#     - week: 4
#       dose: 200.0
#     - week: 5+
#       dose: 400.0

# ─────────────────────────────────────────────────────────────────────────────
# DOSE MODIFICATION RULES
# ─────────────────────────────────────────────────────────────────────────────
dose_modifications:
  - drug: IXAZOMIB
    mod_rate: 0.15                  # Overall dose modification rate [source: published safety table]
    hold_rules:
      - lab: PLT
        threshold: 75               # ×10⁹/L
        operator: "<"
        action: HOLD
      - grade: 3                    # CTCAE Grade ≥3 non-hematologic AE
        action: HOLD
    reduction_steps:
      - from: 4.0
        to: 3.0
        trigger: PLT_lt_30_or_grade3
      - from: 3.0
        to: 2.3
        trigger: repeat_toxicity
    discontinue_threshold:
      lab: PLT
      value: 25
      operator: "<"

# ─────────────────────────────────────────────────────────────────────────────
# PK MODELS
# ─────────────────────────────────────────────────────────────────────────────
pk_models:
  - drug: IXAZOMIB
    arms: [IRd]                     # Which arms this drug appears in
    model_type: 3cmt_oral           # 1cmt_oral | 2cmt_oral | 3cmt_oral | 1cmt_iv |
                                    # 2cmt_iv | tmdd | irreversible_btk | adc_twolevel
    params:                         # Typical population values [source: Gupta 2017]
      CL: 1.86                      # L/h
      V2: 7.0                       # L (central)
      Q3: 14.4                      # L/h
      V3: 87.3                      # L
      Q4: 0.60                      # L/h
      V4: 448.6                     # L
      Ka: 0.50                      # h⁻¹
      F: 0.58                       # bioavailability
      ALAG: 0.15                    # h, absorption lag
    iiv:                            # Log-normal IIV omegas (not CV%) [source: Gupta 2017]
      CL: 0.36
      V2: 0.28
      V4: 0.45
      Ka: 0.55
    iiv_cholesky:                   # Off-diagonal correlations for Cholesky OMEGA
      CL_V2: 0.30
      CL_V4: 0.20
      V2_V4: 0.25
    iov:
      Ka: 0.25                      # IOV on absorption rate
    covariates:
      - variable: BSA
        parameter: V4
        power: 0.70
        reference: 1.73             # m²
      - variable: CYP3A4_STRONG_INHIB
        parameter: CL
        multiplier: 0.55
      - variable: CYP3A4_STRONG_INDUC
        parameter: CL
        multiplier: 2.0
    lloq: 0.5                       # ng/mL
    units: ng/mL
    residual_error:
      proportional: 0.20
      additive: 0.50
    nca_reference:                  # Published NCA targets [source: Gupta 2017]
      Cmax_median: 41.0             # ng/mL
      AUCinf_median: 1247.0         # ng·h/mL
      t_half: 228.0                 # h

# ─────────────────────────────────────────────────────────────────────────────
# PD MODELS
# ─────────────────────────────────────────────────────────────────────────────
pd_models:
  - model: TwoPopulationMprotein   # Engine class name (see pd_modeling_guide.md)
    arms: [IRd]                    # Active drug arm only (Rd uses fallback drift)
    primary_drug: IXAZOMIB
    pk_input: Cp                   # What PK output feeds this PD model: Cp | AUC | occupancy
    params:                        # [source: Srimani 2022]
      IC50: 3.29                   # ng/mL
      Imax: 0.758
      k_R: 0.206                   # /wk
      Y_SS: 0.143                  # fraction
      k_L: 0.00951                 # /wk
    iiv:
      IC50: 0.42
      k_R: 0.81
      Y_SS: 1.55
    ar1_rho: 0.60                  # AR(1) autocorrelation on M-protein residuals
    biomarker_column: SPEP_MPROT   # Which sdtm_lb column is this PD endpoint

  - model: FribergMyelosuppression
    cell_type: PLT
    arms: [IRd, Rd]
    primary_drug: IXAZOMIB
    pk_input: AUC_cum              # Cumulative AUC drives progressive PLT nadir
    params:
      k_cum: 1.337e-5              # /(ng·h/mL) — cumulative AUC slope
      nadir_day: 13                # Expected nadir day within cycle
      dip_amp_IRd: 0.45
      dip_amp_Rd: 0.47
    iov_omega: 0.03                # Per-cycle PLT nadir IOV
    biomarker_column: PLT

# ─────────────────────────────────────────────────────────────────────────────
# POPULATION (DEMOGRAPHICS + BASELINE)
# ─────────────────────────────────────────────────────────────────────────────
population:
  ndmm: false                      # true = newly diagnosed, false = relapsed/refractory
  median_age: 66                   # [source: Table 1]
  age_sd: 9
  pct_female: 0.43                 # [source: Table 1]
  ecog_0: 0.55                     # Fraction ECOG 0
  ecog_1: 0.40
  ecog_2: 0.05

  # ISS staging probabilities
  iss_probs: [0.27, 0.40, 0.33]    # Stage I, II, III [source: Table 1]

  # Immunoglobulin type distribution
  ig_probs:
    IgG: 0.54
    IgA: 0.27
    IgM: 0.04
    other: 0.15

  # Cytogenetic abnormalities [source: Table 1]
  cytogenetics:
    DEL17P: 0.105
    T414: 0.063                    # t(4;14)
    T1416: 0.030
    T1420: 0.010
    GAIN1Q21: 0.350
    DEL1P32: 0.080
    AMP1Q: 0.080
    high_risk_threshold: 0.01      # Minor allele freq cutoff

  # MVN baseline covariate correlation matrix
  # Variables: [Age, CrCL, Weight, BSA, Mprot, PLT, HGB]
  mvn_means: [66, 65, 72, 1.82, 3.8, 203, 9.85]
  mvn_sds:   [9,  30, 17, 0.24, 2.3,  84, 1.85]
  mvn_correlations:              # Upper triangle only; matrix is symmetric
    Age_CrCL:  -0.45
    Age_Wt:    -0.10
    Age_BSA:   -0.10
    Age_Mprot:  0.10
    Age_PLT:   -0.10
    Age_HGB:   -0.15
    CrCL_Wt:    0.15
    CrCL_BSA:   0.15
    CrCL_Mprot:-0.15
    CrCL_PLT:   0.10
    CrCL_HGB:   0.10
    Wt_BSA:     0.85
    Wt_Mprot:   0.05
    Wt_PLT:     0.05
    Wt_HGB:     0.05
    BSA_Mprot:  0.05
    BSA_PLT:    0.05
    BSA_HGB:    0.05
    Mprot_PLT: -0.20
    Mprot_HGB: -0.30
    PLT_HGB:    0.25

  # CYP3A4 DDI concomitant medications
  cyp3a4_inhibitor_rate: 0.08
  cyp3a4_inducer_rate:   0.03

# ─────────────────────────────────────────────────────────────────────────────
# SURVIVAL
# ─────────────────────────────────────────────────────────────────────────────
survival:
  pfs:
    distribution: weibull
    median_active: 20.6            # months [source: Moreau 2016]
    median_control: 14.7           # months [source: Moreau 2016]
    hr_active_vs_control: 0.74     # [source: Moreau 2016]
    censoring_rate_target: 0.35    # ~35% censored for PFS
  os:
    distribution: weibull
    median_active: 53.6            # months [source: Kumar 2021]
    median_control: 51.6           # months
    censoring_rate_target: 0.55

  # M-protein response → PFS Gaussian copula
  mprot_pfs_copula:
    rho: -0.80                     # Spearman correlation
    arm: IRd                       # Apply only to active arm
    biomarker_cycle: 6             # Use M-protein % change at this cycle

# ─────────────────────────────────────────────────────────────────────────────
# VALIDATION TARGETS
# ─────────────────────────────────────────────────────────────────────────────
validation_targets:
  - metric: PFS_median_IRd
    value: 20.6
    unit: months
    tolerance: 0.05                # ±5% for primary survival endpoint
    source: "Moreau 2016 NEJM 374:1621"

  - metric: PFS_median_Rd
    value: 14.7
    unit: months
    tolerance: 0.05
    source: "Moreau 2016 NEJM 374:1621"

  - metric: ORR_IRd
    value: 78.0
    unit: pct
    tolerance: 0.15
    source: "Moreau 2016 Table 2"

  - metric: Grade3_PLT_IRd
    value: 31.0
    unit: pct
    tolerance: 0.20
    source: "Moreau 2016 Table S5"

  - metric: Ixazomib_Cmax_median
    value: 41.0
    unit: ng/mL
    tolerance: 0.15
    source: "Gupta 2017 Clin Pharmacokinet"

  # Cross-correlation criteria (generated automatically from mvn_correlations section)
  # Add any non-standard cross-correlation targets here:
  - metric: r_AUC_Mprot_C6
    value: 0.0
    direction: flat                # abs(r) < 0.20
    source: "Srimani 2022 CPT:PSP"

  - metric: Cox_HR_deepresp_vs_shallow
    value: 0.26
    range: [0.20, 0.45]
    source: "Srimani 2022 CPT:PSP Table 4"
```

---

## Validation Target Tolerance Conventions

| Endpoint type | Default tolerance | Rationale |
|--------------|-----------------|-----------|
| Primary survival (PFS/OS median) | ±5% | Primary endpoint, tightest |
| Secondary efficacy (ORR, VGPR+) | ±15% | Allows realistic noise |
| Safety (Grade ≥3 rates) | ±20% | Higher variability in real data |
| PK NCA (Cmax, AUCinf) | ±15% | Published analytic variability |
| PK IIV CV% | ±10 percentage points | |
| Demographic % | ±10% | |
| Cross-correlation r | ±0.10 absolute | |
| Cox HR range | target range | Use published 95% CI if available |

---

## Naming Conventions

- Branch name = `trial.name` (e.g., `POLLUX`, `CLL14`, `FLAURA`)
- Data directory = `trials/<trial.name>/`
- Config file = `trials/<trial.name>/config.yaml`
- Generated data = `trials/<trial.name>/data/` (SDTM) and `trials/<trial.name>/adam/`
- Validation report = `trials/<trial.name>/outputs/VALIDATION_REPORT.md`
- Agent sign-off = `trials/<trial.name>/outputs/agent_signoff.md`
