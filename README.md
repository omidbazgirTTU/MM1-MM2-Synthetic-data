# Synthetic TOURMALINE MM1 / MM2 Trial Data

Fully synthetic clinical trial datasets modeled after the TOURMALINE-MM2 (newly diagnosed
MM) and TOURMALINE-MM1 (relapsed/refractory MM) trials. Generated to support development
of a pharmacology-informed causal deep learning framework (PK → PD → Endpoints) for
personalized dosing optimization in Multiple Myeloma.

> Hussain Z, De Brouwer E, et al. "Joint AI-driven event prediction and longitudinal
> modeling in newly diagnosed and relapsed multiple myeloma."
> *npj Digital Medicine* 7, 200 (2024). https://doi.org/10.1038/s41746-024-01189-3

**All values are synthetically generated and do not correspond to any real patients.**

---

## Validation Status

**68 / 68 criteria PASS** against published TOURMALINE trial summary statistics and
mechanistic pharmacometric targets.

| Category | Criteria | Status |
|----------|----------|--------|
| Survival (PFS/OS medians) | 8 | All PASS (±0.3%) |
| Enrollment (N per arm) | 4 | All PASS (exact) |
| Demographics (age, sex) | 4 | All PASS (±2.4%) |
| ISS staging | 6 | All PASS (±8.5%) |
| Cytogenetics + renal function | 8 | All PASS (±11.7%) |
| Efficacy (ORR, VGPR+, CR+) | 12 | All PASS (±11.5%) |
| PK — Ixazomib Cmax (NCA, Cycle 1) | 2 | All PASS (+8–10%) |
| Safety — Grade 3 PLT | 4 | All PASS (±11%) |
| **Baseline covariate cross-correlations** | **8** | **All PASS** |
| **PK NCA cross-correlations** | **4** | **All PASS** |
| **AUC → PLT nadir depth** | **2** | **All PASS** |
| **M-protein Cycle 6 → PFS (Cox HR)** | **2** | **All PASS (HR 0.31–0.33)** |
| **Exposure-efficacy flatness** | **4** | **All PASS** |

Full results: `outputs/tables/validation_summary.csv`

---

## Mechanistic Model

The synthetic data embeds published pharmacokinetic/pharmacodynamic relationships from the
TOURMALINE-MM1 popPK/PD analysis (Srimani 2022, *CPT:PSP*):

| Mechanism | Implementation | Reference |
|-----------|---------------|-----------|
| 3-compartment Ixazomib PK | Gupta 2017 params, BSA on V4, CYP3A4 DDI | Gupta 2017 |
| Linear AUC → PLT dip | `dip_i = dip_pop × AUC_i/AUC_pop` (cap 2.5×) | Srimani 2022 |
| Flat AUC → efficacy | \|r(AUC, M-prot C6)\| < 0.20, ORR p > 0.05 | Srimani 2022 |
| M-protein → PFS link | Gaussian copula (ρ = −0.80), IRd arm | Srimani 2022 |
| 7×7 MVN baseline covariates | Age↔CrCL r=−0.45, M-prot↔HGB r=−0.30 | Published |
| OMEGA Cholesky PK etas | ρ(CL,V2)=0.30, ρ(CL,V4)=0.20, ρ(V2,V4)=0.25 | Gupta 2017 |

See [`MECHANISTIC_MODEL_AND_CROSS_CORRELATIONS.md`](MECHANISTIC_MODEL_AND_CROSS_CORRELATIONS.md)
for the complete model specification.

---

## File Structure

```
Takeda-data/
├── MM2/                          TOURMALINE-MM2 — Newly Diagnosed MM (NDMM), N=705
│   ├── sdtm_dm.csv               Demographics, anthropometrics, cytogenetics
│   ├── sdtm_ex.csv               Dosing (IRd vs Rd, dose modifications, adherence)
│   ├── sdtm_lb.csv               34 longitudinal lab tests per cycle
│   ├── sdtm_ae.csv               12 AEs with CTCAE grading
│   ├── sdtm_ds.csv               Disposition / survival events
│   ├── sdtm_cm.csv               Concomitant medications (CYP3A4, anticoagulants, G-CSF)
│   ├── sdtm_pc.csv               PK concentrations (sparse + dense, Cycles 1 & 3)
│   ├── sdtm_pp.csv               NCA-derived PK parameters (Cmax, AUCinf, t½, CL/F)
│   ├── sdtm_rs.csv               Response assessments (IMWG criteria)
│   ├── adam_adsl.csv             Subject-level: covariates, survival, cytogenetics, R-ISS
│   ├── adam_adtte.csv            Time-to-event (PFS + OS)
│   ├── adam_adlb.csv             Longitudinal lab ADaM (CHG, PCHG, ABLFL)
│   ├── adam_adpc.csv             PK analysis dataset with covariates
│   └── adam_adrs.csv             Response analysis dataset
│
├── MM1/                          TOURMALINE-MM1 — Relapsed/Refractory MM (RRMM), N=722
│   └── (same structure as MM2)
│
├── scripts/
│   ├── generate_v2.py            Main SDTM/ADaM generator (all non-PK domains)
│   │                             Tracks A1, A3–A6: MVN covariates, AUC→response,
│   │                             AUC→PLT, HGB adjustment, Gaussian copula PFS
│   ├── generate_pk_v2.py         PK concentration + NCA generator (Track A2: OMEGA Cholesky)
│   ├── validate_data.py          68-criterion validation script (Tracks A1–A7)
│   ├── pk_vpc_gof.py             VPC and GOF figure generator
│   └── plot_individual_patients.py  Per-patient longitudinal profile figures
│
├── outputs/
│   ├── tables/
│   │   └── validation_summary.csv   Machine-readable 68-row pass/fail table
│   └── validation/               Validation figures (KM, M-protein, PLT, PK VPC,
│                                 covariate distributions, cross-correlation scatters)
│
├── .claude/skills/               Reusable pharmacometric knowledge
│   ├── Cross_Correlations_Synthetic_Data_Guide.md
│   └── Published_Mechanistic_Correlations_TOURMALINE.md
│
├── MECHANISTIC_MODEL_AND_CROSS_CORRELATIONS.md   Full PK/PD model + correlation spec
├── TOURMALINE_Synthetic_Data_Specifications.md   Original data generation spec
├── TOURMALINE-MM1_Schedule_of_Assessments.md
└── TOURMALINE-MM2_Schedule_of_Assessments.md
```

*MM1/ and MM2/ directories are excluded from git (reproducible by running the scripts).*

---

## Quick Start

```bash
# 1. Set up environment
python3 -m venv .venv && source .venv/bin/activate
pip install numpy pandas scipy lifelines matplotlib

# 2. Generate data (both studies, ~3–5 minutes)
python3 scripts/generate_v2.py
python3 scripts/generate_pk_v2.py

# 3. Validate (should reach 68/68 PASS)
python3 scripts/validate_data.py
```

Seeds are fixed (`MM2=42, MM1=43, SURV_RNG=77`) — all outputs are fully reproducible.

---

## SDTM Domains

### DM — Demographics
- N = 705 (MM2, NDMM) / 722 (MM1, RRMM)
- Age, sex, race, ECOG, Ig type, light chain type, Durie-Salmon stage
- Anthropometrics: weight (kg), height (cm), BMI (kg/m²), BSA (Mosteller, m²)
- ISS staging (I/II/III) per published proportions
- **Baseline covariates drawn from 7×7 MVN** (Age, CrCL, Weight, BSA, M-protein, PLT, HGB)
  with published physiological correlation structure

### EX — Exposure / Dosing
- **IRd arm**: Ixazomib 4 mg (Days 1, 8, 15), Lenalidomide 25 mg (Days 1–21), Dexamethasone 40 mg (Days 1, 8, 15, 22) — 28-day cycles
- **Rd arm**: Placebo (Days 1, 8, 15), Lenalidomide 25 mg, Dexamethasone 40 mg
- Dose modifications (~15%), adherence flags, dose holiday durations, modification reasons
- Up to 26 cycles (~2 years follow-up)

### LB — Laboratory Results (34 tests, longitudinal)
- **Hematology** (7): HGB, HCT, NEUT, PLT, LYMPH, MONO, WBC
- **Chemistry** (19): ALBUMIN, ALP, ALT, AST, BILI, B2MG, BUN, CA, CA_CORR, CL, CO2, CREAT, GFR, GLOB, GLUC, LDH, MG, NA, PHOS, PROT, URATE
- **Disease biomarkers** (8): IGA, IGG, IGM, KAPPA_LAMBDA, SPEP_GAMMA, SPEP_KAPPA, SPEP_LAMBDA, SPEP_MPROT, UPEP_MPROT, BMPC
- PLT nadir at Days 11–15 per cycle (mechanistic; Srimani 2022)
- AUC-proportional PLT dip per patient (linear model; not Emax)
- ~5–10% missing at random per visit

### AE — Adverse Events
12 AEs (CTCAE grading), rates matched to published TOURMALINE incidence:
- Hematologic: Neutropenia, Thrombocytopenia
- Non-hematologic: Diarrhea (prior IMiD covariate), Nausea, Vomiting, Peripheral Neuropathy, Rash (race covariate), Acute Renal Failure, Cardiac Arrhythmias, Heart Failure, Hypotension, Liver Impairment

### DS — Disposition / Survival
- PFS and OS events with Weibull distributions (calibrated KM medians)
- IRd vs Rd PFS HR ≈ 0.74
- M-protein Cycle 6 response linked to PFS via Gaussian copula (ρ = −0.80, IRd arm only)
- Cox PH HR for ≥75% responders: 0.31–0.33 (target 0.20–0.45)

### CM — Concomitant Medications
- CYP3A4 strong inhibitors (~8%): clarithromycin, itraconazole, ketoconazole
- CYP3A4 strong inducers (~3%): rifampin, carbamazepine
- Anticoagulants (~85%): aspirin / LMWH (lenalidomide DVT prophylaxis)
- Supportive care: G-CSF (~40%), EPO (~20%), transfusions

### PC / PP — PK Concentrations and NCA Parameters
- Ixazomib: 3-compartment oral (Gupta 2017), full IIV via OMEGA Cholesky
- Lenalidomide: 1-compartment oral (Chen 2012), CrCL covariate
- Dexamethasone: 1-compartment oral
- Sparse PK substudy: ~150 subjects/study, sampled at Cycles 1 & 3
- NCA in PP: Cmax, AUCinf, t½, CL/F (Cycle 1 scheduled timepoints only)
- Individual `IXAZ_CL_I` in `adam_adsl.csv` links PK to PD responses

---

## ADaM Datasets

### ADSL — Subject Level
Baseline covariates, survival endpoints, treatment flags, risk stratification:
- PFS/OS duration (months) + censoring flags (CNSR=0 event, CNSR=1 censored)
- `IXAZ_CL_I`: per-patient Ixazomib CL (L/h) — shared with PK generator
- `BASE_SPEP_MPROT_MVN`: M-protein from MVN draw (g/dL, correlated with HGB)
- `BASE_PLT_MVN`, `BASE_HGB_MVN`: PLT/HGB from MVN draw
- Cytogenetics: del(17p), t(4;14), t(14;16), t(14;20), gain(1q21), del(1p32), AMP1Q
- R-ISS stage (I/II/III): derived from ISS + del(17p)/t(4;14) + LDH
- CrCL (Cockcroft-Gault), renal group (RENGRP)
- CYP3A4 inhibitor/inducer flags, DDI CL multiplier

### ADTTE — Time-to-Event
Long format (2 records per subject: PFS + OS). M-protein → PFS copula applied in IRd arm.

### ADLB — Longitudinal Lab ADaM
CHG, PCHG, ABLFL, ANL01FL per lab test per visit.

### ADPC — PK Analysis Dataset
PC records merged with ADSL covariates (BSA, CrCL, CYP3A4 flags, Cmax, AUCinf).

### ADRS — Response Analysis Dataset
IMWG best-response assignments per patient.

---

## Simulation Design

### Seeds and Reproducibility
| Scope | Seed |
|-------|------|
| MM2 generation | 42 |
| MM1 generation | 43 |
| Survival (SURV_RNG) | 77 |

### PK Models
| Drug | Model | Key parameters |
|------|-------|----------------|
| Ixazomib | 3-cmt oral (Gupta 2017) | CL=1.86 L/h, Vss=543 L, t½=228h, F=58%, Ka=0.5 h⁻¹ |
| Lenalidomide | 1-cmt oral (Chen 2012) | CL/F=8.94 L/h, CrCL covariate (power 0.60) |
| Dexamethasone | 1-cmt oral | CL/F=16 L/h |

### Response Calibration
IMWG response rates calibrated using bimodal phenotype mixture model:
- IRd: ORR ~82% (MM2), ~78% (MM1)
- Rd: ORR ~75% (MM2), ~72% (MM1)
- AUC→response link is weak (flat E-R within therapeutic range per Srimani 2022)

---

## Usage Notes

- **Time unit**: ADY=1 = Day 1 of first dose; VISITNUM = cycle number
- **Censoring**: CNSR=0 → event observed; CNSR=1 → censored
- **Response**: per-patient best PCHG from SPEP_MPROT (ORR ≤−50%, VGPR+ ≤−90%, CR+ ≤−99%)
- **PK Cmax**: use `sdtm_pp.csv` (NCA, Cycle 1), not raw max from `sdtm_pc.csv`
- **M-protein cross-correlation**: use `BASE_SPEP_MPROT_MVN` for HGB/PLT correlation analyses
- **Missing data**: use LOCF or model-based imputation; ABLFL flags baseline records

---

## Citation

If using this synthetic dataset, please contact omidbazgir00@gmail.com or ranadip.pal@ttu.edu
