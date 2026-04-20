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

**48 / 48 criteria PASS** against published TOURMALINE trial summary statistics.

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

Full results: `outputs/VALIDATION_REPORT.md` · Machine-readable: `outputs/tables/validation_summary.csv`

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
│   ├── generate_pk_v2.py         PK concentration + NCA generator (3-cmt Ixazomib)
│   ├── validate_data.py          48-criterion validation script
│   ├── pk_vpc_gof.py             VPC and GOF figure generator
│   └── plot_individual_patients.py  Per-patient longitudinal profile figures
│
├── outputs/
│   ├── VALIDATION_REPORT.md      Comprehensive validation report (48 criteria + PK VPC/GOF)
│   ├── tables/
│   │   └── validation_summary.csv   Machine-readable 48-row pass/fail table
│   └── figures/
│       ├── MM2/                  705 individual patient figures
│       ├── MM1/                  722 individual patient figures
│       ├── pk_vpc_MM2.png        VPC — Ixazomib / Lenalidomide / Dexamethasone (MM2)
│       ├── pk_vpc_MM1.png        VPC — MM1
│       ├── pk_gof_MM2.png        GOF panels — Cmax, AUCinf, t½, covariates, CWRES (MM2)
│       └── pk_gof_MM1.png        GOF panels — MM1
│
└── .claude/skills/               Reusable pharmacometric knowledge for this project
```

---

## SDTM Domains

### DM — Demographics
- N = 705 (MM2, NDMM) / 722 (MM1, RRMM)
- Age, sex, race, ECOG, Ig type, light chain type, Durie-Salmon stage
- Anthropometrics: weight (kg), height (cm), BMI (kg/m²), BSA (Mosteller, m²)
- ISS staging (I/II/III) per published proportions

### EX — Exposure / Dosing
- **IRd arm**: Ixazomib 4 mg (Days 1, 8, 15), Lenalidomide 25 mg (Days 1–21), Dexamethasone 40 mg (Days 1, 8, 15, 22) — 28-day cycles
- **Rd arm**: Placebo (Days 1, 8, 15), Lenalidomide 25 mg, Dexamethasone 40 mg
- Dose modifications (~15%), adherence flags, dose holiday durations, modification reasons (AE / protocol / other)
- Up to 26 cycles (~2 years follow-up)

### LB — Laboratory Results (34 tests, longitudinal)
- **Hematology** (7): HGB, HCT, NEUT, PLT, LYMPH, MONO, WBC
- **Chemistry** (19): ALBUMIN, ALP, ALT, AST, BILI, B2MG, BUN, CA, CA_CORR, CL, CO2, CREAT, GFR, GLOB, GLUC, LDH, MG, NA, PHOS, PROT, URATE
- **Disease biomarkers** (8): IGA, IGG, IGM, KAPPA_LAMBDA, SPEP_GAMMA, SPEP_KAPPA, SPEP_LAMBDA, SPEP_MPROT, UPEP_MPROT, BMPC
- Physiologically plausible trajectories: M-protein ↓ with treatment (deeper nadir in IRd), HGB recovery, renal function linked to disease burden
- ~5–10% missing at random per visit (3% for primary efficacy endpoints)

### AE — Adverse Events
12 AEs (CTCAE grading), rates matched to published TOURMALINE incidence:
- Hematologic: Neutropenia, Thrombocytopenia (Grade 3 PLT: IRd 25%/32%, Rd 14%/16% for NDMM/RRMM)
- Non-hematologic: Diarrhea, Nausea, Vomiting, Peripheral Neuropathy, Rash, Acute Renal Failure, Cardiac Arrhythmias, Heart Failure, Hypotension, Liver Impairment

### DS — Disposition / Survival
- PFS and OS events with censoring flags
- Survival calibrated to published KM medians (PFS: MM2 ≈18 mo, MM1 ≈12 mo; OS: MM2 ≈40 mo, MM1 ≈30 mo)
- IRd vs Rd PFS HR ≈ 0.74

### CM — Concomitant Medications
- CYP3A4 strong inhibitors (~8% of IRd patients): clarithromycin, itraconazole, ketoconazole
- CYP3A4 strong inducers (~3%): rifampin, carbamazepine
- Anticoagulants (~85%): aspirin / LMWH (lenalidomide DVT prophylaxis)
- Supportive care: G-CSF (~40%), EPO (~20%), transfusions

### PC / PP — PK Concentrations and NCA Parameters
Three drugs: Ixazomib (3-compartment, Gupta 2017), Lenalidomide (1-compartment, Chen 2012), Dexamethasone (1-compartment).
- Sparse PK substudy: ~150 subjects per study, sampled at Cycles 1 & 3 (dense) + trough cycles
- Full IIV on all PK parameters (log-normal); proportional + additive residual error; BLQ flags
- Covariates: BSA on Ixazomib V4 (power 0.70), CrCL on Lenalidomide CL (power 0.60), CYP3A4 DDI on Ixazomib CL (×0.55 inhibitor / ×2.0 inducer)
- NCA parameters in PP: Cmax, AUCinf, t½, CL/F (Cycle 1 scheduled timepoints)

---

## ADaM Datasets

### ADSL — Subject Level
Baseline covariates, survival endpoints, treatment flags, risk stratification:
- PFS/OS duration (months) + censoring flags (CNSR=0 event, CNSR=1 censored)
- Treatment flags: ITTFL, SAFFL, PPSFL
- Cytogenetics: del(17p), t(4;14), t(14;16), t(14;20), gain(1q21), del(1p32), AMP1Q
- R-ISS stage (I/II/III): derived from ISS + del(17p)/t(4;14) + LDH
- CrCL (Cockcroft-Gault), renal group (RENGRP)
- CYP3A4 inhibitor/inducer flags, DDI CL multiplier

### ADTTE — Time-to-Event
Long format (2 records per subject: PFS + OS), compatible with standard KM/Cox analyses.

### ADLB — Longitudinal Lab ADaM
CHG, PCHG, ABLFL, ANL01FL per lab test per visit. Used for IMWG response computation
(ORR, VGPR+, CR+ from SPEP_MPROT PCHG).

### ADPC — PK Analysis Dataset
PC records merged with ADSL covariates (BSA, CrCL, CYP3A4 flags, CMAX_PP, AUCINF_PP).

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

Per-study reseeding prevents cross-study RNG contamination: changing MM2 generation code
does not alter MM1 results.

### PK Models
| Drug | Model | Key parameters |
|------|-------|----------------|
| Ixazomib | 3-compartment oral (Gupta 2017) | CL=1.86 L/h, Vss=543 L, t½=228h, F=58%, Ka=0.5 h⁻¹ |
| Lenalidomide | 1-compartment oral (Chen 2012) | CL/F=8.94 L/h, CrCL covariate |
| Dexamethasone | 1-compartment oral | CL/F=16 L/h |

### Response Calibration
IMWG response rates (ORR, VGPR+, CR+) calibrated using a bimodal phenotype mixture
model with three mechanisms: guaranteed nadir override at Cycle 9 for all responders,
disease-biomarker miss_rate=3%, and CR tier starting at resp_rate=0.993 to avoid
boundary noise misclassification.

---

## Usage Notes

- **Time unit**: ADY=1 = Day 1 of first dose; VISITNUM = cycle number
- **Censoring**: CNSR=0 → event observed; CNSR=1 → censored
- **Response computation**: per-patient best PCHG from SPEP_MPROT in ADLB
  (ORR ≤ −50%, VGPR+ ≤ −90%, CR+ ≤ −99%)
- **PK Cmax**: use `sdtm_pp.csv` (NCA, Cycle 1) not raw max from `sdtm_pc.csv`
  (which contains multi-dose accumulated concentrations)
- **Missing data**: use LOCF or model-based imputation; ABLFL flags baseline records

---

## Citation

If using this synthetic dataset, please reach out to omidbazgir00@gmail.com or ranadip.pal@ttu.edu