"""
Generate a data dictionary and README for the synthetic TOURMALINE datasets
"""
import pandas as pd
import os

out_dir = "/mnt/user-data/outputs/tourmaline_synthetic"

# ──────────────────────────────────────────────────────────────────────────────
# DATA DICTIONARY
# ──────────────────────────────────────────────────────────────────────────────

dict_rows = {
    "SDTM_DM": [
        ("STUDYID","Study Identifier","CHAR","Study name"),
        ("DOMAIN","Domain Abbreviation","CHAR","DM"),
        ("USUBJID","Unique Subject Identifier","CHAR","Study-Subject ID"),
        ("SUBJID","Subject Identifier","CHAR","4-digit zero-padded integer"),
        ("RFSTDTC","Subject Reference Start Date","CHAR","ISO 8601 date"),
        ("SITEID","Study Site Identifier","CHAR","Site code"),
        ("AGE","Age at Enrollment","NUM","Years"),
        ("AGEU","Age Units","CHAR","YEARS"),
        ("SEX","Sex","CHAR","M/F"),
        ("RACE","Race","CHAR","CTCAE race terms"),
        ("ETHNIC","Ethnicity","CHAR","Hispanic/Not Hispanic"),
        ("COUNTRY","Country","CHAR","3-letter ISO country"),
        ("ARMCD","Planned Arm Code","CHAR","IRd/Rd"),
        ("ARM","Description of Planned Arm","CHAR","Full treatment description"),
        ("IGTYPE","Immunoglobulin Type","CHAR","IgG/IgA/IgD/IgE/IgM/Biclonal/No Heavy Chain"),
        ("LCTYPE","Light Chain Type","CHAR","KAPPA/LAMBDA/BICLONAL"),
        ("ISSSTAGE","ISS Stage at Study Entry","NUM","1/2/3"),
        ("ECOG","ECOG Performance Status","NUM","0/1/2"),
        ("DXMONTHS","Time Since Diagnosis","NUM","Months"),
        ("POPULATION","Study Population","CHAR","NDMM/RRMM"),
    ],
    "SDTM_EX": [
        ("STUDYID","Study Identifier","CHAR",""),
        ("USUBJID","Unique Subject Identifier","CHAR",""),
        ("EXSEQ","Sequence Number","NUM","Auto-incremented per subject"),
        ("EXTRT","Name of Treatment","CHAR","IXAZOMIB/LENALIDOMIDE/DEXAMETHASONE/PLACEBO"),
        ("EXDOSE","Dose per Administration","NUM","mg"),
        ("EXDOSU","Dose Units","CHAR","mg"),
        ("EXDOSFRM","Dose Form","CHAR","TABLET"),
        ("EXDOSFRQ","Dosing Frequency","CHAR","Schedule description"),
        ("EXROUTE","Route of Administration","CHAR","ORAL"),
        ("EXSTDTC","Start Date of Treatment","CHAR","ISO 8601"),
        ("EXENDTC","End Date of Treatment","CHAR","ISO 8601"),
        ("VISITNUM","Visit Number","NUM","Cycle number"),
        ("VISIT","Visit Name","CHAR","CYCLE N"),
        ("EXDOSMOD","Dose Modification","CHAR","DOSE REDUCTION/blank"),
    ],
    "SDTM_LB": [
        ("STUDYID","Study Identifier","CHAR",""),
        ("USUBJID","Unique Subject Identifier","CHAR",""),
        ("LBSEQ","Sequence Number","NUM",""),
        ("LBCAT","Category","CHAR","HEMATOLOGY/CHEMISTRY/SERUM IMMUNOGLOBULINS"),
        ("LBTESTCD","Lab Test Short Name","CHAR","See lab codebook"),
        ("LBTEST","Lab Test Long Name","CHAR","Full descriptive name"),
        ("LBORRES","Result in Original Units","CHAR","Raw value"),
        ("LBORRESU","Units","CHAR","Original units"),
        ("LBSTRESN","Numeric Standardized Result","NUM",""),
        ("LBSTRESU","Standard Units","CHAR",""),
        ("LBNRLO","Normal Range Lower Limit","NUM",""),
        ("LBNRHI","Normal Range Upper Limit","NUM",""),
        ("LBNRIND","Normal Range Indicator","CHAR","LOW/NORMAL/HIGH"),
        ("VISITNUM","Visit Number","NUM","Cycle number (1=baseline)"),
        ("LBDTC","Date of Lab Test","CHAR","ISO 8601"),
        ("LBDY","Study Day","NUM","Day 1 = study start"),
        ("EPOCH","Epoch","CHAR","BASELINE/TREATMENT"),
    ],
    "SDTM_AE": [
        ("STUDYID","Study Identifier","CHAR",""),
        ("USUBJID","Unique Subject Identifier","CHAR",""),
        ("AESEQ","Sequence Number","NUM",""),
        ("AETERM","Reported Term","CHAR","Verbatim AE term"),
        ("AEDECOD","Dictionary-Derived Term","CHAR","MedDRA preferred term"),
        ("AEBODSYS","Body System","CHAR","MedDRA SOC"),
        ("AESEV","Severity","CHAR","MILD/MODERATE/SEVERE/LIFE THREATENING/FATAL"),
        ("AETOXGR","CTCAE Toxicity Grade","CHAR","1-5"),
        ("AESER","Serious Event","CHAR","Y/N (≥Grade 3)"),
        ("AEREL","Causality","CHAR","RELATED/NOT RELATED/POSSIBLY RELATED"),
        ("AEACN","Action Taken","CHAR","DOSE REDUCTION/WITHDRAWN/NOT CHANGED/INTERRUPTED"),
        ("AESTDTC","AE Start Date","CHAR","ISO 8601"),
        ("AEENDTC","AE End Date","CHAR","ISO 8601"),
        ("AEDY","AE Start Day","NUM","Study day"),
        ("AEOUT","Outcome","CHAR","Resolved/Recovering/Not resolved"),
        ("EPOCH","Epoch","CHAR","TREATMENT"),
    ],
    "SDTM_DS": [
        ("STUDYID","Study Identifier","CHAR",""),
        ("USUBJID","Unique Subject Identifier","CHAR",""),
        ("DSSEQ","Sequence Number","NUM","1=PFS, 2=OS"),
        ("DSDECOD","Standardized Disposition Term","CHAR","DISEASE PROGRESSION/DEATH/CENSORED"),
        ("DSSCAT","Disposition Subcategory","CHAR","PROGRESSION FREE SURVIVAL/OVERALL SURVIVAL"),
        ("DSSTDTC","Date of Disposition","CHAR","ISO 8601"),
        ("EVENTFL","Event Flag","CHAR","Y=event occurred, N=censored"),
        ("CNSR","Censoring Indicator","NUM","0=event, 1=censored"),
        ("AVAL","Analysis Value","NUM","Time in months"),
        ("PARAMCD","Parameter Code","CHAR","PFS/OS"),
        ("PARAM","Parameter","CHAR","Full description"),
    ],
    "ADaM_ADSL": [
        ("STUDYID","Study Identifier","CHAR",""),
        ("USUBJID","Unique Subject Identifier","CHAR",""),
        ("TRT01P","Planned Treatment for Period 1","CHAR","Full arm description"),
        ("TRT01PN","Planned Treatment Period 1 (N)","NUM","1=IRd, 2=Rd"),
        ("ITTFL","Intent-to-Treat Population Flag","CHAR","Y/N"),
        ("SAFFL","Safety Population Flag","CHAR","Y/N"),
        ("PFSDUR","PFS Duration","NUM","Months"),
        ("PFSCNSR","PFS Censoring Indicator","NUM","0=event, 1=censored"),
        ("OSDUR","OS Duration","NUM","Months"),
        ("OSCNSR","OS Censoring Indicator","NUM","0=event, 1=censored"),
        ("RISKGR","Risk Group","CHAR","LOW/INTERMEDIATE/HIGH (ISS-based)"),
        ("CYTOGR","Cytogenetics Group","CHAR","HIGH RISK/STANDARD RISK"),
        ("DEL17P","del(17p) Abnormality","CHAR","Y/N"),
        ("T414","t(4;14) Abnormality","CHAR","Y/N"),
        ("AMP1Q","1q Amplification","CHAR","Y/N"),
        ("POLYCC","CC Polymorphism","NUM","0/1"),
        ("POLYCG","CG Polymorphism","NUM","0/1"),
        ("POLYGG","GG Polymorphism","NUM","0/1"),
        ("DSSTAGE","Durie-Salmon Stage","CHAR","IA/IIA/IIB/IIIA/IIIB"),
        ("BASE_CREAT","Baseline Creatinine","NUM","umol/L"),
        ("BASE_CREACL","Baseline Creatinine Clearance","NUM","mL/min"),
        ("BASE_ALB","Baseline Albumin","NUM","g/L"),
        ("BASE_B2MG","Baseline Beta-2 Microglobulin","CHAR","Categorical: <3.5/3.5-5.5/>5.5"),
        ("BASE_FLC","Baseline Involved Free Light Chain","NUM","mg/L"),
        ("BASE_BMPC","Bone Marrow Plasma Cells","NUM","%"),
        ("MEAS_DISFL","Measurable Disease Flag","CHAR","Measurement type used"),
    ],
    "ADaM_ADTTE": [
        ("STUDYID","Study Identifier","CHAR",""),
        ("USUBJID","Unique Subject Identifier","CHAR",""),
        ("PARAMCD","Parameter Code","CHAR","PFS/OS"),
        ("PARAM","Parameter","CHAR","Full description"),
        ("AVAL","Analysis Value","NUM","Time in months (primary analysis unit)"),
        ("AVALU","Analysis Value Unit","CHAR","MONTHS"),
        ("CNSR","Censoring Indicator","NUM","0=event, 1=censored"),
        ("EVNTDESC","Event Description","CHAR","DISEASE PROGRESSION/DEATH/CENSORED"),
        ("STARTDT","Start Date","CHAR","ISO 8601 reference start"),
        ("STRATFL1","Stratification Factor 1","CHAR","ISS Stage"),
        ("STRATFL2","Stratification Factor 2","CHAR","Cytogenetics group"),
        ("STRATFL3","Stratification Factor 3","CHAR","Ig type"),
        ("ARMCD","Treatment Arm Code","CHAR","IRd/Rd"),
        ("TRT01PN","Planned Treatment (Numeric)","NUM","1/2"),
    ],
    "ADaM_ADLB": [
        ("STUDYID","Study Identifier","CHAR",""),
        ("USUBJID","Unique Subject Identifier","CHAR",""),
        ("PARAMCD","Parameter Code","CHAR","Lab test code"),
        ("PARAM","Parameter","CHAR","Lab test description"),
        ("AVAL","Analysis Value","NUM","Lab result"),
        ("AVALU","Analysis Value Unit","CHAR","Test-specific unit"),
        ("BASE","Baseline Value","NUM","Value at EPOCH=BASELINE"),
        ("CHG","Change from Baseline","NUM","AVAL - BASE"),
        ("PCHG","Percent Change from Baseline","NUM","(CHG/BASE)*100"),
        ("ABLFL","Baseline Record Flag","CHAR","Y=baseline record"),
        ("ANL01FL","Analysis Flag 01","CHAR","Y=included in primary analysis"),
        ("LBCAT","Lab Category","CHAR","HEMATOLOGY/CHEMISTRY/SERUM IMMUNOGLOBULINS"),
        ("LBNRIND","Normal Range Indicator","CHAR","LOW/NORMAL/HIGH"),
        ("VISITNUM","Visit Number","NUM","Cycle number"),
        ("EPOCH","Epoch","CHAR","BASELINE/TREATMENT"),
        ("ARMCD","Treatment Arm Code","CHAR","IRd/Rd"),
        ("IGTYPE","Immunoglobulin Type","CHAR","Myeloma subtype"),
    ],
}

# Write Excel data dictionary
with pd.ExcelWriter(f"{out_dir}/data_dictionary.xlsx", engine="openpyxl") as writer:
    for sheet, rows in dict_rows.items():
        df = pd.DataFrame(rows, columns=["Variable","Label","Type","Notes/Values"])
        df.to_excel(writer, sheet_name=sheet[:31], index=False)

print("✓ Data dictionary written.")

# ──────────────────────────────────────────────────────────────────────────────
# README
# ──────────────────────────────────────────────────────────────────────────────

readme = """# Synthetic TOURMALINE MM1 / MM2 Trial Data
## Overview
This dataset contains **fully synthetic** clinical trial data modeled after the
TOURMALINE-MM1 (relapsed/refractory MM) and TOURMALINE-MM2 (newly diagnosed MM)
trials, as described in:

> Hussain Z, De Brouwer E, et al. "Joint AI-driven event prediction and longitudinal
> modeling in newly diagnosed and relapsed multiple myeloma."
> *npj Digital Medicine* 7, 200 (2024). https://doi.org/10.1038/s41746-024-01189-3

**All values are synthetically generated and do not correspond to any real patients.**

---
## File Structure
```
tourmaline_synthetic/
├── MM2/                         ← Newly Diagnosed MM (NDMM), N=703
│   ├── sdtm_dm.csv              SDTM Demographics
│   ├── sdtm_ex.csv              SDTM Exposure (Dosing)
│   ├── sdtm_lb.csv              SDTM Laboratory Results
│   ├── sdtm_ae.csv              SDTM Adverse Events
│   ├── sdtm_ds.csv              SDTM Disposition
│   ├── adam_adsl.csv            ADaM Subject-Level Analysis Dataset
│   ├── adam_adtte.csv           ADaM Time-to-Event Analysis Dataset
│   └── adam_adlb.csv            ADaM Longitudinal Lab Analysis Dataset
├── MM1/                         ← Relapsed/Refractory MM (RRMM), N=720
│   └── (same structure)
├── data_dictionary.xlsx         Variable-level data dictionary
└── dataset_summary.csv          Row counts per dataset
```

---
## SDTM Domains

### DM — Demographics
- N = 703 (MM2) / 720 (MM1), matching published Table 1
- Demographics: age (median 73/66), sex (~50%/43% female), race
- Ig type: IgG ~57%/54%, IgA ~20%/17%, etc.
- ISS Stage at entry: I/II/III per published proportions
- Randomization: ~50/50 IRd vs Rd (actual and planned)
- Additional: ECOG, Durie-Salmon stage, Ig type, light chain type

### EX — Exposure / Dosing
- **IRd arm**: Ixazomib 4mg (Days 1,8,15), Lenalidomide 25mg (Days 1–21),
  Dexamethasone 40mg (Days 1,8,15,22) — per 28-day cycle
- **Rd arm**: Placebo (Days 1,8,15), Lenalidomide 25mg, Dexamethasone 40mg
- Dose modifications (~15%) encoded as EXDOSMOD
- Up to 26 cycles (≈2 years of follow-up)

### LB — Laboratory Results (Longitudinal Biomarkers)
34 lab tests across 3 categories, assessed each 28-day cycle:
- **Hematology** (8): Hemoglobin, Hematocrit, Neutrophils, Platelets, Lymphocytes,
  Monocytes, Leukocytes
- **Chemistry** (18): Albumin, Creatinine, GFR, Calcium, Electrolytes, LFTs,
  BUN, LDH, etc.
- **Serum Immunoglobulins** (8): IgG, IgA, IgM, SPEP M-protein, Kappa/Lambda FLC,
  UPEP M-protein — **key disease biomarkers for SCOPE model**

Lab trajectories are physiologically simulated:
- Disease markers (M-protein, IgA/G) ↓ with treatment (deeper nadir in IRd)
- Hemoglobin low at baseline, gradual recovery
- Creatinine/GFR may deteriorate with disease

### AE — Adverse Events
12 AEs tracked (from paper Figure 2b), with CTCAE grading:
- Hematologic (≥Grade 3): Neutropenia, Thrombocytopenia
- Non-hematologic (≥Grade 2): Acute Renal Failure, Cardiac Arrhythmias,
  Diarrhea, Heart Failure, Hypotension, Liver Impairment, Nausea,
  Peripheral Neuropathies, Rash, Vomiting
- Includes: AETOXGR, AESER, AEREL, AEACN, onset/end dates, outcome

### DS — Disposition / Survival
- PFS events (disease progression) and OS events (death) with censoring
- Simulated using Weibull survival: median PFS ≈18mo (MM2), 12mo (MM1)
- Includes CNSR flag for use in Cox/KM analyses

---
## ADaM Datasets

### ADSL — Subject Level Analysis Dataset
Extends DM with derived analysis variables:
- PFS/OS duration (months) + censoring flags
- Treatment flags: ITTFL, SAFFL, PPSFL
- Risk stratification: ISS-based RISKGR, CYTOGR (del17p, t(4;14), amp1q)
- Polymorphisms: CC/CG/GG binary flags (used in paper's subgroup discovery)
- Baseline labs: creatinine, creatinine clearance, albumin, β2-microglobulin, FLC
- Durie-Salmon stage, measurable disease flag, bone lesions, extramedullary disease

### ADTTE — Time-to-Event Analysis Dataset
Long format (2 records per subject: PFS + OS):
- AVAL in months (primary analysis unit per paper)
- Stratification factors: ISS stage, cytogenetics, Ig type
- Compatible with standard PROC LIFETEST / survfit() analyses

### ADLB — Longitudinal Lab Analysis Dataset
Extends LB with ADaM-standard derived variables:
- BASE: baseline value
- CHG: change from baseline
- PCHG: percent change from baseline
- ABLFL: baseline flag
- ANL01FL: analysis inclusion flag
- Treatment arm and Ig subtype for subgroup analyses

---
## Clinical Variables Reference (from Methods section)

### Key Disease Biomarkers (for SCOPE model input)
| Variable | Category | Clinical Significance |
|----------|----------|----------------------|
| SPEP_MPROT | Serum Immunoglobulins | Primary disease burden marker |
| IGA / IGG | Serum Immunoglobulins | Dominant heavy chain overproduction |
| SPEP_KAPPA / SPEP_LAMBDA | Serum Immunoglobulins | Light chain disease marker |
| HGB | Hematology | Anemia (CRAB criterion) |
| CREAT / GFR | Chemistry | Renal function (CRAB criterion) |
| CA / CA_CORR | Chemistry | Hypercalcemia (CRAB criterion) |

### Treatment Regimens
| Arm | Drugs | Schedule |
|-----|-------|----------|
| IRd | Ixazomib 4mg + Lenalidomide 25mg + Dexamethasone 40mg | 28-day cycles |
| Rd  | Placebo + Lenalidomide 25mg + Dexamethasone 40mg | 28-day cycles |

---
## Usage Notes
- Time unit: 1 treatment period = 28 days ≈ 1 month (matches paper convention)
- Censoring: CNSR=0 → event observed; CNSR=1 → censored (ADTTE/DS)
- Missing lab values (~15%) reflect real-world missingness; use LOCF or model-based imputation
- For SCOPE model replication: use LB time series + DM baseline covariates + EX dosing

---
## Citation
If using this synthetic dataset, please cite the original paper:
```
Hussain Z, De Brouwer E, et al. (2024). Joint AI-driven event prediction and
longitudinal modeling in newly diagnosed and relapsed multiple myeloma.
npj Digital Medicine, 7, 200. https://doi.org/10.1038/s41746-024-01189-3
```
"""

with open(f"{out_dir}/README.md", "w") as f:
    f.write(readme)

print("✓ README written.")
