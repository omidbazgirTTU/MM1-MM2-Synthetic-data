---
name: tourmaline-data-generation-workflow
description: >
  End-to-end workflow for generating, validating, and extending the TOURMALINE
  MM1/MM2 synthetic CDISC datasets (SDTM + ADaM). Use this when adding new
  domains, re-running generation, or checking which script does what.
applies_when:
  - Adding a new SDTM or ADaM domain
  - Re-running data generation after code changes
  - Looking up which file owns a given variable or domain
  - Checking generation order or command syntax
  - Understanding CDISC conventions used in this project
keywords:
  - generate_v2.py, generate_pk_v2.py, validate_data.py
  - SDTM, ADaM, CDISC, sdtm_dm, sdtm_ex, sdtm_lb, adam_adsl, adam_adlb, adam_adtte
  - MM1, MM2, TOURMALINE, IRd, Rd, NDMM, RRMM
load_cost: low   # lightweight procedural reference
---

## Level 2 — Instructions

### Commands (run in order)
```bash
source .venv/bin/activate
python3 scripts/generate_v2.py MM2 MM1     # SDTM + ADaM (both studies)
python3 scripts/validate_data.py MM2 MM1   # figures + pass/fail table
python3 scripts/plot_individual_patients.py MM2 MM1  # 1427 patient figures
```
`generate_pk_v2.py` is called internally; run separately only when tuning PK params.

### Domain Ownership
| Domain | Function | Key outputs |
|--------|----------|-------------|
| DM | `make_dm()` | weight, height, BMI, BSA, age, sex, race, ISS, ECOG |
| EX | `make_ex()` | dose records, dose mods, adherence |
| CM | `make_cm()` | CYP3A4 drugs, anticoagulants, G-CSF |
| LB | `make_lb()` | 34 lab tests × all cycles × 3 weekly timepoints |
| DS/ADTTE | `make_ds()`, `make_adtte()` | PFS/OS (Weibull, SURV_RNG=77) |
| ADSL | `make_adsl()` | ISS, R-ISS, cytogenetics, CrCL, BSA, baseline labs |
| ADLB | `make_adlb()` | CHG, PCHG, ABLFL, EPOCH from sdtm_lb |
| PC/PP | `generate_pk_v2.py` | PK concentrations + NCA |

### CDISC Conventions
- **USUBJID**: `TOURMALINE-MM2-0001` (no site prefix)
- **EPOCH**: `"BASELINE"` for cycle index 0, `"TREATMENT"` for all others
- **VISITNUM**: cycle index + 1 (1-indexed; 1=baseline, 2=cycle 1, …)
- **WEEKNUM**: 1/2/3 = Day 1/8/15 within cycle (hematology/chemistry only)
- **ABLFL**: `"Y"` on the baseline record per PARAMCD per USUBJID
- **CNSR**: 0 = event, 1 = censored (ADTTE/DS)
- **ADY**: study day; ADY=1 = day of first dose

### Assessment Schedule (per SoA)
| Test type | When assessed |
|-----------|---------------|
| SPEP/Ig/UPEP | Day 1 of each cycle (every 4 weeks) |
| B2MG | Every 12 weeks (every 3rd cycle Day 1) |
| Hematology (PLT, HGB, NEUT, WBC) | Days 1, 8, 15 every cycle |
| Chemistry | Days 1, 8, 15 every cycle |

### Arm Randomization (exact counts)
```python
arm_arr = np.array(["IRd"] * n_ird + ["Rd"] * n_rd, dtype=object)
RNG.shuffle(arm_arr)
# MM2: n_ird=351, n_rd=354  |  MM1: n_ird=360, n_rd=362
```

### Current Validation Status (46/48 PASS — 2026-04)
Failures: Ixazomib Cmax +30% (MM2) / +38% (MM1) — Gupta 2017 model overprediction;
not fixable without ODE refactor.

---

## Level 3 — Resources

### Key file locations
| File | Purpose |
|------|---------|
| `scripts/generate_v2.py` | Main generator — all domains |
| `scripts/generate_pk_v2.py` | PK concentrations + NCA |
| `scripts/validate_data.py` | Validation figures + pass/fail |
| `scripts/plot_individual_patients.py` | Per-patient longitudinal plots |
| `outputs/tables/validation_summary.csv` | Current pass/fail per metric |
| `outputs/tables/response_rates.csv` | ORR/VGPR/CR by arm × study |

### Published targets reference
- MM2 (NDMM): Moreau et al. 2019 NEJM (NCT01850524)
- MM1 (RRMM): Moreau et al. 2016 Lancet Oncol (NCT01564537)
- PK: Gupta et al. 2017 Clin Pharmacokinet (3-cmt Ixazomib population PK)
- PD validation method: Hussain Z et al. 2024 npj Digital Medicine

### Target summary
| Metric | MM2 IRd | MM2 Rd | MM1 IRd | MM1 Rd |
|--------|---------|--------|---------|--------|
| PFS median (mo) | 35.3 | 21.8 | 20.6 | 14.7 |
| OS median (mo) | 60.0 | 48.0 | 53.6 | 51.6 |
| ORR (%) | 82 | 75 | 78 | 72 |
| VGPR+ (%) | 63 | 55 | 48 | 39 |
| CR+ (%) | 28 | 14 | 12 | 7 |
| Grade 3 PLT (%) | 25 | 14 | 31 | 16 |
