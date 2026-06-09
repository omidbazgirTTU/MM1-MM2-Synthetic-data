# FLAURA — Synthetic Data Branch

| Field | Value |
|-------|-------|
| Branch | `FLAURA` |
| Status | PENDING |
| NCT | NCT02296125 |
| Sponsor | AstraZeneca |
| Indication | Treatment-naïve EGFR-mutated advanced NSCLC (del19 or L858R) |
| N | 556 (279 osimertinib / 277 std EGFR TKI) |
| Arms | Osimertinib 80mg QD vs Standard EGFR TKI (gefitinib or erlotinib) |
| Primary endpoint | PFS (RECIST 1.1) |
| PFS median (active) | 18.9 mo |
| PFS HR | 0.46 (95% CI 0.37–0.57) |
| Primary pub | Soria 2018 NEJM |

## Mechanism
3rd-generation irreversible EGFR TKI: covalent binding to C797S+T790M; superior
CNS penetration vs 1st-generation TKIs (CNS/plasma ratio ~15×).

## PK Model
2-compartment oral + CNS compartment (BBB penetration ODE).
Dickinson 2016 CPT:PSP; Vishwanathan 2019; FDA briefing document.

## PD Models
- `TumorGrowthInhibition` (primary lesion, Simeoni 2004) — EGFR occupancy drives TGI
- Separate CNS-TGI module: CNS drug concentration drives intracranial lesion control
- Resistance ODE: C797S mutation emergence rate ∝ cumulative AUC_osimertinib

## Unique Mechanistic Features
- CNS compartment is mechanistically required: osimertinib controls brain mets (1st-gen TKIs do not)
- EGFR mutation type as covariate: del19 has better prognosis than L858R (lower EC50)
- Acquired resistance: C797S emergence modeled as time-to-event (AUC-dependent rate)
- Skin rash as PD biomarker of EGFR target engagement (rash severity ↔ tumor response)
- QTc prolongation: model as cardiac AE (osimertinib QTc effect, grade-based management)
- CNS progression as competing event for PFS analysis

## Dosing Schedule
- Osimertinib: 80 mg QD continuous (dose reduction to 40 mg for toxicity)
- Standard TKI arm: gefitinib 250 mg QD or erlotinib 150 mg QD (investigator choice)

## Multi-Agent Status
- Literature Search: PENDING
- QSP Scientist: PENDING
- Drug Developer: PENDING
- Oncologist: PENDING
- Config: PENDING
- Validation: PENDING
