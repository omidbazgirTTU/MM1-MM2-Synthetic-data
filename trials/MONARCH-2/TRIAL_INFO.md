# MONARCH-2 — Synthetic Data Branch

| Field | Value |
|-------|-------|
| Branch | `MONARCH-2` |
| Status | PENDING |
| NCT | NCT02107703 |
| Sponsor | Eli Lilly |
| Indication | HR+ HER2− metastatic breast cancer (endocrine-resistant) |
| N | 669 (446 abemaciclib / 223 placebo) [2:1 randomization] |
| Arms | Abemaciclib+Fulvestrant vs Placebo+Fulvestrant |
| Primary endpoint | PFS (RECIST 1.1) |
| PFS median (active) | 16.4 mo |
| PFS HR | 0.553 (95% CI 0.449–0.681) |
| Primary pub | Sledge 2017 JCO |

## Mechanism
Selective CDK4/6 inhibitor → G1 phase cell cycle arrest in Rb-positive tumor cells.
Abemaciclib is CONTINUOUSLY dosed (unlike palbociclib/ribociclib with 7-day breaks).

## PK Model
2-compartment oral (abemaciclib + active metabolites M2 and M20 modeled explicitly).
Tate 2018 CPT:PSP; Bhansali 2019.

## PD Models
- `TumorGrowthInhibition` (Simeoni 2004) — CDK4/6 inhibition → G1 arrest
- `FribergMyelosuppression` (ANC) — CDK4/6-specific nadir timing (Day 14–21)

## Unique Mechanistic Features
- Active metabolites M2 and M20: combined parent+metabolite PK drives PD
- Continuous BID dosing (150mg BID): no drug holiday unlike palbociclib/ribociclib
- Diarrhea dominant early toxicity (Cycles 1–2, up to 86% any grade): model as GI AE
- ANC nadir Day 14–21 (later than proteasome inhibitor Day 11–15)
- Response endpoint: RECIST 1.1 SLD% change (not M-protein)
- ESR1 mutation: acquired resistance biomarker (optional, second-line context)
- Visceral involvement (liver/lung mets) at baseline as prognostic covariate

## Dosing Schedule
- Abemaciclib: 150 mg BID continuous (28-day cycles, no off-days)
- Fulvestrant: 500 mg IM D1+D15 Cycle 1, then D1 each cycle

## Multi-Agent Status
- Literature Search: PENDING
- QSP Scientist: PENDING
- Drug Developer: PENDING
- Oncologist: PENDING
- Config: PENDING
- Validation: PENDING
