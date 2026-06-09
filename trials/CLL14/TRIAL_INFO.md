# CLL14 — Synthetic Data Branch

| Field | Value |
|-------|-------|
| Branch | `CLL14` |
| Status | PENDING |
| NCT | NCT02242942 |
| Sponsor | AstraZeneca / Roche |
| Indication | Treatment-naïve CLL (unfit: CIRS >6 or CrCL 30–69 mL/min) |
| N | 432 (216 Ven+Obi / 216 Clb+Obi) |
| Arms | Venetoclax+Obinutuzumab vs Chlorambucil+Obinutuzumab |
| Primary endpoint | PFS (iwCLL 2018) |
| PFS median (active) | Not reached |
| PFS HR | 0.35 (95% CI 0.23–0.53) |
| Primary pub | Fischer 2019 NEJM |

## Mechanism
BCL-2 inhibition → disruption of intrinsic apoptosis pathway → CLL cell death

## PK Model
1-compartment oral with mandatory ramp-up (20→50→100→200→400 mg, weekly escalation).
Salem 2017 CPT:PSP; Mensah 2018 CPT.

## PD Models
- `BCL2OccupancyTLS` — BCL-2 occupancy drives CLL cell kill + tumor lysis syndrome risk
- Lymphocyte kinetics ODE (ALC dynamics with treatment)

## Unique Mechanistic Features
- Fixed-duration therapy: 12 cycles (venetoclax 12 mo + obinutuzumab 6 mo)
- Ramp-up schedule is MECHANISTICALLY MANDATORY (TLS prophylaxis): simulate step-up doses
- TLS risk stratified by lymph node size + ALC at baseline
- TLS prophylaxis (allopurinol, hydration) mandatory → sdtm_cm
- Response endpoint: iwCLL criteria (ALC + lymph node SPD + cytopenias) — NOT M-protein
- MRD negativity as secondary endpoint (bone marrow + peripheral blood)

## Dosing Schedule
- Venetoclax ramp-up (6 weeks): 20 mg × 1wk → 50 mg × 1wk → 100 mg × 1wk →
  200 mg × 1wk → 400 mg starting week 5, continued to cycle 12
- Obinutuzumab IV: Cycle 1 (100mg D1, 900mg D2, 1000mg D8 and D15),
  then 1000mg D1 Cycles 2–6

## Multi-Agent Status
- Literature Search: PENDING
- QSP Scientist: PENDING
- Drug Developer: PENDING
- Oncologist: PENDING
- Config: PENDING
- Validation: PENDING
