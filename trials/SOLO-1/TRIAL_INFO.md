# SOLO-1 — Synthetic Data Branch

| Field | Value |
|-------|-------|
| Branch | `SOLO-1` |
| Status | PENDING |
| NCT | NCT01844986 |
| Sponsor | AstraZeneca |
| Indication | BRCA1/2-mutated advanced ovarian cancer (1st-line maintenance post-platinum) |
| N | 391 (260 olaparib / 131 placebo) [2:1 randomization] |
| Arms | Olaparib 300mg BID maintenance vs Placebo |
| Primary endpoint | PFS |
| PFS median (active) | Not reached (51.8 mo at 7-yr follow-up) |
| PFS HR | 0.33 (95% CI 0.25–0.43) |
| Primary pub | Moore 2018 NEJM |

## Mechanism
PARP1/2 inhibition → synthetic lethality in BRCA1/2-deficient tumor cells.
BRCA-mutated cells cannot repair double-strand breaks via HR; PARP trapping causes replication
fork collapse → tumor-selective cell death.

## PK Model
2-compartment oral (tablet, 300mg BID). Plummer 2020 CPT:PSP; Menear 2008.

## PD Models
- `SyntheticLethalityER` — PARP occupancy × BRCA status → effective tumor kill rate
  (BRCA mut: EC50 ↓ 10×; BRCA WT: minimal effect)
- `FribergMyelosuppression` (ANC + HGB) — anemia is common Grade ≥3

## Unique Mechanistic Features
- Maintenance setting: patients in CR/PR post-platinum; baseline SLD near 0 or undetectable
- BRCA1/2 status is the defining biomarker: germline vs somatic must be tracked separately
- HRD score as optional continuous covariate within BRCA-WT patients
- Synthetic lethality: BRCA status gates the EC50 — BRCA-mut patients respond; BRCA-WT do not
- MDS/AML risk: rare late-onset AE (~0.5%); model as low-rate competing event
- Anemia often Grade ≥3 (model separately from standard ANC Friberg)
- Dose reductions: 300 mg → 250 mg → 200 mg BID (primarily for hematologic toxicity)

## Dosing Schedule
- Olaparib tablet: 300 mg BID continuous (starting ≤8 weeks after last platinum dose)
- Duration: until disease progression or unacceptable toxicity (max 2 years in original protocol)

## Multi-Agent Status
- Literature Search: PENDING
- QSP Scientist: PENDING
- Drug Developer: PENDING
- Oncologist: PENDING
- Config: PENDING
- Validation: PENDING
