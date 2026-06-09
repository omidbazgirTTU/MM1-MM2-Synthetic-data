# RESONATE — Synthetic Data Branch

| Field | Value |
|-------|-------|
| Branch | `RESONATE` |
| Status | PENDING |
| NCT | NCT01578707 |
| Sponsor | Janssen (J&J) / Pharmacyclics |
| Indication | Relapsed/Refractory CLL or SLL |
| N | 391 (195 ibrutinib / 196 ofatumumab) |
| Arms | Ibrutinib 420mg QD (continuous) vs Ofatumumab IV |
| Primary endpoint | PFS (iwCLL) |
| PFS median (active) | Not reached (crossover at interim) |
| PFS HR | 0.22 (95% CI 0.15–0.32) |
| Primary pub | Byrd 2014 NEJM |

## Mechanism
Irreversible covalent BTK inhibitor: ibrutinib binds C481 of BTK, blocking BCR signaling.
ALC redistribution is the pathognomonic PK-PD signature.

## PK Model
1-compartment oral (Ka ~0.5 h⁻¹, CL/F ~60 L/h, Vd/F ~10,000 L, t½ ~4–6h, F~3%).
de Zwart 2016 CPT:PSP; Lippert 2017 CPT:PSP.

## PD Models
- `BTKOccupancyALC` — irreversible binding ODE: dR/dt = ksyn − kdeg×R − kon×C×R
  (R = free BTK, C = ibrutinib plasma); % occupancy drives ALC redistribution
- Two-compartment lymphocyte model: tissue ↔ blood redistribution

## Unique Mechanistic Features
- ALC redistribution MUST be reproduced: ALC rises 50–100% weeks 1–8, then declines
  (ibrutinib mobilizes CLL cells from lymph nodes to blood before killing them)
- Continuous oral dosing (QD no cycle structure) — SDTM EX must encode continuous therapy
- Ofatumumab control arm: IV mAb, 8 doses over 24 weeks (complex schedule in sdtm_ex)
- IGHV mutation status as key prognostic covariate (unmutated = worse prognosis)
- del(17p) / TP53 mutation: ibrutinib activity maintained (unlike chemotherapy)

## Dosing Schedule
- Ibrutinib: 420 mg QD continuous (no cycle structure)
- Ofatumumab: 300 mg D1, then 2000 mg weekly ×7 (Cycle 1–2), then monthly ×4

## Multi-Agent Status
- Literature Search: PENDING
- QSP Scientist: PENDING
- Drug Developer: PENDING
- Oncologist: PENDING
- Config: PENDING
- Validation: PENDING
