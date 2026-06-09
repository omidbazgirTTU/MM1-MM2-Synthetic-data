# POLLUX — Synthetic Data Branch

| Field | Value |
|-------|-------|
| Branch | `POLLUX` |
| Status | PENDING |
| NCT | NCT02076009 |
| Sponsor | Janssen (J&J) |
| Indication | Relapsed/Refractory Multiple Myeloma (RRMM) |
| N | 569 (286 DRd / 283 Rd) |
| Arms | DRd (daratumumab+lenalidomide+dex) vs Rd |
| Primary endpoint | PFS |
| PFS median (active) | Not reached (>44 mo) |
| PFS HR | 0.37 (95% CI 0.27–0.52) |
| Primary pub | Dimopoulos 2016 NEJM |

## Mechanism
CD38-targeted monoclonal antibody; target-mediated disposition (TMDD)

## PK Model
2-compartment IV + TMDD QSS (Gibiansky 2010). Free drug (not total) drives PD.

## PD Models
- `TMDDModel` — daratumumab binding to CD38+ tumor cells
- `TwoPopulationMprotein` — free drug drives M-protein ODE
- NK cell depletion arm (CD38 also expressed on NK cells)

## Key Published PK Parameters
- Xu 2017 CPT:PSP (daratumumab popPK in MM); Kd ~10 nM
- CD38 target baseline with IIV as patient-level covariate

## Unique Mechanistic Features
- Free drug vs total drug distinction — TMDD is mandatory, not optional
- First infusion reaction rate ~50%: encode as Cycle 1 Day 1 AE
- Pre-medication (antihistamine, antipyretic, dexamethasone) → sdtm_cm
- CD38 expression as IIV covariate on Kd (binding affinity)
- CD38 depletion on NK cells: immune suppression side effect

## Dosing Schedule
- Daratumumab IV: 16 mg/kg weekly Cycles 1–2, biweekly Cycles 3–6, monthly thereafter
- Lenalidomide 25 mg Days 1–21, Dexamethasone 40 mg weekly (28-day cycles)

## Multi-Agent Status
- Literature Search: PENDING
- QSP Scientist: PENDING
- Drug Developer: PENDING
- Oncologist: PENDING
- Config: PENDING
- Validation: PENDING
