# RATIFY — Synthetic Data Branch

| Field | Value |
|-------|-------|
| Branch | `RATIFY` |
| Status | PENDING |
| NCT | NCT00651261 |
| Sponsor | Novartis |
| Indication | Newly diagnosed FLT3-positive AML (age ≤59, fit for intensive chemotherapy) |
| N | 717 (360 midostaurin+7+3 / 357 placebo+7+3) |
| Arms | Midostaurin+7+3 chemo vs Placebo+7+3 chemo |
| Primary endpoint | OS |
| OS median (active) | 74.7 mo |
| OS HR | 0.78 (95% CI 0.63–0.96) |
| Primary pub | Stone 2017 NEJM |

## Mechanism
Multi-kinase inhibitor (FLT3, PKC, VEGFR): primarily FLT3-ITD/TKD inhibition in AML
blasts. Midostaurin does not achieve complete FLT3 inhibition alone; chemo backbone dominates.

## PK Model
2-compartment oral with active metabolite CGP52421 (add separate 1-cpt metabolite PK).
Larson 2016 Blood; FDA clinical pharmacology review (NDA 207997); Weisberg 2017.

## PD Models
- `FLT3KinaseInhibition` — FLT3 occupancy (midostaurin + CGP52421) → blast kill rate
- `FribergMyelosuppression` (ANC) — chemotherapy-dominant; midostaurin adds modest increment
- FLT3-ITD allelic ratio (AR) as continuous covariate on EC50

## Unique Mechanistic Features
- Midostaurin given Days 8–21 ONLY (NOT days 1–7): avoids PK interaction with ara-C
- FLT3-ITD allelic ratio: key covariate (high AR = more FLT3-dependent → lower EC50)
- FLT3-TKD (D835) mutation: different IC50 from FLT3-ITD — must be simulated separately
- Response endpoint: AML CR (blast% <5% + ANC ≥1.0×10⁹/L + PLT ≥100×10⁹/L) — not M-protein
- SCT (stem cell transplant): ~50% proceed to transplant; model as competing event
- ANC nadir dominated by chemo (7+3 induction Day 14–21 neutropenia); midostaurin adds small increment
- Cycle structure: Induction (1–2 cycles 7+3+midostaurin) → Consolidation (HiDAC+midostaurin) → Maintenance (midostaurin monotherapy × 12 cycles)

## Dosing Schedule
- Induction: Cytarabine 200 mg/m² CIV D1–7, Daunorubicin 60 mg/m² IV D1–3,
  Midostaurin 50 mg BID PO **Days 8–21** (2 induction cycles max)
- Consolidation: HiDAC (cytarabine 3 g/m² Q12h D1,3,5) + Midostaurin 50 mg BID D8–21 (4 cycles)
- Maintenance: Midostaurin 50 mg BID D1–28 × 12 cycles (for patients not proceeding to SCT)

## Multi-Agent Status
- Literature Search: PENDING
- QSP Scientist: PENDING
- Drug Developer: PENDING
- Oncologist: PENDING
- Config: PENDING
- Validation: PENDING
