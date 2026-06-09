# KEYNOTE-189 — Synthetic Data Branch

| Field | Value |
|-------|-------|
| Branch | `KEYNOTE-189` |
| Status | PENDING |
| NCT | NCT02578680 |
| Sponsor | Merck |
| Indication | Metastatic non-squamous NSCLC without EGFR/ALK alteration (1st-line) |
| N | 616 (410 pembrolizumab+chemo / 206 placebo+chemo) [2:1 randomization] |
| Arms | Pembrolizumab+Pemetrexed+Platinum vs Placebo+Pemetrexed+Platinum |
| Primary endpoint | OS + PFS (co-primary) |
| PFS median (active) | 9.0 mo |
| OS median (active) | 22.0 mo |
| PFS HR | 0.48 (95% CI 0.40–0.58) |
| Primary pub | Gandhi 2018 NEJM |

## Mechanism
PD-1 checkpoint blockade → T-cell re-activation → immune-mediated tumor kill.
Delayed kinetics (immune activation lag) and PD-L1 expression govern response depth.

## PK Model
2-compartment IV, Q3W flat dosing (200 mg) — receptor-mediated clearance.
Li 2017 CPT:PSP; Gibiansky 2020; Freshwater 2017.

## PD Models
- `TumorGrowthInhibition` (TGI) with immunotherapy lag: Tlag ≈ 6–12 weeks
- T-cell activation ODE: PD-1 occupancy → effector T-cell expansion → tumor kill rate
- Chemotherapy component (pemetrexed+platinum) adds direct tumor kill (parallel model)

## Unique Mechanistic Features
- Delayed response: immune activation lag ~6–12 weeks before tumor regression begins
- PD-L1 TPS as continuous covariate: TPS ≥50% → lower EC50 → deeper/faster response
- Hyperprogression: ~10% patients show early PD (kgrowth × 2 in first 2 cycles)
- irAE (immune-related AEs): pneumonitis, colitis, hepatitis — steroid management required
- Pembrolizumab has NO dose modifications (only discontinuation for severe irAEs)
- PD-L1 score must be generated as continuous covariate correlated with TPS endpoint
- Combination arm: two PK/PD models active simultaneously (additive kill assumption)

## Dosing Schedule
- Pembrolizumab: 200 mg IV Q3W (flat dosing, 35 cycles max = ~2 years)
- Pemetrexed: 500 mg/m² IV Q3W (premedicated: folic acid, B12, dexamethasone)
- Platinum: Carboplatin AUC5 or cisplatin 75 mg/m² Q3W × 4 cycles then stop

## Multi-Agent Status
- Literature Search: PENDING
- QSP Scientist: PENDING
- Drug Developer: PENDING
- Oncologist: PENDING
- Config: PENDING
- Validation: PENDING
