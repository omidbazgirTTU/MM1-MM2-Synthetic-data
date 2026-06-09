# EMILIA — Synthetic Data Branch

| Field | Value |
|-------|-------|
| Branch | `EMILIA` |
| Status | PENDING |
| NCT | NCT00829166 |
| Sponsor | Genentech / Roche |
| Indication | HER2+ metastatic breast cancer (2nd-line, post trastuzumab+taxane) |
| N | 991 (495 T-DM1 / 496 capecitabine+lapatinib) |
| Arms | T-DM1 3.6 mg/kg IV Q3W vs Capecitabine+Lapatinib |
| Primary endpoint | PFS + OS (co-primary) |
| PFS median (active) | 9.6 mo |
| OS median (active) | 30.9 mo |
| PFS HR | 0.650 (95% CI 0.55–0.77) |
| Primary pub | Verma 2012 NEJM |

## Mechanism
HER2-targeted antibody-drug conjugate (ADC): trastuzumab antibody delivers DM1
(microtubule inhibitor) via lysosomal degradation after HER2-mediated internalization.

## PK Model
Two-level PK system:
1. T-DM1 plasma: TMDD (HER2 binding in tumor) + FcRn recycling
2. Intracellular DM1 payload: release compartment (lysosomal degradation → free DM1)
Bender 2012 CPT; Girish 2012; Lu 2016.

## PD Models
- `ADCPayloadRelease` — payload release model (DAR degradation, DM1 intracellular conc)
- `FribergMyelosuppression` (PLT) — DM1 direct megakaryocyte precursor toxicity, **Day 8 nadir**
- Peripheral neuropathy ODE — cumulative DM1, worsened by prior taxane

## Unique Mechanistic Features
- Two-level PK: T-DM1 (antibody PK) + DM1 payload (small molecule PK) — both required
- DAR (drug-antibody ratio): T-DM1 average DAR = 3.5; simulate DAR deconjugation kinetics
- PLT nadir Day 8 (NOT Day 11–15): DM1 direct MK precursor mechanism (faster than PI)
- LVEF monitoring every 3 cycles (cardiotoxicity from HER2 blockade)
- Prior taxane therapy as covariate on peripheral neuropathy baseline risk
- HER2 IHC score (3+ vs 2+/FISH+) as PK covariate on TMDD Kd

## Dosing Schedule
- T-DM1: 3.6 mg/kg IV over 90 min (D1) Q21D
- Capecitabine: 1000 mg/m² BID D1–D14 Q21D; Lapatinib 1250 mg QD continuous

## Multi-Agent Status
- Literature Search: PENDING
- QSP Scientist: PENDING
- Drug Developer: PENDING
- Oncologist: PENDING
- Config: PENDING
- Validation: PENDING
