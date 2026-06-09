---
name: trial-registry
description: Master registry of all 11 clinical trials in the multi-agent synthetic data framework. Includes branch name, status, mechanism, published sources, and validation target summary.
metadata:
  type: reference
---

# Trial Registry

## Framework Overview
11 Phase 3 oncology trials. MM1 and MM2 share the same mechanism (ixazomib) and are
retained as the validated seed pair. All other 9 trials use mechanistically distinct
drug classes. Main branch is the merge target; each trial has its own branch.

```
main ← MM1-MM2 (complete, 68/68 PASS)
     ← POLLUX
     ← CLL14
     ← RESONATE
     ← MONARCH-2
     ← EMILIA
     ← FLAURA
     ← SOLO-1
     ← KEYNOTE-189
     ← RATIFY
```

---

## Trial 1 & 2 — TOURMALINE-MM1 / TOURMALINE-MM2

| Field | MM1 | MM2 |
|-------|-----|-----|
| Branch | `MM1-MM2` | `MM1-MM2` |
| Status | **COMPLETE — 68/68 PASS** | **COMPLETE — 68/68 PASS** |
| NCT | NCT01564537 | NCT01850524 |
| Sponsor | Takeda | Takeda |
| Indication | RRMM | NDMM |
| N | 722 (360/362) | 705 (351/354) |
| Arms | IRd vs Rd | IRd vs Rd |
| Primary endpoint | PFS | PFS |
| PFS median (active) | 20.6 mo | 35.3 mo |
| PFS HR | 0.74 | 0.83 |
| Primary pub | Moreau 2016 NEJM | Kumar 2023 Leukemia |

**Mechanism**: Reversible 20S proteasome inhibition (oral, boronic acid)
**PK model**: 3-compartment oral (Gupta 2017)
**PD models**: TwoPopulationMprotein (Srimani 2022) + FribergMyelosuppression (PLT cumulative AUC)
**Unique features**: AR(1) M-protein residuals (ρ=0.60), PLT nadir IOV, Ka IOV, pk_series pipeline

---

## Trial 3 — POLLUX

| Field | Value |
|-------|-------|
| Branch | `POLLUX` |
| Status | PENDING |
| NCT | NCT02076009 |
| Sponsor | Janssen (J&J) |
| Indication | RRMM |
| N | 569 (286/283) |
| Arms | DRd (daratumumab+Rd) vs Rd |
| Primary endpoint | PFS |
| PFS median (active) | Not reached (>44 mo) |
| PFS HR | 0.37 (95% CI 0.27–0.52) |
| Primary pub | Dimopoulos 2016 NEJM |

**Mechanism**: CD38-targeted monoclonal antibody; target-mediated disposition (TMDD)
**PK model**: 2-compartment IV + TMDD QSS (Gibiansky 2010)
**PD models**: TMDDModel + TwoPopulationMprotein (free drug drives M-protein ODE) + NK cell depletion arm
**Key published PK**: Xu 2017 CPT:PSP (daratumumab popPK in MM); Kd ~10 nM, CD38 target baseline IIV
**Unique features**:
- Free drug (not total) is the active PK driver — must implement TMDD
- First infusion reaction rate ~50% (encode as AE Cycle 1 Day 1)
- Pre-medication (antihistamine, antipyretic, dexamethasone) mandatory → sdtm_cm
- CD38 expression as patient-level IIV covariate on Kd

---

## Trial 4 — CLL14

| Field | Value |
|-------|-------|
| Branch | `CLL14` |
| Status | PENDING |
| NCT | NCT02242942 |
| Sponsor | AstraZeneca / Roche |
| Indication | Treatment-naïve CLL (unfit patients, CIRS >6 or CrCL 30–69) |
| N | 432 (216/216) |
| Arms | Ven+Obi (venetoclax+obinutuzumab) vs Clb+Obi (chlorambucil+obinutuzumab) |
| Primary endpoint | PFS (iwCLL 2018) |
| PFS median (active) | Not reached |
| PFS HR | 0.35 (95% CI 0.23–0.53) |
| Primary pub | Fischer 2019 NEJM |

**Mechanism**: BCL-2 inhibition → intrinsic apoptosis pathway activation
**PK model**: 1-compartment oral with ramp-up schedule (weekly dose escalation)
**PD models**: BCL2OccupancyTLS + lymphocyte kinetics ODE
**Key published PK**: Salem 2017 CPT:PSP; Mensah 2018 CPT (venetoclax popPK in CLL)
**Unique features**:
- Fixed-duration therapy: 12 cycles total (unlike ibrutinib which is continuous)
- Mandatory ramp-up (20→50→100→200→400mg weekly): mechanistically driven by TLS risk
- TLS prophylaxis (allopurinol, hydration) → sdtm_cm
- Response endpoint: iwCLL criteria (ALC + lymph node SPD + cytopenias, not M-protein)
- MRD negativity rate as secondary endpoint

---

## Trial 5 — RESONATE

| Field | Value |
|-------|-------|
| Branch | `RESONATE` |
| Status | PENDING |
| NCT | NCT01578707 |
| Sponsor | Janssen (J&J) / Pharmacyclics |
| Indication | Relapsed/Refractory CLL or SLL |
| N | 391 (195/196) |
| Arms | Ibrutinib vs Ofatumumab |
| Primary endpoint | PFS (iwCLL) |
| PFS HR | 0.22 (95% CI 0.15–0.32) |
| Primary pub | Byrd 2014 NEJM |

**Mechanism**: Irreversible covalent BTK inhibitor; BTK occupancy model + ALC redistribution
**PK model**: 1-compartment oral (ibrutinib absorbed rapidly; t½ ~4–6h; high first-pass)
**PD models**: BTKOccupancyALC (irreversible binding ODE + two-compartment lymphocyte)
**Key published PK**: de Zwart 2016 CPT:PSP; Lippert 2017 CPT:PSP
**Unique features**:
- ALC redistribution MUST be modeled: ALC rises 50–100% in weeks 1–8, then declines
- Continuous oral dosing (420mg QD, no cycle structure) — unique SDTM EX encoding
- Ofatumumab control arm: IV mAb, 8 doses over 24 weeks (complex schedule)
- IGHV mutation status as key prognostic covariate

---

## Trial 6 — MONARCH-2

| Field | Value |
|-------|-------|
| Branch | `MONARCH-2` |
| Status | PENDING |
| NCT | NCT02107703 |
| Sponsor | Eli Lilly |
| Indication | HR+ HER2− metastatic breast cancer (endocrine-resistant) |
| N | 669 (446/223) [2:1 randomization] |
| Arms | Abemaciclib+Fulvestrant vs Placebo+Fulvestrant |
| Primary endpoint | PFS (RECIST 1.1) |
| PFS median (active) | 16.4 mo |
| PFS HR | 0.553 (95% CI 0.449–0.681) |
| Primary pub | Sledge 2017 JCO |

**Mechanism**: Selective CDK4/6 inhibitor → G1 cell cycle arrest in tumor cells
**PK model**: 2-compartment oral (abemaciclib + active metabolites M2, M20) — Tate 2018 CPT:PSP
**PD models**: TumorGrowthInhibition (Simeoni 2004) + FribergMyelosuppression (ANC, CDK4/6-specific)
**Key published PK**: Tate 2018 CPT:PSP; Bhansali 2019
**Unique features**:
- Abemaciclib dosed CONTINUOUSLY (150mg BID, no days off) unlike palbociclib/ribociclib
- Active metabolites M2 and M20: combined parent+metabolite drives PD
- Diarrhea is dominant early toxicity (Cycles 1–2, up to 86% any grade) — must model GI AE
- ANC nadir timing different from MM drugs: Day 14–21 (not Day 11–15)
- ESR1 mutation as acquired resistance biomarker (optional second-line context)

---

## Trial 7 — EMILIA

| Field | Value |
|-------|-------|
| Branch | `EMILIA` |
| Status | PENDING |
| NCT | NCT00829166 |
| Sponsor | Genentech / Roche |
| Indication | HER2+ metastatic breast cancer (second-line, post trastuzumab+taxane) |
| N | 991 (495/496) |
| Arms | T-DM1 vs Capecitabine+Lapatinib |
| Primary endpoint | PFS + OS (co-primary) |
| PFS median (active) | 9.6 mo |
| OS median (active) | 30.9 mo |
| PFS HR | 0.650 (95% CI 0.55–0.77) |
| Primary pub | Verma 2012 NEJM |

**Mechanism**: HER2-targeted ADC; antibody delivers DM1 (microtubule inhibitor) intracellularly
**PK model**: TMDD for T-DM1 antibody + intracellular DM1 payload compartment
**PD models**: ADCPayloadRelease + FribergMyelosuppression (PLT Day 8 nadir) + peripheral neuropathy
**Key published PK**: Bender 2012 CPT; Girish 2012; Lu 2016
**Unique features**:
- Two-level PK: plasma T-DM1 (TMDD, HER2 binding) + intracellular DM1 (payload release)
- DAR (drug-antibody ratio): T-DM1 average DAR = 3.5 DM1 molecules per antibody
- PLT nadir Day 8 (faster than PI Day 11–15) — DM1 direct MK precursor toxicity
- LVEF monitoring every 3 cycles (T-DM1 has mild cardiotoxicity from HER2 blockade)
- Peripheral neuropathy: cumulative DM1 component, worsens with prior taxane

---

## Trial 8 — FLAURA

| Field | Value |
|-------|-------|
| Branch | `FLAURA` |
| Status | PENDING |
| NCT | NCT02296125 |
| Sponsor | AstraZeneca |
| Indication | Treatment-naïve EGFR-mutated advanced NSCLC (del19 or L858R) |
| N | 556 (279/277) |
| Arms | Osimertinib 80mg QD vs Standard EGFR TKI (gefitinib or erlotinib) |
| Primary endpoint | PFS (RECIST 1.1) |
| PFS median (active) | 18.9 mo |
| PFS HR | 0.46 (95% CI 0.37–0.57) |
| Primary pub | Soria 2018 NEJM |

**Mechanism**: 3rd-generation EGFR TKI with CNS penetration + T790M coverage
**PK model**: 2-compartment oral with high BBB penetration (CNS compartment)
**PD models**: TumorGrowthInhibition (primary lesion) + separate CNS-TGI + resistance ODE
**Key published PK**: Dickinson 2016 CPT:PSP; Vishwanathan 2019; FDA briefing document
**Unique features**:
- CNS compartment: osimertinib achieves CNS/plasma ratio ~15× higher than 1st-gen TKIs
- EGFR mutation type as covariate: del19 vs L858R (del19 has better prognosis)
- Acquired resistance ODE: C797S mutation emergence rate ∝ cumulative AUC_osimertinib
- Skin rash as PD biomarker of EGFR target engagement (higher rash ↔ better tumor response)
- QTc monitoring: osimertinib prolongs QTc (model as separate cardiac AE)

---

## Trial 9 — SOLO-1

| Field | Value |
|-------|-------|
| Branch | `SOLO-1` |
| Status | PENDING |
| NCT | NCT01844986 |
| Sponsor | AstraZeneca |
| Indication | BRCA1/2-mutated advanced ovarian cancer (first-line maintenance post-platinum) |
| N | 391 (260/131) [2:1 randomization] |
| Arms | Olaparib 300mg BID vs Placebo (maintenance) |
| Primary endpoint | PFS |
| PFS median (active) | Not reached (51.8 mo at 7-yr follow-up) |
| PFS HR | 0.33 (95% CI 0.25–0.43) |
| Primary pub | Moore 2018 NEJM |

**Mechanism**: PARP1/2 inhibition → synthetic lethality in BRCA1/2-deficient tumor cells
**PK model**: 2-compartment oral (tablet formulation 300mg BID); Plummer 2020 CPT:PSP
**PD models**: SyntheticLethalityER (BRCA-gated PARP inhibition → tumor growth suppression)
**Key published PK**: Plummer 2020 CPT:PSP; Menear 2008 JMedChem
**Unique features**:
- Maintenance setting: patients are in CR/PR after platinum; baseline SLD may be near 0
- BRCA status (germline vs somatic) is the defining biomarker: modifies EC50 in simulation
- HRD score as optional continuous covariate within BRCA-WT patients
- MDS/AML risk: small but serious; model as rare late-onset AE (rate ~0.5%)
- Dose reductions (300→250→200mg BID): primarily hematologic (anemia, ANC)
- Anemia is common and often Grade ≥3: model separately from oncology-typical ANC nadir

---

## Trial 10 — KEYNOTE-189

| Field | Value |
|-------|-------|
| Branch | `KEYNOTE-189` |
| Status | PENDING |
| NCT | NCT02578680 |
| Sponsor | Merck |
| Indication | Metastatic non-squamous NSCLC without EGFR/ALK alteration (first-line) |
| N | 616 (410/206) [2:1 randomization] |
| Arms | Pembrolizumab+Pemetrexed+Platinum vs Placebo+Pemetrexed+Platinum |
| Primary endpoint | OS + PFS (co-primary) |
| PFS median (active) | 9.0 mo |
| OS median (active) | 22.0 mo |
| PFS HR | 0.48 (95% CI 0.40–0.58) |
| Primary pub | Gandhi 2018 NEJM |

**Mechanism**: PD-1 checkpoint blockade → T-cell re-activation; delayed kinetics + PD-L1 covariate
**PK model**: 2-compartment IV (pembrolizumab Q3W dosing); Li 2017 CPT:PSP
**PD models**: TumorGrowthInhibition (TGI with immunotherapy lag) + T-cell activation ODE
**Key published PK**: Li 2017 CPT:PSP; Gibiansky 2020; Freshwater 2017
**Unique features**:
- Delayed response: implement Tlag ≈ 6–12 weeks before immune-mediated tumor kill begins
- PD-L1 TPS as continuous covariate: TPS ≥50% → lower EC50 → deeper/faster response
- Hyperprogression: ~10% of patients show early PD (kgrowth × 2 in first 2 cycles)
- irAE (immune-related AEs): pneumonitis, colitis, hepatitis — steroid management
- Combination with pemetrexed+platinum: chemotherapy adds direct tumor kill (separate arm)
- No dose modifications for pembrolizumab (discontinue only)

---

## Trial 11 — RATIFY

| Field | Value |
|-------|-------|
| Branch | `RATIFY` |
| Status | PENDING |
| NCT | NCT00651261 |
| Sponsor | Novartis |
| Indication | Newly diagnosed FLT3-positive AML (age ≤59, fit for intensive chemotherapy) |
| N | 717 (360/357) |
| Arms | Midostaurin+7+3 chemo vs Placebo+7+3 chemo |
| Primary endpoint | OS |
| OS median (active) | 74.7 mo |
| OS HR | 0.78 (95% CI 0.63–0.96) |
| Primary pub | Stone 2017 NEJM |

**Mechanism**: Multi-kinase inhibitor (FLT3, PKC, VEGFR) — primarily FLT3-ITD/TKD inhibition in AML blasts
**PK model**: 2-compartment oral with active metabolite CGP52421; Larson 2016 Blood
**PD models**: FLT3KinaseInhibition + FribergMyelosuppression (ANC, chemotherapy-dominant)
**Key published PK**: Weisberg 2017; FDA clinical pharmacology review (NDA 207997)
**Unique features**:
- Given Days 8–21 only (NOT on chemo days 1–7 to avoid PK interaction with ara-C)
- FLT3-ITD allelic ratio (AR): key continuous covariate on EC50 (high AR = more FLT3-dependent)
- FLT3-TKD mutations also included (D835): different IC50 from FLT3-ITD
- Response endpoint: AML CR rate (blast% < 5% + ANC ≥1.0 + PLT ≥100) — not M-protein or SLD
- SCT (stem cell transplant): ~50% proceed; model as competing event or censoring at SCT
- ANC toxicity dominated by chemotherapy component (not midostaurin): important for calibration
- Cycle structure: induction (1–2 cycles 7+3+midostaurin) → consolidation (HiDAC + midostaurin) → maintenance (midostaurin monotherapy 12 cycles)

---

## Mechanism Tree (Quick Reference)

```
ORAL SMALL MOLECULES
├── Reversible 20S proteasome inhibitor    → MM1/MM2  (Takeda)
├── BTK irreversible covalent              → RESONATE (J&J)
├── CDK4/6 kinase inhibitor                → MONARCH-2 (Lilly)
├── EGFR TKI + CNS penetration            → FLAURA   (AZ)
├── PARP inhibitor / synthetic lethality  → SOLO-1   (AZ)
└── FLT3 multi-kinase inhibitor            → RATIFY   (Novartis)

IV MONOCLONAL ANTIBODIES
├── CD38 TMDD (naked mAb)                 → POLLUX   (J&J)
└── PD-1 checkpoint + T-cell activation  → KEYNOTE-189 (Merck)

ANTIBODY-DRUG CONJUGATE
└── HER2 TMDD + intracellular payload    → EMILIA   (Genentech)

BCL-2 FAMILY INHIBITOR (ORAL)
└── BCL-2 occupancy + TLS                → CLL14    (AZ/Roche)
```

## Company Coverage
| Company | Trials |
|---------|--------|
| Takeda | MM1, MM2 |
| J&J (Janssen) | POLLUX, RESONATE |
| AstraZeneca | CLL14, FLAURA, SOLO-1 |
| Eli Lilly | MONARCH-2 |
| Genentech/Roche | EMILIA |
| Merck | KEYNOTE-189 |
| Novartis | RATIFY |

## Indication Coverage
| Indication | Trials |
|-----------|--------|
| Multiple Myeloma | MM1, MM2, POLLUX |
| CLL/SLL | CLL14, RESONATE |
| HR+ Breast Cancer | MONARCH-2 |
| HER2+ Breast Cancer | EMILIA |
| NSCLC | FLAURA, KEYNOTE-189 |
| Ovarian Cancer | SOLO-1 |
| AML | RATIFY |
