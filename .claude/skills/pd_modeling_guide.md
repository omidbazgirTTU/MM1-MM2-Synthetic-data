---
name: pd-modeling-guide
description: Reference for all PD model classes in the engine library. For each model: ODE structure, parameters, mechanistic dependencies, biomarker output, and which trials use it.
metadata:
  type: reference
---

# PD Modeling Guide — Engine Library

## Design Principles
1. Every PD model takes PK output (Cp, AUC, occupancy, cell count) as its mechanistic input.
2. Every safety endpoint (PLT, ANC, cardiac) can feed back into dose modification logic.
3. PD models are instantiated from the trial YAML config — no hardcoding in model classes.
4. IIV on PD parameters is log-normal. IOV where applicable.
5. AR(1) autocorrelation on longitudinal biomarker residuals where clinically justified.

---

## Model 1: TwoPopulationMprotein
**Trials**: TOURMALINE-MM1/MM2, POLLUX (daratumumab, same ODE different drug input)
**Endpoint**: Serum M-protein (SPEP_MPROT, g/L), urine UPEP (g/24h)

### ODE System (Srimani 2022)
```
dS/dt = k_R × (Y_SS − S) − [Imax × Cp / (IC50 + Cp)] × S
dR/dt = k_L × R

M-protein(t) ∝ S(t) + R(t)
```

### Parameters
| Parameter | Symbol | Units | Typical Value | IIV (CV%) |
|-----------|--------|-------|--------------|-----------|
| Inhibitory Emax | Imax | dimensionless | 0.758 | — |
| IC50 | IC50 | ng/mL | 3.29 | 42% |
| Sensitive cell turnover | k_R | /wk | 0.206 | 81% |
| Sensitive cell baseline | Y_SS | fraction | 0.143 | 155% |
| Resistant clone growth | k_L | /wk | 0.00951 | — |

### Key Properties
- Near-complete saturation at therapeutic dose (Cp ≈ 6–25× IC50) → flat E-R
- Resistant clone grows independently of Cp → late relapse is drug-independent
- AR(1) residuals: ε_t = 0.60·ε_{t-1} + √(1−0.36)·σ·z_t

---

## Model 2: FribergMyelosuppression
**Trials**: TOURMALINE (PLT), MONARCH-2 (ANC, CDK4/6-specific kinetics), EMILIA (PLT Day 8)
**Endpoint**: PLT (×10⁹/L), ANC (×10⁹/L)

### ODE System (Friberg 2002, modified for cumulative AUC)
```
dProl/dt  = k_prol × Prol × [1 − E(t)] − k_tr × Prol
dTrans1/dt = k_tr × Prol − k_tr × Trans1
dTrans2/dt = k_tr × Trans1 − k_tr × Trans2
dCirc/dt  = k_tr × Trans2 − k_circ × Circ

E(t) = β_inst × Cp(t) + β_cum × AUC_cum(t)   (for proteasome inhibitors)
     = β_slope × Cp(t)                          (for CDK4/6, one term)
```

### Parameters
| Parameter | Symbol | Typical Value (PLT, PI) | Notes |
|-----------|--------|------------------------|-------|
| Slope (instantaneous) | β_inst | from dip_amp calibration | Cp-driven |
| Slope (cumulative) | β_cum | 1.337e-5 /(ng·h/mL) | AUC_cum-driven |
| Transition rate | k_tr | 1.0–2.0 /day | Governs nadir timing |
| Nadir timing | Day 11–15 | PI drugs | Day 8: T-DM1 |
| Per-cycle IOV | ω_IOV | 0.03 | Mean-corrected log-normal |

### Drug-Specific Nadir Timing
- Proteasome inhibitors: Day 11–15 (MK maturation ~7–10 days)
- T-DM1: Day 8 (faster, DM1 payload acts on dividing MK precursors)
- CDK4/6: Day 14–21 (longer MTT due to G1 arrest kinetics)

---

## Model 3: TMDDModel
**Trials**: POLLUX (daratumumab, CD38), EMILIA (T-DM1 antibody component, HER2)
**Endpoint**: Free drug Cp (drives PD), total drug (measured in PC domain)

### ODE System (Gibiansky 2010, QSS approximation)
```
dCp_total/dt = −CL/V × Cp_free − kon × Cp_free × R_free + koff × RC
dR_free/dt   = ksyn_R − kdeg × R_free − kon × Cp_free × R_free + koff × RC
dRC/dt       = kon × Cp_free × R_free − koff × RC − kint × RC

Cp_free = QSS approximation (Gibiansky):
  Cp_free = [(Cp_total − R_total − Kd) + √((Cp_total − R_total − Kd)² + 4Kd·Cp_total)] / 2
```

### Parameters
| Parameter | Symbol | Notes |
|-----------|--------|-------|
| Association rate | kon | Drug-target binding on-rate |
| Dissociation rate | koff | koff/kon = Kd |
| Internalization rate | kint | RC complex internalized |
| Target synthesis | ksyn_R | Baseline target expression |
| Target degradation | kdeg_R | kdeg_R = ksyn_R / R_baseline |
| Target baseline | R0 | Patient covariate with IIV (CD38 expression varies) |

### Notes
- PD model (M-protein ODE) is driven by Cp_free, not Cp_total
- R0 as patient covariate: higher CD38 expression → faster drug binding → higher apparent CL
- NK cell depletion (daratumumab): CD38 also expressed on NK cells; secondary PD arm

---

## Model 4: BTKOccupancyALC
**Trial**: RESONATE (ibrutinib)
**Endpoints**: BTK occupancy (%), ALC (×10⁹/L), lymph node SPD

### ODE System
```
Irreversible BTK binding:
  d[BTK_free]/dt = ksyn_BTK − kdeg_BTK × BTK_free − kinact × Cp × BTK_free

BTK_occupied = BTK_total − BTK_free
BTK_occupancy(%) = 100 × BTK_occupied / BTK_total

ALC redistribution (two-compartment lymphocyte):
  d[ALC_blood]/dt  = k_exit × ALC_node − k_entry × ALC_blood − (ibrutinib effect on k_entry)
  d[ALC_node]/dt   = k_entry × ALC_blood − k_exit × ALC_node − k_kill × BTK_occ × ALC_node

Phase 1 (ibrutinib starts): ALC_blood rises (node → blood redistribution)
Phase 2 (weeks 8–24):       ALC_blood declines as CLL cells die
```

### Critical Implementation Note
ALC redistribution is a mechanistic signature that MUST be reproduced:
- Initial lymphocytosis: ALC rises 50–100% above baseline in first 2–8 weeks
- This is NOT progression — physicians continue ibrutinib through this phase
- Failure to model this makes the synthetic data clinically implausible

---

## Model 5: TumorGrowthInhibition
**Trials**: MONARCH-2 (abemaciclib), FLAURA (osimertinib), KEYNOTE-189 (pembrolizumab)
**Endpoint**: Sum of longest diameters (SLD, mm) — RECIST 1.1

### ODE System (Simeoni 2004)
```
dX1/dt = λ0 × X1 / (1 + (λ0/λ1 × (X1+X2+X3+X4))^ψ)^(1/ψ) − k1 × E(t) × X1
dX2/dt = k1 × E(t) × X1 − k2 × X2
dX3/dt = k2 × X2 − k2 × X3
dX4/dt = k2 × X3 − k2 × X4

E(t) = Emax × Cp(t) / (EC50 + Cp(t))   [standard]
     = delayed_Emax(t)                   [pembrolizumab — with Tlag]

SLD(t) = SLD_baseline × (X1+X2+X3+X4) / (X1+X2+X3+X4)_baseline
```

### RECIST 1.1 Classification
```
CR:  SLD = 0 (complete disappearance)
PR:  SLD decrease ≥30% from baseline
SD:  neither PR nor PD criteria met
PD:  SLD increase ≥20% AND ≥5mm from nadir
```

### Drug-Specific Modifications

**Osimertinib (FLAURA)**:
- Separate CNS compartment: E_CNS = Emax × Cp_CNS / (EC50_CNS + Cp_CNS)
  where EC50_CNS < EC50_primary (osimertinib has high CNS penetration)
- Resistance ODE: d[R_clone]/dt = k_res × AUC_cum × S(t)

**Pembrolizumab (KEYNOTE-189)**:
- T-cell activation delay: E(t) = E_max × (1 − exp(−kact × max(0, t − Tlag)))
- PD-L1 TPS modifies EC50: EC50_i = EC50_pop × exp(−β_PDL1 × TPS_i)
- Hyperprogression: ~10% of patients assigned k_growth × 2 in first 2 cycles

---

## Model 6: ADCPayloadRelease
**Trial**: EMILIA (T-DM1)
**Endpoints**: Tumor SLD (RECIST), PLT (Friberg), peripheral neuropathy

### Two-Level PK/PD
```
Level 1 — Antibody PK (plasma):
  Standard TMDD (see Model 3): total antibody → free antibody

Level 2 — Intracellular payload (tumor):
  d[DM1_tumor]/dt = kint × RC_tumor × N_payload/mol − kdeg_DM1 × DM1_tumor

where:
  RC_tumor      = antibody-receptor complex in tumor
  N_payload/mol = average DM1 molecules per antibody (~3.5 DAR)
  kint          = HER2 internalization rate
  kdeg_DM1      = intracellular DM1 degradation rate

Tumor kill: d[TumorCell]/dt = −(kgrowth − E_DM1(t)) × TumorCell
  E_DM1 = kmax × DM1_tumor / (KM + DM1_tumor)   [Michaelis-Menten in tumor]

Bystander effect: fraction of adjacent HER2-low cells killed by released DM1
  E_bystander = fbystander × E_DM1
```

---

## Model 7: BCL2OccupancyTLS
**Trial**: CLL14 (venetoclax)
**Endpoints**: ALC (×10⁹/L), lymph node SPD (iwCLL), TLS risk score

### ODE System
```
BCL-2 occupancy:
  Occ(t) = Cp(t)^n / (IC50^n + Cp(t)^n)   [Hill equation, n≈1 for venetoclax]

CLL cell apoptosis:
  d[CLL]/dt = kgrowth × CLL − kaptosis_base × CLL − (kaptosis_max × Occ(t)) × CLL

ALC(t) = CLL_blood(t) + Non-CLL lymphocytes (approximately constant)

TLS risk score:
  TLS_risk = CLL_burden × BCL2_Occ(t)   [threshold model for TLS events]
  TLS event probability at ramp-up step = logistic(TLS_risk − TLS_threshold)
```

### Ramp-Up Encoding
Venetoclax ramp-up MUST be mechanistically encoded:
- Week 1 (20mg): AUC low → low occupancy → low TLS risk → safe
- Each dose step: TLS risk re-evaluated before dose escalation
- ~5% of patients have ramp-up delay due to TLS risk

---

## Model 8: SyntheticLethalityER
**Trial**: SOLO-1 (olaparib)
**Endpoints**: PFS (time-to-event), CA-125 (secondary), ANC/HGB (toxicity)

### PD Model
```
PARP inhibition:
  PARP_occ(t) = Cp(t) / (IC50_PARP + Cp(t))   [competitive inhibition]

Synthetic lethality (BRCA-dependent):
  E_kill(t) = PARP_occ(t) × BRCA_sensitivity

where BRCA_sensitivity:
  BRCA1/2-mutated:     1.0   (full synthetic lethality)
  HRD+/BRCA-WT:        0.40  (partial)
  HRD−/BRCA-WT:        0.10  (minimal)

Tumor growth:
  d[Tumor]/dt = kgrowth × (1 − E_kill(t)) × Tumor − (natural apoptosis)

PFS:
  Event when Tumor(t) > Tumor_progression_threshold
  Time-to-progression derived from Tumor ODE, then Gaussian copula → survival
```

### BRCA Covariate
BRCA mutation status is the dominant biomarker:
- BRCA-mut: low IC50 (high sensitivity)
- HRD+/BRCA-WT: intermediate
- HRD−: high IC50, little benefit
IIV on IC50 within each BRCA category

---

## Model 9: TCellEngagementCRS (placeholder — future trial)
**Trial**: GO40516 or RATIFY-adjacent T-cell engager
**Endpoints**: CRS grade, B-cell depletion, tumor response (Lugano SPD)

```
T-cell bridging:
  Bridge(t) = kon_TCE × Cp(t) × T_cell(t) × Tumor_CD20(t)

T-cell activation:
  d[Tact]/dt = kact × Bridge(t) − kdec × Tact

Cytokine release:
  d[IL6]/dt = kIL6 × Tact × Tumor_cell − kelim_IL6 × IL6
  d[IFNg]/dt = kIFNg × Tact − kelim_IFNg × IFNg

CRS grade:
  Grade 0: IL6 < threshold1
  Grade 1: threshold1 ≤ IL6 < threshold2
  Grade 2: threshold2 ≤ IL6 < threshold3, no organ dysfunction
  Grade 3: IL6 ≥ threshold3 OR organ dysfunction

B-cell kill:
  d[B_cell]/dt = kgrowth_B − (kIL6 × CRS_activity) × B_cell

Step-up dosing mechanistically reduces C1D1 Cmax → lower Bridge(t) → lower CRS probability
```

---

## Model 10: FLT3KinaseInhibition
**Trial**: RATIFY (midostaurin)
**Endpoints**: CR rate, AML blast %, FLT3 allelic ratio dynamics

```
FLT3 kinase occupancy:
  Occ_FLT3(t) = Cp(t) / (IC50_FLT3 + Cp(t))   [competitive, reversible for midostaurin]

Blast kill:
  d[Blast]/dt = kprolif × (1 − Eff_chemo(t)) × Blast
                − kkill_FLT3 × Occ_FLT3(t) × Blast × FLT3_ITD_frac

FLT3 allelic ratio (AR) as covariate:
  AR = FLT3_ITD_alleles / (FLT3_WT_alleles + FLT3_ITD_alleles)
  High AR (≥0.5): more FLT3-dependent → higher response to midostaurin
  Low AR (<0.5):  less FLT3-dependent → smaller benefit

Chemotherapy effect:
  Eff_chemo(t) = Emax_chemo × 1{t ∈ [Day1, Day7]}   [7+3 induction window]

CR definition (AML-specific):
  CR when: Blast% < 5 AND ANC ≥ 1.0×10⁹/L AND PLT ≥ 100×10⁹/L
```

---

## Cross-Model Dependencies

| Downstream | Upstream | Mechanism | Direction |
|-----------|---------|-----------|-----------|
| Dose modification | PLT nadir | PLT < threshold → hold/reduce drug | ↓ dose |
| PK exposure | Dose | Lower dose → lower Cp → lower AUC | ↓ efficacy |
| M-protein response | PK (Cp, AUC) | AUC drives TwoPopulation ODE | ↑ AUC → ↑ suppression |
| PFS | M-protein Cycle 6 | Copula: deeper response → longer PFS | ↓ M-prot → ↑ PFS |
| TLS event | Venetoclax AUC + tumor burden | High burden × high Occ → TLS | Risk ∝ AUC×burden |
| CRS grade | C1D1 Cmax + T-cell:tumor ratio | Step-up dose reduces Cmax → ↓ CRS | ↓ Cmax → ↓ CRS |
| Resistance | Cumulative AUC | Acquired resistance rate ∝ AUC_cum | ↑ AUC → faster C797S |
