# TOURMALINE Trials - Comprehensive Synthetic Data Generation Specifications

**Version:** 1.0  
**Date:** April 13, 2026  
**Purpose:** Generate realistic synthetic datasets for PK→PD→Endpoint deep learning modeling

---

## Table of Contents
1. [Trial Overviews](#trial-overviews)
2. [Population Characteristics](#population-characteristics)
3. [Dosing Schedules](#dosing-schedules)
4. [Pharmacokinetic Parameters](#pharmacokinetic-parameters)
5. [Pharmacodynamic Biomarkers](#pharmacodynamic-biomarkers)
6. [Clinical Outcomes](#clinical-outcomes)
7. [Adverse Events](#adverse-events)
8. [Data Sampling Strategy](#data-sampling-strategy)
9. [Baseline Covariate Distributions](#baseline-covariate-distributions)
10. [Synthetic Data Generation Guide](#synthetic-data-generation-guide)

---

## 1. Trial Overviews

### TOURMALINE-MM1 (NCT01564537)
**Design:** Global, phase 3, randomized, double-blind, placebo-controlled  
**Population:** Relapsed and/or Refractory Multiple Myeloma  
**Sample Size:** N = 722 (360 Ixazomib-Rd, 362 Placebo-Rd)  
**Stratification Factors:**
- Number of prior therapies (1 vs 2 or 3)
- Previous proteasome inhibitor exposure (Yes vs No)
- ISS disease stage (I or II vs III)

**Primary Endpoint:** Progression-Free Survival (PFS) by independent review committee  
**Result:** ✅ PRIMARY ENDPOINT MET
- Median PFS: 20.6 months (Ixazomib-Rd) vs 14.7 months (Placebo-Rd)
- Hazard Ratio: 0.74 (95% CI: 0.59–0.94)
- P-value: 0.01

**Key Secondary Endpoints:**
- Overall Survival: Median 53.6 vs 51.6 months (HR 0.939, p=0.495)
- Overall Response Rate: Similar between arms
- High-risk cytogenetics: HR 0.870 for OS

**FDA Status:** APPROVED for RRMM (≥1 prior therapy)

### TOURMALINE-MM2 (NCT01850524)
**Design:** International, randomized, double-blind, placebo-controlled, phase 3  
**Population:** Newly Diagnosed Multiple Myeloma (Transplant-Ineligible)  
**Sample Size:** N = 705 (351 Ixazomib-Rd, 354 Placebo-Rd)

**Primary Endpoint:** Progression-Free Survival (PFS)  
**Result:** ❌ PRIMARY ENDPOINT NOT MET
- Median PFS: 35.3 months (Ixazomib-Rd) vs 21.8 months (Placebo-Rd)
- Hazard Ratio: 0.83 (95% CI: not reported)
- P-value: 0.073 (NOT statistically significant)
- Clinical benefit: 13.5 month improvement (clinically meaningful but not significant)

**Key Characteristics:**
- >40% of patients ≥75 years old
- ~40% high-risk cytogenetics
- 58% CrCl >60 ml/min

**FDA Status:** NOT APPROVED for NDMM

---

## 2. Population Characteristics

### A. TOURMALINE-MM1 (Relapsed/Refractory MM)

| Characteristic | Ixazomib-Rd (N=360) | Placebo-Rd (N=362) | Notes |
|---|---|---|---|
| **Demographics** ||||
| Median Age | ~66 years | ~66 years | Range: 23-91 |
| Age >65 years | ~50% | ~50% | |
| Male | ~57% | ~57% | |
| White race | ~80% | ~80% | |
| **Prior Therapies** ||||
| 1 prior line | 62% | 60% | |
| 2 prior lines | 27% | 31% | |
| 3 prior lines | 11% | 9% | |
| Prior PI exposure | 69% | 70% | Stratification factor |
| Prior IMiD exposure | ~55% | ~55% | |
| **ISS Stage** ||||
| Stage I | 63% | 64% | |
| Stage II | 25% | 24% | |
| Stage III | 12% | 12% | Stratification factor |
| **Disease Characteristics** ||||
| IgG myeloma | ~60% | ~60% | |
| IgA myeloma | ~25% | ~25% | |
| Light chain only | ~15% | ~15% | |
| High-risk cytogenetics‡ | ~20% | ~20% | del(17p), t(4;14), t(14;16) |
| del(17p) | ~10% | ~10% | >5% positive cells |
| t(4;14) | ~8% | ~8% | |
| t(14;16) | ~3% | ~3% | |
| +1q21 amplification | ~35% | ~35% | Expanded high-risk |
| **Renal Function** ||||
| CrCl ≤60 ml/min | ~30% | ~30% | Requires dose adjustment |
| CrCl >60 ml/min | ~70% | ~70% | |
| **Baseline Labs** ||||
| M-protein (mean ± SD) | 3.2 ± 2.1 g/dL | 3.1 ± 2.0 g/dL | SPEP |
| β2-microglobulin | Elevated | Elevated | Correlates with ISS |
| Platelet count (mean ± SD) | 220 ± 80 ×10⁹/L | 225 ± 82 ×10⁹/L | |
| Hemoglobin (mean ± SD) | 10.2 ± 1.8 g/dL | 10.1 ± 1.7 g/dL | |
| **ECOG Performance Status** ||||
| 0-1 | ~90% | ~90% | |
| 2 | ~10% | ~10% | |

‡ High-risk cytogenetics defined as presence of del(17p), t(4;14), and/or t(14;16)

### B. TOURMALINE-MM2 (Newly Diagnosed MM)

| Characteristic | Ixazomib-Rd (N=351) | Placebo-Rd (N=354) | Notes |
|---|---|---|---|
| **Demographics** ||||
| Median Age | ~73 years | ~73 years | Transplant-ineligible |
| Age ≥75 years | >40% | >40% | Older than MM1 |
| Male | ~55% | ~55% | |
| **Disease Characteristics** ||||
| High-risk cytogenetics | ~40% | ~40% | Higher than MM1 |
| ISS Stage I | ~35% | ~35% | Estimated |
| ISS Stage II | ~35% | ~35% | |
| ISS Stage III | ~30% | ~30% | |
| **Renal Function** ||||
| CrCl >60 ml/min | 58% | 58% | Lower than MM1 |
| CrCl ≤60 ml/min | 42% | 42% | More renal impairment |
| **Baseline Labs** ||||
| M-protein (mean ± SD) | 3.8 ± 2.3 g/dL | 3.7 ± 2.2 g/dL | Higher tumor burden |
| Platelet count (mean ± SD) | 200 ± 85 ×10⁹/L | 205 ± 83 ×10⁹/L | |
| Hemoglobin (mean ± SD) | 9.8 ± 1.9 g/dL | 9.9 ± 1.8 g/dL | More anemic |

**Key Differences MM1 vs MM2:**
- **Age:** MM2 older (median 73 vs 66 years)
- **Prior therapy:** MM1 has 1-3 prior lines; MM2 is treatment-naive
- **Renal function:** MM2 worse (42% CrCl ≤60 vs 30%)
- **Tumor burden:** MM2 higher baseline M-protein
- **High-risk cytogenetics:** MM2 higher prevalence (40% vs 20%)

---

## 3. Dosing Schedules

### Standard Regimen (Both Trials)
**28-Day Cycles, Treatment Until Disease Progression or Unacceptable Toxicity**

#### Ixazomib
- **Dose:** 4.0 mg oral capsule
- **Schedule:** Days 1, 8, 15 of each 28-day cycle
- **Administration:** Taken at approximately the same time on dosing days
- **Food:** At least 1 hour before or at least 2 hours after food
- **Dose Modifications:**
  - Start at 3 mg if moderate/severe hepatic impairment
  - Start at 3 mg if severe renal impairment (CrCl <30) or ESRD on dialysis
  - Reduce to 3 mg for Grade 3 thrombocytopenia or Grade 2/3 PN
  - Hold for Grade 4 thrombocytopenia or Grade 4 PN
  - Resume at 2.3 mg after recovery

#### Lenalidomide
- **Dose:** 25 mg oral capsule
- **Schedule:** Days 1-21 of each 28-day cycle
- **Dose Adjustments:**
  - 10 mg if CrCl ≤60 ml/min (based on local guidelines)
  - 10 mg if CrCl ≤50 ml/min (stricter sites)
  - Reduce for neutropenia/thrombocytopenia per label

#### Dexamethasone
- **Dose:** 40 mg oral
- **Schedule:** Days 1, 8, 15, 22 of each 28-day cycle
- **Dose Adjustments:**
  - 20 mg for patients >75 years (common practice)
  - Reduce for severe AEs

### Placebo Arm
- **Placebo for ixazomib:** Days 1, 8, 15
- **Lenalidomide:** 25 mg Days 1-21
- **Dexamethasone:** 40 mg Days 1, 8, 15, 22

### Treatment Duration
- **MM1 (RRMM):**
  - Median cycles Ixazomib-Rd: 18 cycles (range 1-70+)
  - Median cycles Placebo-Rd: 16 cycles
- **MM2 (NDMM):**
  - Longer treatment expected due to newly diagnosed setting
  - Median cycles not fully reported

---

## 4. Pharmacokinetic Parameters

### Published Two-Compartment Model (3-Compartment with Peripheral)

**Structural Model:**
- **Model Type:** Three-compartment with first-order absorption and linear elimination
- **Absorption:** Rapid (Ka = 0.5 h⁻¹)
- **Bioavailability (F):** 58-60% (oral)
- **Tmax:** 0.5-1.0 hour (median)
- **Dose proportionality:** Linear (dose-independent PK)

**Population PK Parameters:**

| Parameter | Population Estimate | Unit | Notes |
|---|---|---|---|
| **Clearance (CL)** | 2.0 L/h | L/h | Alternatively reported as 44.2 L/h in some models |
| **Central Volume (V2)** | 543 L | L | Steady-state Vd |
| **Peripheral Volume (V3)** | Not specified | L | First peripheral compartment |
| **Peripheral Volume (V4)** | ~1560 L | L | Second peripheral compartment |
| **Absorption Rate (Ka)** | 0.5 h⁻¹ | h⁻¹ | Rapid absorption |
| **Bioavailability (F)** | 0.58 | - | 58% |
| **Inter-compartmental Clearance (Q3)** | Not specified | L/h | |
| **Inter-compartmental Clearance (Q4)** | 17.1 L/h | L/h | |

**Covariate Effects:**
- **Body Surface Area (BSA):** Small effect on V4 (12.9% reduction in BSV)
  - No effect on CL
  - Fixed dosing (4 mg) justified
- **Creatinine Clearance:** No significant effect on CL
  - CrCl range tested: 22-213.7 ml/min
  - No dose adjustment needed for mild/moderate renal impairment
- **Age:** No significant effect
  - Age range: 23-86 years
- **Weight/BSA:** Minimal effect
  - BSA range: 1.3-2.6 m²

**Inter-Individual Variability (IIV):**
- CL: ~30-40% CV
- V2: ~25-35% CV
- V4: ~45-58% CV (reduced to 45.2% with BSA covariate)

**Residual Variability:**
- Proportional error: ~20-30% CV

**Exposure Metrics (4 mg dose):**
- **Cmax (geometric mean):** 89.1 ng/mL (CV 62.3%)
- **AUC0-312h:** 1180 ng·h/mL (CV 46.0%)
- **Renal clearance:** 0.119 L/h
- **Urinary recovery:** 3.23% (2.13% SD) of dose as unchanged drug over 168h

**Elimination:**
- **Routes:** Primarily hepatic metabolism (non-CYP mediated)
- **Urine:** 62% (mostly metabolites)
- **Feces:** 22%
- **Half-life:** Multi-compartment (distribution + elimination)

**PK Sampling Strategy for Synthetic Data (4 measurements per dosing interval):**

For Days 1, 8, 15 of each cycle (ixazomib dosing days):

| Time Point | Hours Post-Dose | Purpose | Expected Concentration (ng/mL) |
|---|---|---|---|
| **Pre-dose** | 0 | Trough from prior dose (Day 8, 15 only) | <5 (mostly eliminated) |
| **Peak** | 1.0 | Cmax | 70-110 (CV ~60%) |
| **Mid-interval** | 4.0 | Post-distribution | 30-50 |
| **Late** | 8.0 | Pre-trough | 10-20 |

**Note:** Ixazomib has minimal accumulation due to rapid clearance. Day 1 has no pre-dose sample (first dose).

---

## 5. Pharmacodynamic Biomarkers

### A. M-Protein (Primary PD Marker)

**Measurement:** Serum Protein Electrophoresis (SPEP)
- **Frequency:** Day 1 of each cycle (every 4 weeks)
- **Assay:** Electrophoresis with immunofixation
- **Units:** g/dL

**Baseline Values:**

| Trial | Treatment | Mean ± SD (g/dL) | Range |
|---|---|---|---|
| MM1 | Ixazomib-Rd | 3.2 ± 2.1 | 0.5-9.0 |
| MM1 | Placebo-Rd | 3.1 ± 2.0 | 0.5-8.8 |
| MM2 | Ixazomib-Rd | 3.8 ± 2.3 | 0.8-10.0 |
| MM2 | Placebo-Rd | 3.7 ± 2.2 | 0.8-9.5 |

**Expected Trajectories (Mean values by cycle):**

**MM1 Ixazomib-Rd (Responders):**
- Cycle 1: 3.2 g/dL (baseline)
- Cycle 2: 2.4 g/dL (-25%)
- Cycle 3: 1.8 g/dL (-44%)
- Cycle 4: 1.4 g/dL (-56% - PR achieved)
- Cycle 6: 0.9 g/dL (-72%)
- Cycle 9: 0.4 g/dL (-88%)
- Cycle 12: 0.2 g/dL (-94% - VGPR)
- Cycle 18+: 0.0-0.1 g/dL (CR in ~15-20%)

**MM1 Placebo-Rd (Responders):**
- Slower decline, shallower depth
- Cycle 4: 1.8 g/dL (-42%)
- Cycle 12: 0.8 g/dL (-74%)
- CR rate: ~10%

**MM2 Ixazomib-Rd (Responders):**
- Similar pattern but starting from higher baseline
- Deeper responses due to treatment-naive status
- VGPR/CR rates higher (~25-30% CR)

**Response Categories (IMWG Criteria):**
- **CR (Complete Response):** M-protein 0 (undetectable)
- **VGPR (Very Good PR):** ≥90% reduction
- **PR (Partial Response):** ≥50% reduction
- **SD (Stable Disease):** <50% reduction, <25% increase
- **PD (Progressive Disease):** ≥25% increase OR new lesions

**MM1 Best Response Rates:**
- ORR (≥PR): ~78% (Ixazomib-Rd) vs ~72% (Placebo-Rd)
- VGPR/CR: ~48% vs ~39%
- CR: ~12% vs ~7%

### B. Platelet Dynamics (Dose-Limiting Toxicity)

**Measurement:** Automated hematology analyzer
- **Frequency:**
  - Cycles 1-6: Days 1, 8, 15 of each cycle
  - Cycles 7+: Day 1 (reduced frequency)
- **Units:** ×10⁹/L

**Baseline Values:**

| Trial | Treatment | Mean ± SD (×10⁹/L) |
|---|---|---|
| MM1 | Ixazomib-Rd | 220 ± 80 |
| MM1 | Placebo-Rd | 225 ± 82 |
| MM2 | Ixazomib-Rd | 200 ± 85 |
| MM2 | Placebo-Rd | 205 ± 83 |

**Kinetic Pattern (Ixazomib-Rd):**

Within each 28-day cycle:
- **Day 1 (pre-dose):** Baseline or recovered from prior cycle
  - Cycle 1: ~220 ×10⁹/L
  - Subsequent cycles: ~180-200 (slight cumulative effect)
- **Day 8 (pre-dose):** Beginning of decline
  - ~180-190 ×10⁹/L
- **Day 11-15:** NADIR (predictable timing)
  - ~110-130 ×10⁹/L (mean nadir)
  - Grade 3 (<50): ~25-30% of patients
  - Grade 4 (<25): ~5-8% of patients
- **Day 21:** Recovery begins
  - ~140-160 ×10⁹/L
- **Day 28:** Near-complete recovery
  - ~170-190 ×10⁹/L

**Placebo-Rd Pattern:**
- Less severe nadirs (~140-160 ×10⁹/L)
- Grade 3 thrombocytopenia: ~14-16% vs 25-31% with ixazomib

**Dose Modifications:**
- Hold dose if platelets <30 ×10⁹/L
- Reduce dose if platelets 30-75 ×10⁹/L
- Resume at lower dose after recovery

### C. Hemoglobin

**Measurement:** Automated hematology
- **Frequency:** Days 1, 8, 15 (Cycles 1-6), then Day 1
- **Units:** g/dL

**Baseline:**

| Trial | Treatment | Mean ± SD (g/dL) |
|---|---|---|
| MM1 | Ixazomib-Rd | 10.2 ± 1.8 |
| MM1 | Placebo-Rd | 10.1 ± 1.7 |
| MM2 | Ixazomib-Rd | 9.8 ± 1.9 |
| MM2 | Placebo-Rd | 9.9 ± 1.8 |

**Trajectories:**
- **Responders:** Gradual increase over 6-12 months
  - By Cycle 6: 11.0-11.5 g/dL
  - By Cycle 12: 12.0-12.5 g/dL (normalization)
- **Non-responders:** Stable or declining
  - Risk of transfusion dependence
- **Grade 3 anemia (Hgb <8):** ~25-30% of patients

### D. Serum Free Light Chains

**Measurement:** Serum FLC assay (Freelite or equivalent)
- **Frequency:** Every 4 weeks
- **Units:** mg/L

**Normal Range:**
- Kappa: 3.3-19.4 mg/L
- Lambda: 5.7-26.3 mg/L
- Ratio (κ/λ): 0.26-1.65

**Baseline (for patients with measurable FLC disease, ~15%):**

**Kappa myeloma:**
- Kappa: 500-15,000 mg/L (highly elevated)
- Lambda: 10-50 mg/L (suppressed)
- Ratio: 20-1000 (abnormal)

**Lambda myeloma:**
- Kappa: 5-30 mg/L (suppressed)
- Lambda: 400-10,000 mg/L (elevated)
- Ratio: 0.01-0.10 (abnormal)

**Response Pattern (Kappa myeloma example):**
- Baseline ratio: 50
- Cycle 3: 20 (improvement)
- Cycle 6: 8
- Cycle 12: 2.5
- CR: 0.5-1.5 (normalization)

### E. Other Biomarkers

**β2-Microglobulin:**
- Baseline: Elevated (correlates with ISS stage)
  - ISS I: <3.5 mg/L
  - ISS II: 3.5-5.5 mg/L
  - ISS III: ≥5.5 mg/L
- Response: Decreases with tumor reduction

**LDH:**
- Baseline: Often elevated
- High LDH = aggressive disease, poor prognosis

**Bone Marrow Plasma Cells:**
- Measured at baseline and select timepoints
- Baseline: 30-90% plasma cells
- CR requires <5%

---

## 6. Clinical Outcomes

### A. Progression-Free Survival (PFS)

**MM1 Results:**

| Outcome | Ixazomib-Rd | Placebo-Rd | HR (95% CI) | P-value |
|---|---|---|---|---|
| Median PFS | 20.6 months | 14.7 months | 0.74 (0.59-0.94) | 0.01 |
| 1-year PFS | ~68% | ~57% | | |
| 2-year PFS | ~42% | ~28% | | |

**PFS by Subgroup (Ixazomib-Rd favored in all):**
- 1 prior therapy: HR 0.70
- 2-3 prior therapies: HR 0.83
- Prior PI exposure: HR 0.74
- No prior PI: HR 0.74
- High-risk cytogenetics: HR 0.54
- Standard risk: HR 0.76

**MM2 Results:**

| Outcome | Ixazomib-Rd | Placebo-Rd | HR (95% CI) | P-value |
|---|---|---|---|---|
| Median PFS | 35.3 months | 21.8 months | 0.83 (~0.68-1.02) | 0.073 NS |
| Clinical benefit | +13.5 months | | | |

### B. Overall Survival (OS)

**MM1 Final Analysis (85 months median follow-up):**

| Outcome | Ixazomib-Rd | Placebo-Rd | HR (95% CI) | P-value |
|---|---|---|---|---|
| Median OS | 53.6 months | 51.6 months | 0.939 (0.784-1.124) | 0.495 NS |
| 5-year OS | ~40% | ~38% | | |

**OS Subgroups with Benefit (HR <1):**
- Refractory to any prior line: HR 0.794
- Refractory to last line: HR 0.742
- Age >65-75 years: HR 0.757
- ISS Stage III: HR 0.779
- 2-3 prior therapies: HR 0.845
- High-risk cytogenetics: HR 0.870

**MM2 OS:**
- Not yet mature at primary analysis

### C. Response Rates

**MM1 Best Response (IRC-assessed):**

| Response | Ixazomib-Rd | Placebo-Rd |
|---|---|---|
| ORR (≥PR) | 78% | 72% |
| VGPR + CR | 48% | 39% |
| CR | 12% | 7% |
| VGPR | 36% | 32% |
| PR | 30% | 33% |

**Time to Response:**
- Median time to first response: ~2-3 months (Cycle 2-3)
- Median time to best response: ~6-9 months

**Duration of Response:**
- Median DoR: ~21 months (Ixazomib-Rd) vs ~16 months (Placebo-Rd)

---

## 7. Adverse Events

### A. Common Adverse Events (Any Grade)

**MM1 Data:**

| Adverse Event | Ixazomib-Rd (N=360) | Placebo-Rd (N=362) |
|---|---|---|
| **Hematologic** | | |
| Thrombocytopenia† | 85% | 67% |
| Neutropenia† | 74% | 70% |
| Anemia† | 75% | 73% |
| **Gastrointestinal** | | |
| Diarrhea | 52% | 43% |
| Nausea | 29% | 22% |
| Vomiting | 23% | 12% |
| Constipation | 34% | 26% |
| **Neurologic** | | |
| Peripheral neuropathy | 27-28% | 21-22% |
| **Other** | | |
| Peripheral edema | 28% | 20% |
| Back pain | 23% | 17% |
| Fatigue | 30% | 25% |
| Rash | 35% | 21% |

† Pooled from adverse events and laboratory data

### B. Grade 3-4 Adverse Events

**MM1 Severe AEs:**

| Adverse Event | Ixazomib-Rd Grade 3/4 | Placebo-Rd Grade 3/4 |
|---|---|---|
| **Any Grade ≥3 AE** | 74% | 69% |
| **Serious AEs** | 47% | 49% |
| **Hematologic** | | |
| Thrombocytopenia | 31% / 5% | 16% / 2% |
| Neutropenia | 33% / 12% | 31% / 9% |
| Anemia | 29% / 4% | 27% / 3% |
| **Non-Hematologic** | | |
| Pneumonia | 8% | 6% |
| Diarrhea | 6% | 4% |
| Peripheral neuropathy | 2% | 2% |
| Peripheral edema | 2% | 1% |
| Acute renal failure | 8% | 10% |
| Heart failure | 4% | 3% |

### C. Discontinuation and Dose Modifications

**Treatment Discontinuation:**
- Due to AEs (Ixazomib-Rd): 4%
- Due to AEs (Placebo-Rd): 6%

**Dose Reductions:**
- Ixazomib dose reduced: ~20-25% of patients
- Most common reasons: Thrombocytopenia, GI events, PN

**Death on Study:**
- Similar rates between arms (~5-6%)
- Most deaths due to disease progression

### D. Chinese Population Data (Higher Exposures)

From China Continuation Study:
- **Systemic exposures:** 80% higher than White patients
- **Thrombocytopenia Grade 3/4:** 18%/7% (ixazomib) vs 14%/5% (placebo)
- **Neutropenia Grade 3/4:** 19%/5% vs 19%/2%
- **Overall AE profile:** Similar despite higher exposures

---

## 8. Data Sampling Strategy

### Overview
Generate synthetic datasets with temporal resolution sufficient for PK→PD→Endpoint modeling while matching clinical trial assessment schedules.

### A. PK Sampling (4 measurements per dosing interval)

**Ixazomib Dosing Days: 1, 8, 15 of each cycle**

For EACH dosing day, collect 4 PK samples:

| Sample | Time Post-Dose | Purpose | Notes |
|---|---|---|---|
| Pre-dose | 0 hours | Trough (Days 8, 15 only) | Day 1 Cycle 1: N/A |
| Peak | 1 hour | Cmax | Median Tmax = 0.5-1h |
| Mid | 4 hours | Post-distribution | Early elimination phase |
| Late | 8 hours | Pre-trough | Late elimination |

**Total PK samples per cycle:** 3 dosing days × 4 samples = 12 samples/cycle  
(11 samples for Cycle 1, since Day 1 has no pre-dose)

**Additional Sparse Sampling (Optional):**
- Day 2: Single trough (24h post-dose)
- Day 4: Single sample (to capture elimination)

### B. PD Biomarker Sampling

**Aligned with Schedule of Assessments:**

**Cycle 1 (Intensive Monitoring):**
- Days 1, 8, 15, 22: Full panel
  - M-protein (SPEP)
  - Platelets, Hgb, ANC
  - Chemistry (CrCl, LFTs)
  - FLC (if applicable)

**Cycle 2-3:**
- Days 1, 8, 15: Full panel
- Day 22: Not required (unless clinically indicated)

**Cycle 4+ (Maintenance):**
- Day 1: Full panel every 4 weeks
- Days 8, 15: Safety labs only (platelets, Hgb) if clinically stable

**Response Assessment:**
- Every 4 weeks (Day 1 of each cycle)
- M-protein measurement
- Clinical assessment for PD

### C. Outcome Data

**Progression-Free Survival:**
- **Event:** Disease progression (IMWG criteria) or death
- **Censoring:** Last assessment date if no event
- **Assessment frequency:** Every 4 weeks

**Generate for each patient:**
- PFS time (days or months)
- PFS event indicator (0 = censored, 1 = progressed/died)
- Date of progression (if event occurred)
- Reason (PD vs death)

**Overall Survival:**
- **Event:** Death from any cause
- **Follow-up:** Every 12 weeks after discontinuation
- **Median follow-up:** 85 months for MM1

**Response:**
- Best response achieved (CR, VGPR, PR, SD, PD)
- Time to first response
- Time to best response
- Duration of response (for responders)

**Adverse Events:**
- For each AE type:
  - Onset day/cycle
  - Maximum grade (1-5)
  - Resolution day (if resolved)
  - Action taken (dose reduction, hold, discontinue)
  - Relationship to study drug

---

## 9. Baseline Covariate Distributions

### A. Continuous Variables

| Variable | MM1 Mean (SD) | MM1 Range | MM2 Mean (SD) | MM2 Range | Distribution |
|---|---|---|---|---|---|
| Age (years) | 66 (10) | 23-91 | 73 (7) | 47-93 | Normal (truncated) |
| Weight (kg) | 75 (18) | 40-140 | 72 (17) | 40-130 | Log-normal |
| BSA (m²) | 1.85 (0.25) | 1.3-2.6 | 1.82 (0.24) | 1.3-2.5 | Normal |
| CrCl (ml/min) | 75 (32) | 22-214 | 65 (30) | 25-180 | Log-normal |
| M-protein (g/dL) | 3.2 (2.1) | 0.5-9.0 | 3.8 (2.3) | 0.8-10.0 | Gamma |
| Platelets (×10⁹/L) | 222 (81) | 50-500 | 203 (84) | 45-480 | Log-normal |
| Hemoglobin (g/dL) | 10.15 (1.75) | 6.5-15.0 | 9.85 (1.85) | 6.0-14.0 | Normal |
| β2-microglobulin (mg/L) | 4.5 (3.2) | 1.5-25.0 | 5.2 (3.8) | 1.8-28.0 | Gamma |
| LDH (U/L) | 190 (85) | 100-800 | 210 (95) | 100-850 | Log-normal |

### B. Categorical Variables

**Prior Therapies (MM1 only):**
- 1 prior: 61% (Binomial)
- 2 prior: 29%
- 3 prior: 10%

**ISS Stage:**
- MM1: I (63%), II (25%), III (12%)
- MM2: I (35%), II (35%), III (30%) [estimated]
- Use multinomial distribution

**Cytogenetics (Binary for each):**
- del(17p): 10% (MM1), 15% (MM2)
- t(4;14): 8% (MM1), 12% (MM2)
- t(14;16): 3% (MM1), 5% (MM2)
- +1q21: 35% (MM1), 42% (MM2)
- High-risk (any of del17p, t(4;14), t(14;16)): 20% (MM1), 40% (MM2)

**Myeloma Type:**
- IgG: 60%
- IgA: 25%
- Light chain only: 15%

**Prior PI Exposure (MM1):**
- Yes: 70%
- No: 30%

**Prior IMiD Exposure (MM1):**
- Yes: 55%
- No: 45%

**ECOG Performance Status:**
- 0-1: 90%
- 2: 10%

**Sex:**
- Male: 57% (MM1), 55% (MM2)
- Female: 43% (MM1), 45% (MM2)

**Race:**
- White: 80%
- Asian: 10%
- Black: 5%
- Other: 5%

### C. Covariate Correlations

**Important correlations to preserve:**

1. **Age and Renal Function:**
   - Older age → Lower CrCl
   - Correlation: r = -0.45

2. **ISS Stage and Biomarkers:**
   - ISS stage determined by β2M and albumin
   - Higher ISS → Higher β2M, lower albumin

3. **M-protein and Disease Burden:**
   - Higher M-protein → More anemia
   - Correlation: r = -0.30

4. **Cytogenetics and Outcomes:**
   - High-risk cytogenetics → Worse PFS/OS
   - Multiplicative effect on hazard

5. **Treatment Arm Balance:**
   - Stratification ensures balance on:
     - Number of prior therapies
     - Prior PI exposure
     - ISS stage
   - Other covariates should be similar between arms

---

## 10. Synthetic Data Generation Guide

### Step-by-Step Process

#### Step 1: Generate Baseline Cohort

**For MM1 (N=722):**

```python
# Treatment assignment (1:1 randomization)
treatment_arm = random.choice(['Ixazomib-Rd', 'Placebo-Rd'], n=722)

# Stratification factors (balanced across arms)
num_prior_therapies = stratified_sample([1, 2, 3], p=[0.61, 0.29, 0.10])
prior_pi_exposure = stratified_sample([True, False], p=[0.70, 0.30])
iss_stage = stratified_sample([1, 2, 3], p=[0.63, 0.25, 0.12])

# Demographics
age = truncated_normal(mean=66, sd=10, min=23, max=91)
sex = random.choice(['M', 'F'], p=[0.57, 0.43])
race = random.choice(['White', 'Asian', 'Black', 'Other'], p=[0.80, 0.10, 0.05, 0.05])

# Anthropometrics
weight = lognormal(mean=75, sd=18)
bsa = calculate_bsa(weight, height)

# Renal function (correlated with age)
crcl = lognormal(mean=75 - 0.5*(age-66), sd=32)

# Disease characteristics
m_protein = gamma(shape=2.3, scale=1.4)  # Mean 3.2, SD 2.1
platelets = lognormal(mean=222, sd=81)
hemoglobin = normal(mean=10.15, sd=1.75)
b2m = gamma(shape=1.9, scale=2.4)
ldh = lognormal(mean=190, sd=85)

# Cytogenetics
del17p = random.choice([True, False], p=[0.10, 0.90])
t4_14 = random.choice([True, False], p=[0.08, 0.92])
t14_16 = random.choice([True, False], p=[0.03, 0.97])
amp1q21 = random.choice([True, False], p=[0.35, 0.65])
high_risk = del17p OR t4_14 OR t14_16

# Myeloma type
mm_type = random.choice(['IgG', 'IgA', 'LC'], p=[0.60, 0.25, 0.15])
```

**For MM2 (N=705):**
- Adjust distributions for older, frailer, treatment-naive population
- Higher age (mean 73 vs 66)
- Worse renal function (42% CrCl ≤60)
- Higher high-risk cytogenetics (40% vs 20%)
- NO prior therapy variables

#### Step 2: Generate PK Parameters

For each patient, predict individual PK parameters:

```python
# Population PK parameters
CL_pop = 2.0  # L/h
V2_pop = 543  # L
V4_pop = 1560  # L
Q4_pop = 17.1  # L/h
Ka = 0.5  # h^-1
F = 0.58

# Individual parameters (with IIV)
eta_CL = random.normal(0, 0.35)  # 35% CV
eta_V2 = random.normal(0, 0.30)  # 30% CV
eta_V4 = random.normal(0, 0.45)  # 45% CV (reduced to 0.45 with BSA covariate)

CL_ind = CL_pop * exp(eta_CL)
V2_ind = V2_pop * exp(eta_V2)
V4_ind = V4_pop * exp(eta_V4 + 0.13 * (bsa - 1.85))  # BSA effect

# Store for PK simulation
patient_pk_params = {
    'CL': CL_ind,
    'V2': V2_ind,
    'V4': V4_ind,
    'Q4': Q4_pop,
    'Ka': Ka,
    'F': F
}
```

#### Step 3: Simulate Dosing and PK

```python
# For treatment arm patients only
if treatment_arm == 'Ixazomib-Rd':
    for cycle in range(1, max_cycles+1):
        # Dosing days: 1, 8, 15
        dose_days = [1, 8, 15]
        
        for day in dose_days:
            dose_mg = 4.0  # Standard dose
            
            # Apply dose modifications based on AEs
            if platelet_nadir_prior_cycle < 30:
                dose_mg = 0  # Hold dose
            elif platelet_nadir_prior_cycle < 75:
                dose_mg = 3.0  # Reduce dose
            
            # Simulate 4 PK samples
            if day > 1 or cycle > 1:
                C_predose = simulate_concentration(time=0, patient_pk_params)
            
            C_peak = simulate_concentration(time=1, patient_pk_params, dose=dose_mg)
            C_mid = simulate_concentration(time=4, patient_pk_params, dose=dose_mg)
            C_late = simulate_concentration(time=8, patient_pk_params, dose=dose_mg)
            
            # Add measurement error (~20-30% CV)
            C_peak_measured = C_peak * lognormal(1, 0.25)
            # ... repeat for other timepoints
```

#### Step 4: Simulate PD Biomarkers

```python
# M-protein dynamics
for cycle in range(1, max_cycles+1):
    # Calculate drug effect on tumor cell kill
    AUC_cycle = calculate_auc_cycle(pk_concentrations)
    
    # Patient-specific sensitivity
    IC50 = lognormal(mean=50, sd=25)  # ng/mL, varies by cytogenetics
    if high_risk:
        IC50 *= 1.5  # Less sensitive
    
    k_kill = k_kill_max * (AUC_cycle / (IC50 + AUC_cycle))
    
    # M-protein turnover
    M_protein_new = M_protein_prev * exp(-k_kill * 28) + k_prod * tumor_burden
    
    # Add measurement noise
    M_protein_measured = M_protein_new * lognormal(1, 0.10)
    
# Platelet dynamics
for day in [1, 8, 15] of each cycle:
    # Megakaryocyte inhibition from ixazomib
    inhibition = E_max * (C_ixazomib / (EC50_plt + C_ixazomib))
    
    # Platelet production with delay
    MK_production = MK_production_baseline * (1 - inhibition)
    
    # Nadir typically Day 11-15
    if day == 15:
        platelet_count = platelet_baseline * exp(-cumulative_inhibition)
    
    # Recovery by Day 28
    if day == 1 and cycle > 1:
        platelet_count = platelet_prev_nadir * recovery_factor
```

#### Step 5: Generate Clinical Outcomes

```python
# PFS simulation (parametric survival model)

# Baseline hazard depends on:
hazard_base = 0.06  # ~16 months median for placebo-Rd

# Covariate effects on hazard
HR_treatment = 0.74 if arm == 'Ixazomib-Rd' else 1.0
HR_iss3 = 1.5 if iss_stage == 3 else 1.0
HR_high_risk = 1.4 if high_risk else 1.0
HR_age = 1.02 ** (age - 66)  # Per year over 66
HR_prior = 1.15 ** (num_prior_therapies - 1)

# Time-dependent effect from PD
HR_m_protein_slope = exp(-0.5 * m_protein_decline_rate_week_0_12)

# Total hazard
hazard_total = hazard_base * HR_treatment * HR_iss3 * HR_high_risk * HR_age * HR_prior * HR_m_protein_slope

# Simulate time to event
PFS_time = exponential(rate=hazard_total)

# Censor at data cutoff or dropout
censored = PFS_time > followup_time OR dropout_occurred

# Best response
response = determine_response(m_protein_trajectory, flc_trajectory, bone_marrow)
```

#### Step 6: Generate Adverse Events

```python
# Thrombocytopenia (linked to platelet PK/PD model)
for cycle in range(1, max_cycles+1):
    plt_nadir = platelet_dynamics[cycle]['nadir']
    
    if plt_nadir < 50:
        thrombocytopenia_grade = 3
        ae_record = {
            'term': 'Thrombocytopenia',
            'grade': 3,
            'onset_cycle': cycle,
            'onset_day': 15,  # Typical nadir
            'action': 'dose_reduction' if plt_nadir > 30 else 'dose_hold'
        }
    elif plt_nadir < 75:
        thrombocytopenia_grade = 2
        # ... record Grade 2

# Other AEs (probabilistic, not mechanistically linked)
# Diarrhea
if arm == 'Ixazomib-Rd':
    diarrhea_prob = 0.52
    diarrhea_grade3_prob = 0.06
else:
    diarrhea_prob = 0.43
    diarrhea_grade3_prob = 0.04

if random.uniform() < diarrhea_prob:
    grade = 3 if random.uniform() < (diarrhea_grade3_prob / diarrhea_prob) else 1
    # ... record AE

# Peripheral neuropathy (cumulative risk)
pn_risk_cycle = 0.28 / median_cycles  # ~28% overall
pn_cumulative = 1 - (1 - pn_risk_cycle) ** cycle

if random.uniform() < pn_cumulative AND not pn_already_occurred:
    grade = 2 if random.uniform() < 0.92 else 3  # Mostly Grade 1-2
    # ... record AE
```

### Quality Control Checks

After generating synthetic data, verify:

1. **Population balance:**
   - Treatment arms balanced on stratification factors
   - Covariate distributions match published baselines

2. **PK/PD consistency:**
   - PK exposures match published AUC, Cmax
   - Platelet nadirs occur at expected timing
   - M-protein declines align with response rates

3. **Outcome concordance:**
   - Median PFS matches trial results (MM1: 20.6 vs 14.7 months)
   - OS trends reasonable (MM1: 53.6 vs 51.6 months)
   - Response rates match (ORR ~78% vs 72%)

4. **AE rates:**
   - Grade 3 thrombocytopenia: ~31% (ixazomib) vs ~16% (placebo)
   - Discontinuation rates: ~4% vs ~6%

5. **Correlations preserved:**
   - Age vs CrCl (r = -0.45)
   - M-protein slope vs PFS (strong predictor)
   - High-risk cytogenetics vs worse outcomes

---

## Appendices

### Appendix A: Data Dictionary

Full variable specifications for synthetic data generation.

**Patient Identifiers:**
- `patient_id`: Unique identifier
- `trial`: "MM1" or "MM2"
- `site_id`: Randomized site number
- `treatment_arm`: "Ixazomib-Rd" or "Placebo-Rd"

**Demographics:**
- `age`: Years (continuous)
- `sex`: "M" or "F"
- `race`: "White", "Asian", "Black", "Other"
- `weight_kg`: Kilograms
- `bsa_m2`: Body surface area (m²)

**Disease Characteristics:**
- `mm_type`: "IgG", "IgA", "LC"
- `iss_stage`: 1, 2, 3
- `ecog_ps`: 0, 1, 2
- `num_prior_therapies`: 1, 2, 3 (MM1 only)
- `prior_pi`: TRUE/FALSE (MM1 only)
- `prior_imid`: TRUE/FALSE (MM1 only)

**Cytogenetics:**
- `del17p`: TRUE/FALSE
- `t4_14`: TRUE/FALSE
- `t14_16`: TRUE/FALSE
- `amp1q21`: TRUE/FALSE
- `high_risk_cyto`: TRUE/FALSE

**Baseline Labs:**
- `m_protein_baseline`: g/dL
- `plt_baseline`: ×10⁹/L
- `hgb_baseline`: g/dL
- `anc_baseline`: ×10⁹/L
- `crcl_baseline`: ml/min
- `b2m_baseline`: mg/L
- `ldh_baseline`: U/L
- `albumin_baseline`: g/dL

**PK Data (per sample):**
- `patient_id`, `cycle`, `day`, `time_post_dose`
- `ixazomib_concentration_ng_ml`
- `sample_type`: "predose", "peak", "mid", "late"

**PD Data (per timepoint):**
- `patient_id`, `cycle`, `day`
- `m_protein_g_dl`
- `plt_count`
- `hgb_g_dl`
- `anc`
- `flc_kappa_mg_l`, `flc_lambda_mg_l`, `flc_ratio`

**Dosing Data:**
- `patient_id`, `cycle`, `day`, `drug`
- `planned_dose_mg`
- `actual_dose_mg`
- `dose_taken`: TRUE/FALSE
- `reason_not_taken`: If applicable

**Outcomes:**
- `pfs_days`: From randomization
- `pfs_event`: 0 (censored) or 1 (event)
- `os_days`: From randomization
- `os_event`: 0 (censored) or 1 (death)
- `best_response`: "CR", "VGPR", "PR", "SD", "PD"
- `response_date`: Date of best response
- `time_to_response_days`
- `duration_of_response_days`

**Adverse Events:**
- `patient_id`, `ae_id`, `ae_term`
- `grade`: 1-5
- `onset_cycle`, `onset_day`
- `resolution_cycle`, `resolution_day`
- `serious`: TRUE/FALSE
- `related_to_drug`: "ixazomib", "lenalidomide", "dex", "multiple", "none"
- `action_taken`: "none", "dose_reduction", "dose_hold", "discontinuation"

### Appendix B: Software Recommendations

**Python Libraries:**
- `numpy`, `scipy`: Statistical distributions
- `pandas`: Data manipulation
- `PyTorch`: PK/PD ODE solving, neural network models
- `lifelines`: Survival analysis
- `matplotlib`, `seaborn`: Visualization

**R Packages:**
- `simsurv`: Survival data simulation
- `MASS`: Multivariate distributions
- `survival`: Cox models
- `mrgsolve`: PK/PD simulation

### Appendix C: References

1. Moreau P, et al. Oral Ixazomib, Lenalidomide, and Dexamethasone for Multiple Myeloma. N Engl J Med. 2016;374:1621-34.
2. Richardson PG, et al. Final Overall Survival Analysis of TOURMALINE-MM1. J Clin Oncol. 2021;39(24):2430-2442.
3. Facon T, et al. TOURMALINE-MM2: Ixazomib with Lenalidomide-Dexamethasone in Newly Diagnosed MM. Blood. 2020.
4. Gupta N, et al. Population Pharmacokinetic Analysis of Ixazomib. Clin Pharmacokinet. 2017;56(9):1025-1037.
5. Kumar SK, et al. Pharmacodynamic modeling of ixazomib from TOURMALINE-MM1. CPT Pharmacometrics Syst Pharmacol. 2022.

---

**END OF SPECIFICATIONS**

**For questions or clarifications on synthetic data generation, refer to:**
- Schedule of Assessments documents (TOURMALINE-MM1 and MM2)
- PK_PD_Mathematical_Framework.md
- ClinicalTrials.gov: NCT01564537, NCT01850524
