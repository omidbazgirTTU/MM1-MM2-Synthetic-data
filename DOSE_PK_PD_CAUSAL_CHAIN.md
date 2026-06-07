# Dose → PK → PD → Endpoint Causal Chain
## TOURMALINE-MM1/MM2 Ixazomib + Lenalidomide + Dexamethasone

*Sources: Gupta et al. 2017 (PK), Srimani et al. 2022 (PD), Friberg et al. (PLT model)*

---

## Causal Chain Diagram

```
                            ┌─────────────────────────────────────┐
                            │         DOSE (mg, oral)              │
                            │  Ixazomib 4mg weekly (Days 1,8,15)  │
                            └──────────────┬──────────────────────┘
                                           │
                            F=0.58, KA=0.5/h, ALAG=0.15h
                            4-compartment PBPK (Gupta 2017)
                            CL/F=1.86 L/h, V2/F=14.3 L
                            IIV: ω_CL=0.36, ω_KA=0.55
                            DDI (CYP3A4): CL_mult ~ [0.45–2.50]
                                           │
                                           ▼
                            ┌─────────────────────────────────────┐
                            │    Cp(t)  — plasma concentration     │
                            │   ng/mL, free drug (proteasome)      │
                            │   Cmax ≈ 41 ng/mL  (cycle 1)        │
                            │   AUC_inf ≈ 1247 ng·h/mL            │
                            │   t½ ≈ 9.5 days (tissue-bound)      │
                            └──────────┬───────────────┬──────────┘
                                       │               │
                   ┌───────────────────┘               └────────────────────┐
                   │  M-protein arm                       PLT arm            │
                   │  (Srimani 2022)                  (Srimani 2022 /       │
                   │                                   Friberg model)        │
                   ▼                                        ▼
    ┌──────────────────────────┐            ┌───────────────────────────────┐
    │  Two-population indirect │            │  Modified Friberg (no feedback)│
    │  response model          │            │                               │
    │                          │            │  E_IXA = slp_IXA·Cp(t)       │
    │  Sensitive cells:        │            │         + k_IXA·AUC_cum(t)   │
    │   dS/dt = k_R·(Y_SS−S)  │            │                               │
    │          − Imax·Cp/(     │            │  AUC_cum accumulates cycle    │
    │            IC₅₀+Cp)·S   │            │  by cycle → progressive PLT  │
    │                          │            │  nadir deepening             │
    │  Resistant cells:        │            └───────────────┬───────────────┘
    │   dR/dt = k_L·R          │                            │
    │                          │                            ▼
    │  M-protein(t) ∝ S(t)+R(t)│            ┌──────────────────────────────┐
    └──────────┬───────────────┘            │  PLT count (×10⁹/L)          │
               │                            │  Nadir cycle 1: ~100–120     │
               │                            │  Nadir deepens each cycle    │
               ▼                            └──────────────┬───────────────┘
    ┌──────────────────────┐                               │
    │  M-protein response  │             ┌─────────────────┘
    │  CR / VGPR / PR /    │             │  PLT < 75×10⁹/L → dose hold
    │  MR / SD / PD        │             │  PLT < 30×10⁹/L → dose reduce
    └──────────┬───────────┘             │  4mg→3mg→2.3mg steps
               │                         └──────────────┬───────────────────┐
               │                                        │                   │
               │                            ┌───────────┘           ┌───────┘
               │                            ▼                       ▼
               │              Cp_new(t) decreases         M-protein benefit
               │              (lower dose→lower AUC)      attenuated if
               │                            │              reduction sustained
               │                            └───────────────────────┐
               │                                                     │
               └──────────────────────┬──────────────────────────────┘
                                      ▼
                       ┌──────────────────────────────────┐
                       │   PFS (Progression-Free Survival) │
                       │                                   │
                       │   Driven by:                      │
                       │    • Rate of M-protein increase   │
                       │    • Resistant clone fraction     │
                       │    • Cumulative AUC (PLT toxicity)│
                       └──────────────────────────────────┘
```

---

## Quantitative Parameter Table

| Link | Parameter | Value | Source |
|------|-----------|-------|--------|
| **Dose → Cp** | Bioavailability F | 0.58 | Gupta 2017 |
| **Dose → Cp** | Absorption rate KA | 0.5 /h | Gupta 2017 |
| **Dose → Cp** | Absorption lag ALAG | 0.15 h | Gupta 2017 |
| **Dose → Cp** | CL/F (central) | 1.86 L/h | Gupta 2017 |
| **Dose → Cp** | V2/F (central) | 14.3 L | Gupta 2017 |
| **Dose → Cp** | Q3, V3 (peripheral 1) | 8.0 L/h, 200.0 L | Gupta 2017 |
| **Dose → Cp** | Q4, V4 (peripheral 2) | 0.50 L/h, 328.7 L | Gupta 2017 |
| **Dose → Cp** | IIV on CL (ω) | 0.36 (CV≈37%) | Gupta 2017 |
| **Dose → Cp** | IIV on KA (ω) | 0.55 (CV≈58%) | Gupta 2017 |
| **Dose → Cp** | DDI CL multiplier range | 0.45–2.50 | TOURMALINE ADSL |
| **Cp → M-protein** | IC₅₀ (sensitive cells) | 3.29 ng/mL | Srimani 2022 |
| **Cp → M-protein** | Imax | 0.758 | Srimani 2022 |
| **Cp → M-protein** | Sensitive cell turnover k_R | 0.206 /wk | Srimani 2022 |
| **Cp → M-protein** | Sensitive cell baseline Y_SS | 14.3% | Srimani 2022 |
| **Cp → M-protein** | Resistant cell growth k_L | 0.00951 /wk | Srimani 2022 |
| **IIV → M-protein** | IIV on Y_SS | 155% CV | Srimani 2022 |
| **IIV → M-protein** | IIV on k_R | 81% CV | Srimani 2022 |
| **IIV → M-protein** | IIV on IC₅₀ | 42% CV | Srimani 2022 |
| **Cp → PLT nadir** | Slope term (instantaneous) | slp_IXA (Srimani 2022 est.) | Srimani 2022 |
| **Cumulative AUC → PLT** | AUC term (progressive) | k_IXA (Srimani 2022 est.) | Srimani 2022 |
| **PLT → dose mod** | Hold threshold | PLT < 75 ×10⁹/L | TOURMALINE protocol |
| **PLT → dose mod** | Reduce threshold | PLT < 30 ×10⁹/L | TOURMALINE protocol |
| **Dose mod → Cp** | Dose reduction steps | 4mg → 3mg → 2.3mg | TOURMALINE protocol |
| **M-protein → PFS** | Progression definition | ≥25% M-protein increase from nadir | IMWG criteria |
| **Resistance → PFS** | Resistant clone expansion | Exponential k_L; unchecked by Cp | Srimani 2022 |

---

## Three Key Takeaways

### 1. PK → M-protein: Near-Complete Saturation at Therapeutic Dose

At the therapeutic dose of 4mg, peak plasma concentration Cp ≈ 20–80 ng/mL is **6–25× above IC₅₀ = 3.29 ng/mL**.

The Hill inhibition term becomes:

```
Imax · Cp / (IC₅₀ + Cp)  ≈  0.758 · 80 / (3.29 + 80)  ≈  0.73   (at Cmax)
```

This means the proteasome is **85–96% inhibited** throughout most of the dosing interval. As a result:

- A 2× increase in Cp (e.g., due to CYP3A4 inhibition) changes inhibition from ~92% to ~96% — a **4 percentage point difference**
- A 50% drop in Cp (e.g., due to dose reduction from 4mg→2mg) changes inhibition from ~92% to ~85% — a **7 percentage point difference**

**Practical implication**: The exposure-response (E-R) curve for M-protein is **flat at therapeutic exposure**. Classical E-R analysis (AUC vs ΔM-protein) will show weak or no correlation — not because the drug doesn't work, but because everyone is near-maximally inhibited. ML models trained to predict M-protein response from PK will face a fundamental ceiling effect.

---

### 2. PK → PLT: The Strongest Learnable PK-PD Signal

Unlike M-protein, platelet toxicity has **two distinct exposure drivers** that are both learnable:

**Instantaneous (Cp-driven):** Captures nadir depth in a given cycle.
```
E_instantaneous = slp_IXA · Cp(t)
```

**Cumulative (AUC_cum-driven):** Captures cycle-by-cycle deepening of PLT nadir.
```
E_cumulative = k_IXA · AUC_cumulative(t)
```

The cumulative term means that **a patient with stable Cp over cycles 1–6 will show progressively worsening thrombocytopenia** — the PLT nadir in cycle 6 is deeper than in cycle 1, even with the same weekly dose. This cycle-by-cycle signal is:

- **Causally linked to dose modifications** (holds and reductions)
- **Longitudinally structured** (deepens predictably over time)
- **Observable in structured data** (PLT counts recorded each cycle)

**Practical implication**: PK→PLT is the strongest learnable causal link in the TOURMALINE data. An ML model that predicts PLT trajectory from dose history and PK exposure (Cp, AUC_cum) has a real mechanistic signal to learn. This signal can in turn predict which patients will require dose modifications, and when.

---

### 3. Biology, Not PK, Dominates M-protein Variability

The dominant sources of between-patient variability in M-protein response are **biological**, not pharmacokinetic:

| Source | Variability | Driver |
|--------|------------|--------|
| Y_SS (sensitive cell baseline) | **155% CV** | Biology (tumor burden, clone composition) |
| k_R (sensitive cell turnover) | **81% CV** | Biology (cell cycle kinetics) |
| IC₅₀ (drug sensitivity) | **42% CV** | Biology (proteasome mutation status) |
| k_L (resistant clone growth) | Fixed in Srimani 2022 | Biology (clonal selection pressure) |
| CL/F (PK variability) | 37% CV | CYP3A4 genotype, DDI, body composition |
| KA (absorption) | 58% CV | GI motility, food effects |

Because Y_SS IIV (155%) >> CL IIV (37%), **most of the observed difference in M-protein response between patients is explained by their underlying tumor biology, not by differences in drug exposure**.

**Practical implication**: A model trying to predict M-protein response from PK alone (without tumor biology covariates) will have low R². Conversely, a model that captures Y_SS, k_R, and IC₅₀ at baseline (via M-protein kinetics in early cycles) will greatly outperform a PK-only model. In the TOURMALINE context, **early M-protein trajectory (cycles 1–2) is the best predictor of later response** — not the dose or the Cp.

---

## Full Biological Mechanism Summary

### Why Does Resistance Emerge Despite Near-Complete Proteasome Inhibition?

The two-population model reveals the mechanism:

1. **Sensitive cells** (S): Fully inhibited → suppressed; proliferation blocked; steady-state near zero if Cp >> IC₅₀
2. **Resistant cells** (R): Grow exponentially at rate k_L = 0.00951/wk, **independent of Cp**; not subject to Imax or IC₅₀

Even when proteasome inhibition is maximal, the resistant clone expands unimpeded. After approximately **6–12 months** (depending on initial R₀ fraction), the resistant clone dominates the M-protein signal, producing the characteristic **late relapse pattern** seen in TOURMALINE:

```
M-protein(t) = S(t) + R(t)
             ↓ early (drug suppresses S)   ↑ late (R grows unchecked)
                      \_____________________/
                           U-shaped trajectory
                           in late progressors
```

### Why Does Dose Reduction Hurt More Than It Appears?

When PLT toxicity forces a dose reduction (4mg → 3mg → 2.3mg):

1. Cp drops proportionally (linear PK, F·Dose/CL)
2. Inhibition fraction drops from ~92% to ~85% (a 7 pp change — small)
3. **But**: Sensitive cell killing rate drops by the same 7 pp
4. Sensitive cells now re-expand slightly (dS/dt becomes less negative)
5. Meanwhile, resistant clone continues growing at the same k_L
6. Net effect: **The M-protein nadir is less deep, and the time to progression may shorten**

This dose-response linkage is why patients who require dose modifications may have worse PFS outcomes — not because the drug stopped working, but because the ~7 pp reduction in inhibition shifts the balance toward re-expansion of the sensitive clone.

---

## References

1. **Gupta N et al. (2017)** — Population pharmacokinetic analysis of ixazomib, an oral proteasome inhibitor, in patients with hematologic malignancies. *Clin Pharmacokinet* 56(4):415–427. *(4-compartment PK model, DDI characterization, Gupta 2017 parameters)*

2. **Srimani JK et al. (2022)** — Computational modeling of ixazomib pharmacokinetics and pharmacodynamics to characterize drug-induced M-protein and platelet responses in multiple myeloma patients. *CPT Pharmacometrics Syst Pharmacol.* *(Two-population M-protein ODE + modified Friberg PLT model)*

3. **Friberg LE et al.** — Model of chemotherapy-induced myelosuppression with parameter consistency across drugs. *J Clin Oncol* 20(24):4713–4721. *(Original Friberg myelosuppression framework)*

4. **TOURMALINE-MM1 (Moreau et al. 2016)**, **TOURMALINE-MM2 (Richardson et al. 2023)** — Clinical trial design, dose modification rules, and published PK/efficacy reference values used for simulation validation.

---

*Document created: 2026-05-03*
*Project: TOURMALINE synthetic data generation and ML framework (Takeda-data)*
