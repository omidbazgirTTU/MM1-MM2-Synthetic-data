---
name: orchestrator
description: Coordination protocol for the multi-agent synthetic clinical trial generation system. Defines round structure, agent sequencing, end condition, and failure routing.
metadata:
  type: rules
---

# Orchestrator Protocol

## Purpose
The Orchestrator drives one trial config from a bare NCT ID or trial name to a validated,
agent-approved synthetic dataset. It does not do domain reasoning — it sequences agents,
routes failures, and enforces the end condition.

---

## Inputs
- Trial identifier (trial name, NCT ID, or branch name from [[trial-registry]])
- Target: validated synthetic SDTM + ADaM dataset, 68+ criteria PASS, all 4 agents signed off

## Outputs
- Frozen YAML trial config in `trials/<name>/config.yaml`
- Generated data in `trials/<name>/MM*/` (SDTM + ADaM domains)
- Validation report in `trials/<name>/outputs/VALIDATION_REPORT.md`
- Agent sign-off log in `trials/<name>/outputs/agent_signoff.md`

---

## Round Structure

### Round 0 — Branch Setup
1. Confirm the trial branch exists (see [[trial-registry]])
2. Check out the trial branch
3. Confirm `trials/<name>/config.yaml` stub exists; if not, create from schema template
4. Initialize agent sign-off log with all four agents in `PENDING` state

### Round 1 — Literature Grounding
Invoke **Literature Search Agent** with: trial name, NCT ID, primary drug(s), indication.

Expected outputs (structured evidence package):
- Published popPK parameter table (source, year, journal)
- Trial design (N per arm, arms, primary endpoint, dose, schedule)
- Published summary statistics to serve as validation targets
- Any published popPK/PD ODE model with parameters

Block until evidence package is returned. If key parameters are missing, Literature Search
re-queries with more specific terms before proceeding.

### Round 2 — Parallel Domain Review
Invoke all three domain agents simultaneously with the evidence package:

| Agent | Primary task |
|-------|-------------|
| **Oncologist** | Confirm patient population, disease biology, endpoint definitions |
| **QSP Scientist** | Propose PK/PD model structure, parameterize ODEs, identify IIV/covariate effects |
| **Drug Developer** | Confirm trial design, dosing schedule, dose modification rules, arm structure |

Each agent returns a **domain proposal** (structured dict). Conflicts between proposals
are flagged for Round 3.

### Round 3 — Cross-Agent Challenge
Route conflicts identified in Round 2:
- QSP ↔ Oncologist conflicts → QSP and Oncologist exchange directly
- QSP ↔ Drug Developer conflicts → QSP and Drug Developer exchange directly
- Any unresolved conflict → Literature Search re-queries for additional evidence

Typical challenges to expect:
- "Baseline PLT of 200 is too high for heavily pretreated RRMM — use 160"
- "The 20 mg/m² dose requires BSA normalization, not flat dosing"
- "Published IIV on CL is from healthy volunteers, inflate by 30% for cancer patients"
- "This IC50 estimate is from in vitro — adjust for protein binding in vivo"

Iterate Round 3 until no unresolved conflicts remain (max 3 iterations).

### Round 4 — Config Convergence
QSP Scientist assembles the final YAML config from:
- Literature Search evidence package (parameters + validation targets)
- Oncologist population spec
- Drug Developer trial design spec
- Any parameter adjustments from Round 3

Config written to `trials/<name>/config.yaml`. All agents acknowledge the config.

### Round 5 — Generate + Validate (Automated)
```bash
python engine/run_trial.py trials/<name>/config.yaml
python engine/validate_trial.py trials/<name>/config.yaml
```

Parse validation output. If any criterion FAILs:
- Route to QSP Scientist for parameter recalibration
- QSP proposes minimal adjustment (one parameter at a time)
- Re-run generation and validation
- Max 10 recalibration iterations before escalating to human review

### Round 6 — Agent Sign-Off
Each agent independently reviews the validation output and generated data:

**Literature Search** sign-off criteria:
- All simulated summary statistics within published 95% CI or ±15% tolerance
- PK NCA parameters (Cmax, AUCinf, t½) within ±20% of published values
- Response rates within ±15% of published trial primary and secondary endpoints
- Grade ≥3 toxicity rates within ±20% of published safety table

**QSP Scientist** sign-off criteria:
- PK/PD mechanistic dependencies reproduced (correct direction and magnitude)
- IIV CV% within ±10% of published popPK model estimates
- Cross-correlation structure matches published physiological relationships
- Trajectory shapes (M-protein kinetics, PLT nadir timing, etc.) mechanistically consistent
- VPC: ≥80% of observations within 5th–95th prediction interval

**Oncologist** sign-off criteria:
- Patient demographics match published Table 1 (within ±15%)
- Baseline disease characteristics clinically plausible for the indication
- Response trajectory shapes consistent with clinical experience
- Toxicity patterns consistent with the drug's known safety profile
- Response assessment timing matches trial schedule of assessments

**Drug Developer** sign-off criteria:
- Arm structure and N per arm correctly reproduced
- Dosing schedule encoded accurately (days, cycles, route)
- Dose modification rules match protocol (holds, reductions, thresholds)
- Concomitant medication rates clinically plausible
- Exposure separation between arms confirmed

---

## End Condition
The loop closes **only when all four agents return APPROVED**. A single PENDING or REJECTED
from any agent keeps the loop open.

```
Literature Search: APPROVED ✓
QSP Scientist:     APPROVED ✓
Oncologist:        APPROVED ✓
Drug Developer:    APPROVED ✓
→ TRIAL COMPLETE — config frozen, data committed, branch ready for merge
```

## Failure Routing
| Failure type | Routed to |
|-------------|-----------|
| Validation criterion FAIL (quantitative) | QSP Scientist → recalibrate |
| Mechanistic trajectory wrong | QSP Scientist → revisit ODE parameters |
| Patient characteristics off | Oncologist → adjust population spec |
| Dosing schedule error | Drug Developer → correct config |
| Published target not found | Literature Search → re-query |
| Parameter conflict (QSP vs Lit) | Literature Search + QSP → arbitrate |

## What the Orchestrator Never Does
- Does not make pharmacometric judgments (that is QSP's role)
- Does not override agent sign-offs
- Does not modify the config without routing through an agent
- Does not commit data until all four agents APPROVED
