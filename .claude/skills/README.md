# Skills, Agents, and Trials — Multi-Agent Synthetic Clinical Trial Framework

Reusable procedural knowledge and agent definitions for the 11-trial synthetic
oncology data generation framework. MM1/MM2 (Takeda/Ixazomib) are the validated
seed pair. Nine additional Phase 3 trials cover distinct drug mechanisms.

## Loading Levels

Each skill file has **3 levels** of content, loaded progressively as needed:

| Level | Content type | When to load |
|-------|-------------|-------------|
| **1 — Metadata** | YAML frontmatter: `name`, `description`, `applies_when`, `keywords`, `load_cost` | Always: scan to decide relevance before loading body |
| **2 — Instructions** | Main body: workflows, decision guides, calibration rules, best practices | When the skill matches the current task |
| **3 — Resources & Code** | Code blocks, parameter tables, exact file:line references | When actively implementing or debugging |

**Pattern**: Read frontmatter of all skills first (cheap). Load the full body of only
the matching skill(s). Pull Level 3 code only when writing or modifying actual code.

---

## Agent Definitions (`.claude/agents/`)

| Agent | Role | File |
|-------|------|------|
| Orchestrator | Sequences agents, enforces end condition, routes failures | [agents/orchestrator.md](../agents/orchestrator.md) |
| QSP Scientist | PK/PD model selection, ODE parameterization, recalibration | [agents/qsp_scientist.md](../agents/qsp_scientist.md) |
| Oncologist | Patient population, disease biology, clinical plausibility | [agents/oncologist.md](../agents/oncologist.md) |
| Drug Developer | Trial design, dosing schedule, dose modifications, CDISC encoding | [agents/drug_developer.md](../agents/drug_developer.md) |
| Literature Search | Published PK parameters, trial stats, evidence package | [agents/literature_search.md](../agents/literature_search.md) |

## Trial Registry (`.claude/trials/`)

| File | Contents |
|------|---------|
| [trials/registry.md](../trials/registry.md) | All 11 trials: branch, status, mechanism, N, published sources |

## Skills Index

| Skill | Applies when | load_cost |
|-------|-------------|-----------|
| [pd_modeling_guide.md](pd_modeling_guide.md) | Selecting/implementing a PD ODE class; wiring PK→PD dependencies | medium |
| [trial_config_schema.md](trial_config_schema.md) | Writing or validating a trial config.yaml | medium |
| [validation_framework.md](validation_framework.md) | Defining validation criteria, tolerances, agent sign-off checklists | medium |
| [tourmaline-data-generation-workflow.md](tourmaline-data-generation-workflow.md) | Adding domains, re-running generation, CDISC conventions | low |
| [synthetic-clinical-trial-response-calibration.md](synthetic-clinical-trial-response-calibration.md) | ORR/VGPR/CR rates off target | medium |
| [plt-grade3-myelosuppression-calibration.md](plt-grade3-myelosuppression-calibration.md) | Grade 3 PLT off target, dip_amp calibration | medium |
| [pk-validation-cmax-vpc.md](pk-validation-cmax-vpc.md) | PK metrics failing, Cmax wrong, VPC construction | medium |
| [synthetic-data-rng-management.md](synthetic-data-rng-management.md) | Rate drifted after unrelated edit, RNG isolation | low |
| [Cross_Correlations_Synthetic_Data_Guide.md](Cross_Correlations_Synthetic_Data_Guide.md) | MVN baseline correlations, OMEGA Cholesky, IXAZ_CL_I | medium |
| [Published_Mechanistic_Correlations_TOURMALINE.md](Published_Mechanistic_Correlations_TOURMALINE.md) | Srimani 2022 AUC→PLT, M-protein two-population ODE, flat E-R | medium |

---

## Quick Diagnostics

| Symptom | Load this skill |
|---------|----------------|
| VGPR or CR below target | `synthetic-clinical-trial-response-calibration` |
| CR low but VGPR OK | `synthetic-clinical-trial-response-calibration` → CR tier boundary section |
| PLT Grade 3 off target | `plt-grade3-myelosuppression-calibration` |
| Rate drifted after unrelated code edit | `synthetic-data-rng-management` |
| MM1 rates shifted when MM2 code changed | `synthetic-data-rng-management` → per-study reseeding |
| Cmax validation failing | `pk-validation-cmax-vpc` → check if it's the known +30% issue first |
| Adding a new SDTM/ADaM domain | `tourmaline-data-generation-workflow` |
| Probability change produces non-proportional rate change | `synthetic-data-rng-management` → RNG.choice vs RNG.random |
| Need to implement baseline covariate correlations | `Cross_Correlations_Synthetic_Data_Guide` |
| Age↔CrCl or Weight↔BSA correlation wrong | `Cross_Correlations_Synthetic_Data_Guide` → MVN section |
| Deciding linear vs Emax for AUC→PLT model | `Published_Mechanistic_Correlations_TOURMALINE` → Section 1 |
| Linking M-protein slope to PFS hazard | `Published_Mechanistic_Correlations_TOURMALINE` → Section 4 |
| AUC→PFS exposure-response should be flat? | `Published_Mechanistic_Correlations_TOURMALINE` → Section 3 |
| Criteria 49–68 failing in validate_data.py | `Cross_Correlations_Synthetic_Data_Guide` + `Published_Mechanistic_Correlations_TOURMALINE` |

---

## Current Project State (2026-06-08)

- **Validation**: 68/68 PASS (MM1 and MM2, branch `MM1-MM2`)
- **Seeds**: MM2=42, MM1=43, SURV_RNG=77 — frozen, never change
- **Generator**: `scripts/generate_v2.py` (all domains), `scripts/generate_pk_v2.py` (PK)
- **Next phase**: Option A — generalize engine into `engine/` library; write trial configs
- **Pending trials**: POLLUX, CLL14, RESONATE, MONARCH-2, EMILIA, FLAURA, SOLO-1, KEYNOTE-189, RATIFY
- **Architecture docs**: `MECHANISTIC_MODEL_AND_CROSS_CORRELATIONS.md`, `DOSE_PK_PD_CAUSAL_CHAIN.md`
