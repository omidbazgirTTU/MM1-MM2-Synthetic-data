# Skills — TOURMALINE Project

Reusable procedural knowledge for the TOURMALINE MM1/MM2 synthetic data project.

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

## Skill Index

| Skill | Applies when | load_cost |
|-------|-------------|-----------|
| [tourmaline-data-generation-workflow.md](tourmaline-data-generation-workflow.md) | Adding domains, re-running generation, looking up CDISC conventions or output paths | low |
| [synthetic-clinical-trial-response-calibration.md](synthetic-clinical-trial-response-calibration.md) | ORR/VGPR/CR rates off target, CR appearing as VGPR, best-response waterfall issues | medium |
| [plt-grade3-myelosuppression-calibration.md](plt-grade3-myelosuppression-calibration.md) | Grade 3 PLT off target, dip_amp doesn't scale as expected after code change | medium |
| [pk-validation-cmax-vpc.md](pk-validation-cmax-vpc.md) | PK metrics failing, Cmax computation wrong, VPC construction | medium |
| [synthetic-data-rng-management.md](synthetic-data-rng-management.md) | Rate changed after unrelated edit, MM1 shifted after MM2 change, probability changes non-proportional | low |

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

---

## Current Project State (2026-04)

- **Validation**: 48/48 PASS
- **Seeds**: MM2=42, MM1=43, SURV_RNG=77
- **Generator**: `scripts/generate_v2.py` (all domains), `scripts/generate_pk_v2.py` (PK)
- **Next phase**: ML framework development (Neural ODE PK → PD → Endpoints)
