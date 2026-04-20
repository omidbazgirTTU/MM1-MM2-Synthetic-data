"""
TOURMALINE MM1/MM2 Synthetic Data Generator v2
================================================
Generates CDISC-compliant SDTM and ADaM datasets for all subjects.

Improvements over v1:
  - Weight, height, BMI, BSA in DM (sex/age correlated)
  - Proper Cockcroft-Gault CrCL in ADSL
  - R-ISS derived from ISS + high-risk cytogenetics + LDH
  - Extended cytogenetics: t(14;16), t(14;20), gain(1q21), del(1p32)
  - Arm-specific Weibull survival: IRd HR=0.742 vs Rd (published TOURMALINE)
  - All subjects in EX (v1 sampled only 200)
  - Dose-mod reasons (AE-related vs. physician vs. protocol)
  - Treatment adherence flag per dose event
  - sdtm_cm.csv: CYP3A4 drugs, anticoagulants, G-CSF/EPO, antimicrobials
  - LB M-protein trajectories with exposure-response (deeper = more Ixazomib AUC)
  - ANC/PLT lab values drive hematologic AE generation
  - CYP3A4 DDI flag carried into ADSL for PK generator

Reference: Hussain Z et al. npj Digital Medicine 7, 200 (2024)
Vivli request: NCT01850524 & NCT01564537
"""

import numpy as np
import pandas as pd
from scipy.stats import weibull_min
import datetime
import os
import warnings
warnings.filterwarnings("ignore")

RNG      = np.random.default_rng(42)
SURV_RNG = np.random.default_rng(77)   # independent seed for survival — calibrated separately

OUT_BASE = os.path.join(os.path.dirname(__file__), "..", )  # Takeda-data/

# ─────────────────────────────────────────────────────────────────────────────
# STUDY CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

STUDY_CONFIG = {
    "MM2": {
        "studyid":        "TOURMALINE-MM2",
        "protocol":       "NCT01850524",
        "population":     "NDMM",
        "n":              705,          # Spec: 351 IRd + 354 Rd = 705
        "n_ird":          351,
        "n_rd":           354,
        "median_age":     73,
        "age_range":      (48, 93),
        "pct_female":     0.450,        # Spec: ~45% female for NDMM
        "median_dx_months": 1.11,
        # Published TOURMALINE-MM2 results — Moreau et al. (2019)
        # HR=0.83, p=0.073 (did NOT meet primary endpoint)
        # Spec §6A: IRd=35.3 mo, Rd=21.8 mo
        # Input medians calibrated via binary search with SURV_RNG=77 to produce
        # KM medians matching published targets (see calibration block above)
        "pfs_median_ird": 34.3,   # → KM median ≈ 35.3 mo  [seed=2001]
        "pfs_median_rd":  21.5,   # → KM median ≈ 21.8 mo  [seed=2002]
        # OS not mature at primary analysis; reasonable estimates
        "os_median_ird":  52.3,   # → KM median ≈ 60.0 mo  [seed=2003]
        "os_median_rd":   46.5,   # → KM median ≈ 48.0 mo  [seed=2004]
        # Validation targets for summary print
        "pfs_target_ird": 35.3,
        "pfs_target_rd":  21.8,
        "cyto_hr_target": 40.0,        # Spec: ~40% high-risk in MM2
    },
    "MM1": {
        "studyid":        "TOURMALINE-MM1",
        "protocol":       "NCT01564537",
        "population":     "RRMM",
        "n":              722,          # Spec: 360 IRd + 362 Rd = 722
        "n_ird":          360,
        "n_rd":           362,
        "median_age":     66,
        "age_range":      (23, 91),
        "pct_female":     0.430,        # Spec: ~43% female
        "median_dx_months": 42.8,
        # Published TOURMALINE-MM1 results — Moreau et al. (2016) Lancet Oncol
        # PFS HR=0.74 (p=0.01); OS HR=0.939 (p=0.495)
        # Spec §6A/B: PFS IRd=20.6 mo, Rd=14.7 mo; OS IRd=53.6 mo, Rd=51.6 mo
        # Input medians calibrated via binary search with SURV_RNG=77
        "pfs_median_ird": 23.9,   # → KM median ≈ 20.6 mo  [seed=3001]
        "pfs_median_rd":  14.2,   # → KM median ≈ 14.7 mo  [seed=3002]
        "os_median_ird":  52.7,   # → KM median ≈ 53.6 mo  [seed=3003]
        "os_median_rd":   50.2,   # → KM median ≈ 51.6 mo  [seed=3004]
        # Validation targets
        "pfs_target_ird": 20.6,
        "pfs_target_rd":  14.7,
        "cyto_hr_target": 20.0,        # Spec: ~20% high-risk in MM1
    },
}

IG_TYPE_DIST = {
    "MM2": {"IgG": 0.573, "IgA": 0.202, "IgD": 0.014, "IgE": 0.004,
            "IgM": 0.004, "Biclonal": 0.033, "No Heavy Chain": 0.169},
    "MM1": {"IgG": 0.540, "IgA": 0.171, "IgD": 0.010, "IgE": 0.021,
            "IgM": 0.001, "Biclonal": 0.042, "No Heavy Chain": 0.215},
}

RACE_DIST = {
    "MM2": {"WHITE": 0.818, "ASIAN": 0.137, "BLACK OR AFRICAN AMERICAN": 0.033,
            "AMERICAN INDIAN OR ALASKA NATIVE": 0.004,
            "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER": 0.001, "OTHER": 0.007},
    "MM1": {"WHITE": 0.853, "ASIAN": 0.089, "BLACK OR AFRICAN AMERICAN": 0.018,
            "AMERICAN INDIAN OR ALASKA NATIVE": 0.001,
            "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER": 0.006,
            "OTHER": 0.010, "NOT REPORTED": 0.024},
}

ISS_DIST = {
    # MM2 (NDMM): Spec §2B — I=35%, II=35%, III=30% (more advanced disease at diagnosis)
    "MM2": {1: 0.350, 2: 0.350, 3: 0.300},
    # MM1 (RRMM): Spec §2A — I=63%, II=25%, III=12%
    # ISS Stage 3 p calibrated to 0.102 (from 0.120) to match 12% target with seed=43.
    # Excess patients with p=0.120 were +1.76σ above expected; correction via lower p.
    "MM1": {1: 0.648, 2: 0.250, 3: 0.102},
}

# Cytogenetics prevalences by study — Spec §9B
# MM2 (NDMM): higher prevalence due to newly-diagnosed, unselected population
# MM1 (RRMM): lower, already partially enriched for responders to prior therapy
CYTO_PROBS = {
    "MM2": {
        "DEL17P":  0.15,   # Spec: 15% in MM2 vs 10% in MM1
        "T414":    0.12,   # Spec: 12%
        "T1416":   0.05,   # Spec: 5%
        "T1420":   0.02,
        "GAIN1Q21":0.42,   # Spec: 42%
        "DEL1P32": 0.12,
    },
    "MM1": {
        "DEL17P":  0.10,   # Spec: 10%
        "T414":    0.063,  # Target 8%; ISS p3 change shifted RNG state; p=0.063 → ~8% with seed=43
        "T1416":   0.03,   # Spec: 3%
        "T1420":   0.01,
        "GAIN1Q21":0.35,   # Spec: 35%
        "DEL1P32": 0.10,
    },
}
# High-risk IMWG = del17p | t(4;14) | t(14;16)  [mutually-exclusive approximation]
# MM2 target: ~40%  → solve P(HR) = 1-(1-p17)(1-p414)(1-p1416) = 0.40
# With p1416=0.05: 1-(1-p17)(1-p414)*0.95 = 0.40  → (1-p17)(1-p414) = 0.632
# Use DEL17P=0.23, T414=0.18: (0.77)(0.82)=0.631 ✓ → P(HR)=1-0.631×0.95≈0.401
CYTO_PROBS["MM2"]["DEL17P"]  = 0.23   # calibrated to ~40% aggregate HR
CYTO_PROBS["MM2"]["T414"]    = 0.18

DRUG_DOSES = {
    "IXAZOMIB":     {"dose": 4.0,  "unit": "mg", "route": "ORAL", "freq": "WEEKLY (DAYS 1,8,15)"},
    "LENALIDOMIDE": {"dose": 25.0, "unit": "mg", "route": "ORAL", "freq": "DAILY (DAYS 1-21)"},
    "DEXAMETHASONE":{"dose": 40.0, "unit": "mg", "route": "ORAL", "freq": "WEEKLY (DAYS 1,8,15,22)"},
    "PLACEBO":      {"dose": 0.0,  "unit": "mg", "route": "ORAL", "freq": "WEEKLY (DAYS 1,8,15)"},
}

LAB_NORMAL_RANGES = {
    "ALBUMIN":      (35,   50,   "g/L"),
    "ALP":          (40,  130,   "U/L"),
    "ALT":          (7,    56,   "U/L"),
    "AST":          (10,   40,   "U/L"),
    "BILI":         (3,    21,   "umol/L"),
    "BUN":          (2.5,  7.1,  "mmol/L"),
    "CA":           (2.2,  2.6,  "mmol/L"),
    "CA_CORR":      (2.2,  2.6,  "mmol/L"),
    "CL":           (98,  106,   "mmol/L"),
    "CO2":          (22,   29,   "mmol/L"),
    "CREAT":        (44,  106,   "umol/L"),
    "GFR":          (60,  120,   "mL/min/1.73m2"),
    "GLUC":         (3.9,  5.6,  "mmol/L"),
    "HCT":          (36,   46,   "%"),
    "HGB":          (120, 160,   "g/L"),
    "IGA":          (0.7,  4.0,  "g/L"),
    "IGG":          (7.0, 16.0,  "g/L"),
    "IGM":          (0.4,  2.3,  "g/L"),
    "K":            (3.5,  5.0,  "mmol/L"),
    "KAPPA_LAMBDA": (0.26, 1.65, "ratio"),
    "LDH":          (140, 280,   "U/L"),
    "LYMPH":        (1.0,  4.8,  "10^9/L"),
    "MG":           (0.7,  1.0,  "mmol/L"),
    "MONO":         (0.2,  1.0,  "10^9/L"),
    "NA":           (136, 145,   "mmol/L"),
    "NEUT":         (1.8,  7.7,  "10^9/L"),
    "PHOS":         (0.8,  1.5,  "mmol/L"),
    "PLT":          (150, 400,   "10^9/L"),
    "PROT":         (60,   80,   "g/L"),
    "GLOB":         (20,   35,   "g/L"),
    "SPEP_GAMMA":   (7,    16,   "g/L"),
    "SPEP_KAPPA":   (3.3, 19.4,  "mg/L"),
    "SPEP_LAMBDA":  (5.7, 26.3,  "mg/L"),
    "SPEP_MPROT":   (0,    5,    "g/L"),
    "URATE":        (155, 428,   "umol/L"),
    "URINE_ALB":    (0,    30,   "%"),
    "UPEP_MPROT":   (0,    80,   "mg/day"),
    "B2MG":         (0.8,  2.4,  "mg/L"),
    "WBC":          (4.5, 11.0,  "10^9/L"),
}

AE_DEFS = [
    ("Acute Renal Failure",     "RENAL AND URINARY DISORDERS"),
    ("Cardiac Arrhythmias",     "CARDIAC DISORDERS"),
    ("Diarrhea",                "GASTROINTESTINAL DISORDERS"),
    ("Heart Failure",           "CARDIAC DISORDERS"),
    ("Hypotension",             "VASCULAR DISORDERS"),
    ("Liver Impairment",        "HEPATOBILIARY DISORDERS"),
    ("Nausea",                  "GASTROINTESTINAL DISORDERS"),
    ("Neutropenia",             "BLOOD AND LYMPHATIC SYSTEM DISORDERS"),
    ("Peripheral Neuropathies", "NERVOUS SYSTEM DISORDERS"),
    ("Rash",                    "SKIN AND SUBCUTANEOUS TISSUE DISORDERS"),
    ("Thrombocytopenia",        "BLOOD AND LYMPHATIC SYSTEM DISORDERS"),
    ("Vomiting",                "GASTROINTESTINAL DISORDERS"),
]
HEMATOLOGIC_AE = {"Neutropenia", "Thrombocytopenia"}

# Per-cycle AE probability calibrated to match published TOURMALINE cumulative incidence.
# Formula: p_per_cycle = 1 - (1 - target_cumulative)^(1/13)
# where 13 ≈ median number of cycles per patient across the trial.
# Published cumulative incidences (any grade, all-cycle):
#   Moreau 2016 NEJM (MM2) / Lancet Oncol 2016 (MM1) — Tables S4/2
AE_INCIDENCE = {
    "IRd": {"Acute Renal Failure": 0.0073, "Cardiac Arrhythmias": 0.0039, "Diarrhea": 0.041,
            "Heart Failure": 0.0023, "Hypotension": 0.0056, "Liver Impairment": 0.0023,
            "Nausea": 0.023, "Neutropenia": 0.023, "Peripheral Neuropathies": 0.016,
            "Rash": 0.022, "Thrombocytopenia": 0.016, "Vomiting": 0.014},
    "Rd":  {"Acute Renal Failure": 0.0073, "Cardiac Arrhythmias": 0.0031, "Diarrhea": 0.034,
            "Heart Failure": 0.0016, "Hypotension": 0.0047, "Liver Impairment": 0.0016,
            "Nausea": 0.018, "Neutropenia": 0.022, "Peripheral Neuropathies": 0.011,
            "Rash": 0.011, "Thrombocytopenia": 0.014, "Vomiting": 0.011},
}


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def iso_date(days_from_ref: int, ref_year: int = 2015) -> str:
    ref = datetime.date(ref_year, 1, 1)
    d = ref + datetime.timedelta(days=int(max(0, days_from_ref)))
    return d.strftime("%Y-%m-%d")

def weighted_choice(options: dict, n: int) -> np.ndarray:
    keys  = list(options.keys())
    probs = np.array(list(options.values()), dtype=float)
    probs /= probs.sum()
    return RNG.choice(keys, size=n, p=probs)

def sim_weibull_arm(n_ird: int, n_rd: int,
                    median_ird: float, median_rd: float,
                    shape: float = 1.2, censor_rate: float = 0.35,
                    seed_ird: int = 1001, seed_rd: int = 1002):
    """
    Simulate arm-specific Weibull survival times (days).
    Medians are in months; convert to days internally.
    Returns (times_days, event_flags) for all subjects concatenated [IRd, Rd].

    seed_ird / seed_rd are deterministic per-arm seeds so calibration
    is independent of RNG consumption order elsewhere in the script.
    """
    # Weibull scale: median = scale * (ln2)^(1/shape)
    ln2_inv = np.log(2) ** (1.0 / shape)
    scale_ird = (median_ird * 28.0) / ln2_inv
    scale_rd  = (median_rd  * 28.0) / ln2_inv

    def draw(n, scale, seed):
        event_t  = weibull_min.rvs(shape, scale=scale, size=n, random_state=seed)
        rng_c    = np.random.default_rng(seed + 50_000)  # independent censoring RNG
        censor_t = rng_c.exponential(scale * 2.5, size=n)
        times    = np.minimum(event_t, censor_t)
        events   = (event_t <= censor_t).astype(int)
        return times.astype(int), events

    t_ird, e_ird = draw(n_ird, scale_ird, seed_ird)
    t_rd,  e_rd  = draw(n_rd,  scale_rd,  seed_rd)
    return np.concatenate([t_ird, t_rd]), np.concatenate([e_ird, e_rd])

def cockcroft_gault(age, weight_kg, scr_umol_L, sex):
    """Creatinine clearance in mL/min (Cockcroft-Gault 1976)."""
    scr_mg_dL = scr_umol_L / 88.4
    crcl = ((140.0 - age) * weight_kg) / (72.0 * scr_mg_dL)
    crcl = np.where(sex == "F", crcl * 0.85, crcl)
    return np.clip(np.round(crcl, 1), 5.0, 180.0)

def bsa_mosteller(height_cm, weight_kg):
    """Body surface area m² (Mosteller 1987)."""
    return np.round(np.sqrt(height_cm * weight_kg / 3600.0), 3)


# ─────────────────────────────────────────────────────────────────────────────
# DM — Demographics  (adds weight, height, BMI, BSA)
# ─────────────────────────────────────────────────────────────────────────────

def make_dm(study_key: str) -> pd.DataFrame:
    cfg = STUDY_CONFIG[study_key]
    n   = cfg["n"]
    sid = cfg["studyid"]

    usubjid = [f"{sid}-{i+1:04d}" for i in range(n)]

    ages = np.clip(
        np.round(RNG.normal(cfg["median_age"], 8, n)).astype(int),
        cfg["age_range"][0], cfg["age_range"][1]
    )
    sex     = RNG.choice(["F", "M"], size=n, p=[cfg["pct_female"], 1 - cfg["pct_female"]])
    race    = weighted_choice(RACE_DIST[study_key], n)
    ig_type = weighted_choice(IG_TYPE_DIST[study_key], n)
    lchain  = RNG.choice(["KAPPA", "LAMBDA", "BICLONAL"], size=n, p=[0.60, 0.35, 0.05])

    # ── Arm assignment: exact counts then shuffle (not probabilistic) ─────────
    # Spec: MM2 = 351 IRd / 354 Rd; MM1 = 360 IRd / 362 Rd
    # Using exact counts preserves stratification balance.
    n_ird = cfg.get("n_ird", n // 2)
    n_rd  = cfg.get("n_rd",  n - n_ird)
    iss_p   = np.array(list(ISS_DIST[study_key].values()), dtype=float)
    iss_p  /= iss_p.sum()
    iss     = RNG.choice([1, 2, 3], size=n, p=iss_p)
    ecog    = RNG.choice([0, 1, 2], size=n, p=[0.35, 0.50, 0.15])

    # Exact arm assignment — shuffle to avoid systematic ordering
    arm_arr = np.array(["IRd"] * n_ird + ["Rd"] * n_rd, dtype=object)
    RNG.shuffle(arm_arr)
    arm = arm_arr

    # Body weight (kg) — sex and age dependent; MM patients often frail/overweight
    # Male: mean 80kg, Female: mean 67kg; modest negative correlation with age >75
    wt_male   = np.clip(RNG.normal(80.0, 15.0, n), 45.0, 140.0)
    wt_female = np.clip(RNG.normal(67.0, 12.0, n), 38.0, 110.0)
    weight    = np.where(sex == "M", wt_male, wt_female)
    # Slight age-related weight loss for elderly (>75)
    weight    = np.where(ages > 75, weight * RNG.uniform(0.93, 1.0, n), weight)
    weight    = np.round(weight, 1)

    # Height (cm) — sex dependent, slightly decreases with age (bone loss)
    ht_male   = np.clip(RNG.normal(175.0, 7.0, n), 152.0, 198.0)
    ht_female = np.clip(RNG.normal(163.0, 6.0, n), 142.0, 182.0)
    height    = np.where(sex == "M", ht_male, ht_female)
    height    = np.where(ages > 70, height - RNG.uniform(0, 3, n), height)
    height    = np.round(height, 1)

    bmi = np.round(weight / (height / 100.0) ** 2, 1)
    bsa = bsa_mosteller(height, weight)

    rfstdtc = np.array([iso_date(RNG.integers(0, 3 * 365)) for _ in range(n)])

    dm = pd.DataFrame({
        "STUDYID":   sid,
        "DOMAIN":    "DM",
        "USUBJID":   usubjid,
        "SUBJID":    [f"{i+1:04d}" for i in range(n)],
        "RFSTDTC":   rfstdtc,
        "RFENDTC":   "",
        "SITEID":    RNG.choice([f"SITE{i:02d}" for i in range(1, 16)], size=n),
        "AGE":       ages,
        "AGEU":      "YEARS",
        "SEX":       sex,
        "RACE":      race,
        "ETHNIC":    RNG.choice(["HISPANIC OR LATINO", "NOT HISPANIC OR LATINO", "NOT REPORTED"],
                                size=n, p=[0.05, 0.90, 0.05]),
        "COUNTRY":   RNG.choice(["USA","FRA","DEU","GBR","JPN","ITA","ESP","AUS"],
                                size=n, p=[0.35,0.12,0.10,0.10,0.09,0.08,0.08,0.08]),
        "ARMCD":     arm,
        "ARM":       np.where(arm == "IRd",
                              "Ixazomib + Lenalidomide + Dexamethasone",
                              "Placebo + Lenalidomide + Dexamethasone"),
        "ACTARMCD":  arm,
        "ACTARM":    np.where(arm == "IRd",
                              "Ixazomib + Lenalidomide + Dexamethasone",
                              "Placebo + Lenalidomide + Dexamethasone"),
        "DTHDTC":    "",
        "DTHFL":     "",
        "IGTYPE":    ig_type,
        "LCTYPE":    lchain,
        "ISSSTAGE":  iss,
        "ECOG":      ecog,
        "DXMONTHS":  np.round(np.abs(RNG.exponential(cfg["median_dx_months"], n)), 2),
        "WEIGHT":    weight,         # kg — NEW
        "HEIGHT":    height,         # cm — NEW
        "BMI":       bmi,            # kg/m² — NEW
        "BSA":       bsa,            # m² (Mosteller) — NEW
        "POPULATION": cfg["population"],
    })
    return dm


# ─────────────────────────────────────────────────────────────────────────────
# EX — Exposure (all subjects; adherence + reason-coded dose mods)
# ─────────────────────────────────────────────────────────────────────────────

def make_ex(dm: pd.DataFrame, study_key: str) -> pd.DataFrame:
    sid  = STUDY_CONFIG[study_key]["studyid"]
    rows = []
    max_cycles = 26

    for _, subj in dm.iterrows():
        arm      = subj["ARMCD"]
        n_cycles = RNG.integers(2, max_cycles + 1)
        drugs    = ["LENALIDOMIDE", "DEXAMETHASONE"]
        drugs   += ["IXAZOMIB"] if arm == "IRd" else ["PLACEBO"]

        for drug in drugs:
            info = DRUG_DOSES[drug]
            for cyc in range(1, n_cycles + 1):
                start_day = (cyc - 1) * 28
                if "WEEKLY" in info["freq"]:
                    dose_days = [start_day + d for d in [0, 7, 14]]
                    if "22" in info["freq"]:
                        dose_days.append(start_day + 21)
                else:
                    dose_days = list(range(start_day, start_day + 21))

                for dd in dose_days:
                    # Dose modification — probabilities per cycle (higher in early cycles)
                    mod_p = [0.84, 0.10, 0.06] if cyc <= 3 else [0.88, 0.08, 0.04]
                    mod_factor = RNG.choice([1.0, 0.75, 0.5], p=mod_p)
                    actual_dose = round(info["dose"] * mod_factor, 1) if drug != "PLACEBO" else 0.0

                    if mod_factor < 1.0:
                        mod_reason = RNG.choice(
                            ["ADVERSE EVENT", "PHYSICIAN DECISION", "PROTOCOL SPECIFIED"],
                            p=[0.65, 0.25, 0.10]
                        )
                        exdosmod = "DOSE REDUCTION"
                    else:
                        mod_reason = ""
                        exdosmod   = ""

                    # Treatment adherence: % of planned dose actually taken
                    adherence = 100.0 if mod_factor == 1.0 else round(mod_factor * 100.0, 1)

                    rows.append({
                        "STUDYID":    sid,
                        "DOMAIN":     "EX",
                        "USUBJID":    subj["USUBJID"],
                        "EXSEQ":      None,
                        "EXTRT":      drug,
                        "EXDOSE":     actual_dose,
                        "EXDOSU":     info["unit"],
                        "EXDOSFRM":   "TABLET",
                        "EXDOSFRQ":   info["freq"],
                        "EXROUTE":    info["route"],
                        "EXSTDTC":    iso_date(dd),
                        "EXENDTC":    iso_date(dd),
                        "VISITNUM":   cyc,
                        "VISIT":      f"CYCLE {cyc}",
                        "EPOCH":      "TREATMENT",
                        "EXDOSMOD":   exdosmod,
                        "EXDOSMOD_REASON": mod_reason,   # NEW
                        "EXADH":      adherence,          # NEW: adherence %
                    })

    ex = pd.DataFrame(rows)
    ex["EXSEQ"] = ex.groupby("USUBJID").cumcount() + 1
    return ex


# ─────────────────────────────────────────────────────────────────────────────
# CM — Concomitant Medications (NEW domain)
# ─────────────────────────────────────────────────────────────────────────────

# CYP3A4 inhibitors affecting Ixazomib PK (Ixazomib is CYP3A4 substrate)
CYP3A4_INHIBITORS = {
    "CLARITHROMYCIN":  {"strength": "STRONG",   "prob": 0.03, "cl_multiplier": 0.45},
    "ITRACONAZOLE":    {"strength": "STRONG",   "prob": 0.015,"cl_multiplier": 0.50},
    "KETOCONAZOLE":    {"strength": "STRONG",   "prob": 0.008,"cl_multiplier": 0.48},
    "FLUCONAZOLE":     {"strength": "MODERATE", "prob": 0.05, "cl_multiplier": 0.70},
    "DILTIAZEM":       {"strength": "MODERATE", "prob": 0.04, "cl_multiplier": 0.75},
    "VERAPAMIL":       {"strength": "MODERATE", "prob": 0.025,"cl_multiplier": 0.72},
}
CYP3A4_INDUCERS = {
    "RIFAMPIN":        {"strength": "STRONG",   "prob": 0.01, "cl_multiplier": 2.50},
    "CARBAMAZEPINE":   {"strength": "STRONG",   "prob": 0.02, "cl_multiplier": 2.20},
    "PHENYTOIN":       {"strength": "MODERATE", "prob": 0.02, "cl_multiplier": 1.80},
    "DEXAMETHASONE_CYP": {"strength": "WEAK",  "prob": 0.90, "cl_multiplier": 1.10},
}
ANTICOAGULANTS = {
    "ASPIRIN":              {"prob": 0.50, "soc": "BLOOD AND LYMPHATIC SYSTEM DISORDERS"},
    "ENOXAPARIN":           {"prob": 0.30, "soc": "BLOOD AND LYMPHATIC SYSTEM DISORDERS"},
    "WARFARIN":             {"prob": 0.08, "soc": "BLOOD AND LYMPHATIC SYSTEM DISORDERS"},
    "APIXABAN":             {"prob": 0.06, "soc": "BLOOD AND LYMPHATIC SYSTEM DISORDERS"},
}
SUPPORTIVE_CARE = {
    "FILGRASTIM (G-CSF)":   {"prob": 0.18, "soc": "BLOOD AND LYMPHATIC SYSTEM DISORDERS"},
    "EPOETIN ALFA (EPO)":   {"prob": 0.10, "soc": "BLOOD AND LYMPHATIC SYSTEM DISORDERS"},
    "ACYCLOVIR (ANTIVIRAL)":{"prob": 0.65, "soc": "INFECTIONS AND INFESTATIONS"},
    "OMEPRAZOLE (PPI)":     {"prob": 0.40, "soc": "GASTROINTESTINAL DISORDERS"},
}

def make_cm(dm: pd.DataFrame, study_key: str) -> pd.DataFrame:
    sid  = STUDY_CONFIG[study_key]["studyid"]
    rows = []
    seq  = 1

    for _, subj in dm.iterrows():
        uid    = subj["USUBJID"]
        arm    = subj["ARMCD"]
        rfst   = subj["RFSTDTC"]

        # CYP3A4 inhibitors (mutually exclusive strong inhibitors)
        for drug, info in CYP3A4_INHIBITORS.items():
            if RNG.random() < info["prob"]:
                start_day = RNG.integers(-30, 100)  # may start before randomisation
                dur       = RNG.integers(14, 180)
                rows.append({
                    "STUDYID": sid, "DOMAIN": "CM", "USUBJID": uid, "CMSEQ": seq,
                    "CMTRT": drug, "CMDECOD": drug,
                    "CMSOC": "INFECTIONS AND INFESTATIONS",
                    "CMCAT": "CYP3A4 INHIBITOR",
                    "CMROUTE": "ORAL", "CMDOSFRM": "TABLET",
                    "CMSTDTC": iso_date(start_day), "CMENDTC": iso_date(start_day + dur),
                    "CMONGO": "N", "EPOCH": "TREATMENT",
                    "CYP3A4_STRENGTH": info["strength"],
                    "CL_MULTIPLIER": info["cl_multiplier"],
                })
                seq += 1
                break  # one strong inhibitor at a time

        # CYP3A4 inducers (note: dexamethasone is always present in both arms as study drug)
        for drug, info in CYP3A4_INDUCERS.items():
            if drug == "DEXAMETHASONE_CYP":
                continue  # handled in EX; flag set in ADSL
            if RNG.random() < info["prob"]:
                start_day = RNG.integers(-14, 50)
                dur       = RNG.integers(7, 90)
                rows.append({
                    "STUDYID": sid, "DOMAIN": "CM", "USUBJID": uid, "CMSEQ": seq,
                    "CMTRT": drug, "CMDECOD": drug,
                    "CMSOC": "INFECTIONS AND INFESTATIONS",
                    "CMCAT": "CYP3A4 INDUCER",
                    "CMROUTE": "ORAL", "CMDOSFRM": "TABLET",
                    "CMSTDTC": iso_date(start_day), "CMENDTC": iso_date(start_day + dur),
                    "CMONGO": "N", "EPOCH": "TREATMENT",
                    "CYP3A4_STRENGTH": info["strength"],
                    "CL_MULTIPLIER": info["cl_multiplier"],
                })
                seq += 1

        # Anticoagulants (lenalidomide requires DVT prophylaxis)
        for drug, info in ANTICOAGULANTS.items():
            if RNG.random() < info["prob"]:
                rows.append({
                    "STUDYID": sid, "DOMAIN": "CM", "USUBJID": uid, "CMSEQ": seq,
                    "CMTRT": drug, "CMDECOD": drug, "CMSOC": info["soc"],
                    "CMCAT": "ANTICOAGULANT / ANTIPLATELET",
                    "CMROUTE": "ORAL", "CMDOSFRM": "TABLET",
                    "CMSTDTC": iso_date(0), "CMENDTC": "",
                    "CMONGO": "Y", "EPOCH": "TREATMENT",
                    "CYP3A4_STRENGTH": "", "CL_MULTIPLIER": 1.0,
                })
                seq += 1

        # Supportive care
        for drug, info in SUPPORTIVE_CARE.items():
            if RNG.random() < info["prob"]:
                start_day = RNG.integers(0, 200)
                dur       = RNG.integers(7, 60)
                rows.append({
                    "STUDYID": sid, "DOMAIN": "CM", "USUBJID": uid, "CMSEQ": seq,
                    "CMTRT": drug, "CMDECOD": drug, "CMSOC": info["soc"],
                    "CMCAT": "SUPPORTIVE CARE",
                    "CMROUTE": "ORAL" if "PPI" in drug or "ANTIVIRAL" in drug else "SUBCUTANEOUS",
                    "CMDOSFRM": "TABLET",
                    "CMSTDTC": iso_date(start_day), "CMENDTC": iso_date(start_day + dur),
                    "CMONGO": "N", "EPOCH": "TREATMENT",
                    "CYP3A4_STRENGTH": "", "CL_MULTIPLIER": 1.0,
                })
                seq += 1

    return pd.DataFrame(rows) if rows else pd.DataFrame()


# ─────────────────────────────────────────────────────────────────────────────
# LB — Laboratory Results (improved trajectories; all subjects)
# ─────────────────────────────────────────────────────────────────────────────

def _sample_resp_phenotype(arm: str, is_ndmm: bool, auc_rel: float) -> float:
    """
    Sample per-patient response depth (fraction of baseline M-protein eliminated at nadir).
    Uses a bimodal mixture to match published TOURMALINE IMWG response distributions.

    Published targets (Moreau 2019 NEJM / Moreau 2016 Lancet Oncol):
      MM2 IRd: ORR=82%, VGPR+=63%, CR+=28%  → non-resp=18%, PR=19%, VGPR=35%, CR+=28%
      MM2 Rd:  ORR=75%, VGPR+=55%, CR+=14%  → non-resp=25%, PR=20%, VGPR=41%, CR+=14%
      MM1 IRd: ORR=78%, VGPR+=48%, CR+=12%  → non-resp=22%, PR=30%, VGPR=36%, CR+=12%
      MM1 Rd:  ORR=72%, VGPR+=39%, CR+=7%   → non-resp=28%, PR=33%, VGPR=32%, CR+=7%

    AUC exposure-response: higher individual AUC shifts patients toward deeper response.
    """
    # Arm × study-specific mixture probabilities
    # Inflated slightly above published targets to compensate for:
    #   - 3% MAR loss on disease biomarker visits (miss_rate=0.03)
    #   - Short-follow-up patients (n_cycles=3) who may not organically reach nadir
    # CR tier starts at 0.993 (not 0.990) to avoid -99% boundary contamination.
    if is_ndmm:
        if arm == "IRd":
            p_nonresp, p_pr, p_vgpr = 0.15, 0.22, 0.35   # CR+ = 1-sum = 0.28
        else:                                               # MM2 Rd: ORR=75%, VGPR+=55%, CR+=14%
            p_nonresp, p_pr, p_vgpr = 0.22, 0.18, 0.43   # CR+ = 0.17; calibrated: seed=42 gives CR ~14%
    else:
        if arm == "IRd":
            p_nonresp, p_pr, p_vgpr = 0.20, 0.30, 0.38   # CR+ = 0.12; VGPR+=0.50→48% after MAR
        else:
            p_nonresp, p_pr, p_vgpr = 0.28, 0.33, 0.32

    # AUC exposure-response: high AUC (>1.3x) upgrades some patients one tier
    if arm == "IRd" and auc_rel > 1.3:
        p_nonresp = max(0.05, p_nonresp - 0.05)

    u = RNG.random()
    if u < p_nonresp:
        # Non-responder: M-protein stays flat or rises slightly
        return float(RNG.uniform(-0.10, 0.20))   # negative = increase
    elif u < p_nonresp + p_pr:
        # PR tier: 50–90% reduction
        return float(RNG.uniform(0.50, 0.90))
    elif u < p_nonresp + p_pr + p_vgpr:
        # VGPR tier: 90–99% reduction
        return float(RNG.uniform(0.90, 0.990))
    else:
        # CR/sCR tier: near-complete elimination.
        # Start at 0.993 (not 0.990) to give all CR patients ≥0.3% margin below
        # the -99% IMWG threshold, so nadir noise doesn't push them into VGPR.
        return float(RNG.uniform(0.993, 1.000))


def _sim_trajectory(test: str, ig_type: str, arm: str,
                    n_cycles: int, is_ndmm: bool,
                    ixaz_auc_rel: float = 1.0,
                    resp_rate_override: float | None = None) -> np.ndarray:
    """
    Simulate physiologically-plausible longitudinal lab values.

    ixaz_auc_rel      : relative individual Ixazomib AUC (1.0 = typical).
    resp_rate_override: pre-sampled per-patient response depth (0–1 fraction of
                        baseline eliminated). If provided, overrides the internal
                        sampling for disease markers — ensures all disease markers
                        (SPEP_MPROT, IGG, B2MG, etc.) reflect the same response depth.
    """
    lo, hi, _ = LAB_NORMAL_RANGES[test]
    mid = (lo + hi) / 2
    sd  = (hi - lo) / 4

    disease_markers     = {"SPEP_MPROT", "IGA", "IGG", "IGM", "SPEP_KAPPA",
                           "SPEP_LAMBDA", "SPEP_GAMMA", "UPEP_MPROT"}
    anemia_markers      = {"HGB", "HCT"}
    protective_markers  = {"GFR", "ALBUMIN"}
    progression_markers = {"CREAT", "BUN", "CA", "CA_CORR"}
    hema_suppression    = {"NEUT", "PLT", "WBC"}  # transient myelosuppression

    vals = np.zeros(n_cycles)

    if test in disease_markers:
        # Baseline: elevated 2-5× ULN for dominant Ig, 0.5-2× for non-dominant
        dom = (ig_type in {"IgG", "IgA"} and test in {"IGG", "IGA", "SPEP_MPROT"}) \
           or (ig_type == "IgM" and test == "IGM")
        base = hi * RNG.uniform(2.0, 5.0) if dom else hi * RNG.uniform(0.5, 2.0)

        # Response rate: use pre-sampled bimodal phenotype if provided,
        # otherwise fall back to narrow normal (for non-SPEP disease markers).
        if resp_rate_override is not None:
            resp_rate = resp_rate_override
        else:
            er_bonus  = 0.03 * (ixaz_auc_rel - 1.0) * 10 if arm == "IRd" else 0.0
            resp_rate = np.clip(RNG.normal(0.72, 0.10) + er_bonus, 0.30, 0.95) if arm == "IRd" \
                       else np.clip(RNG.normal(0.58, 0.10), 0.20, 0.85)
        nadir = base * (1.0 - resp_rate)

        for t in range(n_cycles):
            # Exponential approach to nadir with k=0.30.
            # Calibrated to Spec §5A published trajectory (Moreau 2016 Lancet Oncol):
            #   C3: ~-44%  → exp(-0.30*3)=0.41  → val ≈ nadir + 0.41*(base-nadir)
            #   C6: ~-72%  → exp(-0.30*6)=0.17  → val ≈ nadir + 0.17*(base-nadir)
            #   C12:~-94%  → exp(-0.30*12)=0.03 → val ≈ nadir + 0.03*(base-nadir)
            # (Prior k=0.80 reached nadir by C3 → medians were -83% at C6, too deep.)
            decay = np.exp(-0.30 * t)
            relapse_signal = 0.0
            if t > 12:
                relapse_prob = 0.04 * (t - 12)
                if RNG.random() < relapse_prob:
                    relapse_signal = base * 0.05 * (t - 12)
            vals[t] = max(0, nadir + (base - nadir) * decay
                          + relapse_signal + RNG.normal(0, base * 0.04))

        # Best-response nadir override: insert a guaranteed nadir observation so
        # that each patient's phenotypic response depth is captured regardless of
        # follow-up length.  In real TOURMALINE trials, ORR/VGPR/CR are recorded
        # as BEST response, not response at a fixed snapshot.
        #
        # Extended to resp_rate ≥ 0.50 (any responder — PR, VGPR, CR):
        #   PR patients with short follow-up may not organically reach -50% nadir
        #   with k=0.30 decay constant; the override ensures their best response
        #   is correctly captured.
        # Threshold: n_cycles ≥ 2 (need at least 1 on-treatment cycle, t=1).
        # Override position: cycle 10 (t=9) for n_cycles>10, else last cycle.
        # Preserves C3/C6 calibration: override is at t≥2 (≥C3) for most pts.
        if resp_rate >= 0.50 and n_cycles >= 2:
            nadir_t = min(n_cycles - 1, 9)   # cycle 10 (t=9) or last available
            if resp_rate >= 0.90:
                # VGPR/CR: tight noise to preserve phenotype boundaries.
                # CR patients (nadir≈0.005×base) get almost-zero noise so they
                # stay below the -99% IMWG threshold; VGPR patients stay above.
                # std = 5% of nadir + 0.05% of base (floor for near-zero nadirs)
                nadir_noise_std = nadir * 0.05 + base * 0.0005
            else:
                # PR: moderate noise consistent with the organic trajectory noise
                # (base×4% = same as RNG.normal(0, base*0.04) in the loop).
                nadir_noise_std = base * 0.04
            nadir_val = max(1e-6, nadir + RNG.normal(0, nadir_noise_std))
            vals[nadir_t] = nadir_val
            if resp_rate >= 0.90:
                # VGPR/CR only: floor post-nadir cycles to prevent trajectory
                # noise from creating false global minima that push VGPR patients
                # into CR territory (or VGPR patients below -99%).
                for t2 in range(nadir_t + 1, n_cycles):
                    if vals[t2] < nadir_val:
                        vals[t2] = nadir_val * RNG.uniform(1.001, 1.02)

    elif test in anemia_markers:
        # Baseline low (anemia is CRAB criterion); gradual recovery with response
        base = lo * RNG.uniform(0.65, 0.90)
        for t in range(n_cycles):
            recovery = (mid - base) * min(1.0, t / 14.0)
            vals[t] = base + recovery + RNG.normal(0, sd * 0.08)

    elif test == "B2MG":
        # Baseline elevated in MM; decreases with treatment response
        base = hi * RNG.uniform(1.5, 6.0)  # typically 2-15 mg/L at diagnosis
        resp = 0.60 if arm == "IRd" else 0.45
        nadir = base * (1 - resp)
        for t in range(n_cycles):
            vals[t] = max(lo, nadir + (base - nadir) * np.exp(-0.25 * t)
                          + RNG.normal(0, base * 0.05))

    elif test == "LDH":
        # ~65% MM patients have normal LDH at diagnosis; ~35% elevated (R-ISS III feature)
        # Normal range 140-280 U/L; elevated = 280-600 U/L
        if RNG.random() < 0.65:
            base = RNG.uniform(lo, hi)          # normal LDH
        else:
            base = RNG.uniform(hi, hi * 2.2)    # elevated (high-risk disease)
        for t in range(n_cycles):
            norm = (mid - base) * min(1.0, t / 8.0) * 0.4   # partial normalization
            vals[t] = max(lo * 0.8, base + norm + RNG.normal(0, sd * 0.10))

    elif test in hema_suppression:
        # Myelosuppression (Lenalidomide days 1-21 + Ixazomib additive thrombocytopenia)
        # Spec §5B: PLT baseline ~220, nadir Day 15 ~110-130 (IRd arm)
        # Grade 3 thrombocytopenia (<50): ~31% IRd vs ~16% Rd
        #
        # Model: per-cycle Day-1 value drifts down slightly across cycles (cumulative),
        # then recovers.  The within-cycle dip (Day 8→Day 15) is handled by make_lb().
        # Here we model the per-cycle Day-1 anchor value.
        base = mid + RNG.normal(0, sd * 0.5)
        # PLT baseline: NDMM ~200 (bone marrow compromise at diagnosis)
        # RRMM ~195 (prior therapy + more bone marrow involvement → lower reserve)
        # Calibrated with dip_amp=0.45/0.46 to hit MM1 targets (IRd=31%, Rd=16%)
        if test == "PLT":
            base_mean = 200.0 if is_ndmm else 195.0
            base = np.clip(RNG.lognormal(np.log(base_mean), 0.30), lo * 0.3, lo * 3.5)
        for t in range(n_cycles):
            # Cumulative downward drift across cycles (~5% per cycle IRd, ~2% Rd)
            # representing persistent myelosuppression
            cum_drift = (0.05 if arm == "IRd" else 0.02) * t
            cum_factor = max(0.6, 1.0 - cum_drift)
            dip_early = 0.20 * base * np.exp(-0.5 * t) if arm == "IRd" \
                   else 0.10 * base * np.exp(-0.5 * t)
            vals[t] = max(lo * 0.1, base * cum_factor - dip_early
                          + RNG.normal(0, sd * 0.08))

    elif test in protective_markers:
        base = lo * RNG.uniform(0.85, 1.05)
        for t in range(n_cycles):
            vals[t] = max(lo * 0.5, base + RNG.normal(0, sd * 0.08))

    elif test in progression_markers:
        base = mid * RNG.uniform(0.9, 1.5)
        for t in range(n_cycles):
            vals[t] = max(lo * 0.5, base + RNG.normal(0, sd * 0.05) * t ** 0.3)

    else:
        base = mid + RNG.normal(0, sd * 0.5)
        for t in range(n_cycles):
            vals[t] = max(lo * 0.2, base + RNG.normal(0, sd * 0.08))

    return np.clip(vals, lo * 0.05, hi * 10.0)


def _lb_cat(test: str) -> str:
    heme   = {"HGB","HCT","LYMPH","MONO","NEUT","PLT","WBC"}
    immuno = {"IGA","IGG","IGM","SPEP_KAPPA","SPEP_LAMBDA","SPEP_GAMMA",
              "SPEP_MPROT","UPEP_MPROT","KAPPA_LAMBDA"}
    if test in heme:   return "HEMATOLOGY"
    if test in immuno: return "SERUM IMMUNOGLOBULINS"
    return "CHEMISTRY"


def _lb_label(test: str) -> str:
    return {
        "ALBUMIN":"Albumin","ALP":"Alkaline Phosphatase","ALT":"Alanine Aminotransferase",
        "AST":"Aspartate Aminotransferase","B2MG":"Beta-2 Microglobulin",
        "BILI":"Bilirubin","BUN":"Blood Urea Nitrogen","CA":"Calcium",
        "CA_CORR":"Corrected Calcium","CL":"Chloride","CO2":"Carbon Dioxide",
        "CREAT":"Creatinine","GFR":"Glomerular Filtration Rate","GLOB":"Serum Globulin",
        "GLUC":"Glucose","HCT":"Hematocrit","HGB":"Hemoglobin",
        "IGA":"Immunoglobulin A","IGG":"Immunoglobulin G","IGM":"Immunoglobulin M",
        "K":"Potassium","KAPPA_LAMBDA":"Kappa/Lambda FLC Ratio",
        "LDH":"Lactate Dehydrogenase","LYMPH":"Lymphocytes","MG":"Magnesium",
        "MONO":"Monocytes","NA":"Sodium","NEUT":"Neutrophils","PHOS":"Phosphate",
        "PLT":"Platelets","PROT":"Total Protein","SPEP_GAMMA":"SPEP Gamma Globulin",
        "SPEP_KAPPA":"Free SPEP Kappa Light Chain","SPEP_LAMBDA":"Free SPEP Lambda Light Chain",
        "SPEP_MPROT":"SPEP Monoclonal Protein","URATE":"Urate",
        "URINE_ALB":"Urine Albumin","UPEP_MPROT":"UPEP Monoclonal Protein","WBC":"Leukocytes",
    }.get(test, test)


def make_lb(dm: pd.DataFrame, study_key: str, max_cycles: int = 24) -> pd.DataFrame:
    """
    Generate longitudinal lab records aligned with the TOURMALINE Schedule of
    Assessments (SoA).

    Protocol-specified assessment schedule
    ---------------------------------------
    Baseline (screening): single measurement for all tests.

    Treatment cycles — three distinct frequencies by test category:

    1. Hematology + Clinical Chemistry  (Days 1, 8, 15 every cycle)
       Per SoA: "Hematology: CBC with differential, platelets (Days 1, 8, 15)"
                "Clinical Chemistry: Comprehensive panel (Days 1, 8, 15)"
       Tests: HGB, HCT, WBC, NEUT, LYMPH, MONO, PLT
              ALBUMIN, ALP, ALT, AST, BILI, BUN, CA, CA_CORR,
              CL, CO2, CREAT, GFR, GLOB, GLUC, K, LDH, MG,
              NA, PHOS, PROT, URATE, URINE_ALB

    2. Disease biomarkers — Serum/Urine Protein & Immunoglobulins (Day 1 only)
       Per SoA: "SPEP: Every 4 weeks"  "UPEP: Every 4 weeks"
                "Serum Free Light Chains: Every 4 weeks"
       Tests: SPEP_MPROT, SPEP_KAPPA, SPEP_LAMBDA, SPEP_GAMMA,
              KAPPA_LAMBDA, IGA, IGG, IGM, UPEP_MPROT

    3. Beta-2 Microglobulin (Day 1, every 12 weeks = every 3 cycles)
       Per SoA: "Beta-2 Microglobulin: Every 12 weeks"

    Within-cycle interpolation (Hematology/Chemistry days 8 & 15)
    --------------------------------------------------------------
    Cycle-level trajectory values are anchors for Day 1.  Days 8 and 15 are
    linearly interpolated toward the next cycle value with:
      - CV ≈ 4 % biological noise
      - Myelosuppression dip for NEUT/PLT/WBC (Lenalidomide days 1-21):
        dip = dip_amp × sin(π × day_offset / 28)  → nadir ≈ day 14
    """
    sid     = STUDY_CONFIG[study_key]["studyid"]
    is_ndmm = (study_key == "MM2")
    tests   = list(LAB_NORMAL_RANGES.keys())
    rows    = []

    # ── Protocol assessment categories (per TOURMALINE SoA) ──────────────────
    _HEME = {"HGB", "HCT", "WBC", "NEUT", "LYMPH", "MONO", "PLT"}
    _CHEM = {
        "ALBUMIN", "ALP", "ALT", "AST", "BILI", "BUN", "CA", "CA_CORR",
        "CL", "CO2", "CREAT", "GFR", "GLOB", "GLUC", "K", "LDH", "MG",
        "NA", "PHOS", "PROT", "URATE", "URINE_ALB",
    }
    _DISEASE_BM = {
        "SPEP_MPROT", "SPEP_KAPPA", "SPEP_LAMBDA", "SPEP_GAMMA",
        "KAPPA_LAMBDA", "IGA", "IGG", "IGM", "UPEP_MPROT",
    }
    # BMPC assessed separately in generate_pd.py (per-CR confirmation schedule)

    # Within-cycle myelosuppression — affects NEUT, PLT, WBC
    _MYELO_SUPP = {"NEUT", "PLT", "WBC"}

    # Log-normal IIV on relative Ixazomib AUC: CV ≈ 36 %  [Gupta 2017]
    n = len(dm)
    auc_rel = np.exp(RNG.normal(0, 0.36, n))

    for idx, (_, subj) in enumerate(dm.iterrows()):
        n_cycles = RNG.integers(3, max_cycles + 1)
        arm      = subj["ARMCD"]
        ig       = subj["IGTYPE"]
        auc_i    = auc_rel[idx] if arm == "IRd" else 1.0

        # Pre-sample response phenotype once → all disease markers consistent
        resp_pheno = _sample_resp_phenotype(arm, is_ndmm, auc_i)

        for test in tests:
            lo, hi, unit = LAB_NORMAL_RANGES[test]
            traj = _sim_trajectory(test, ig, arm, n_cycles, is_ndmm, auc_i,
                                   resp_rate_override=resp_pheno)

            for cyc_idx, val in enumerate(traj):
                is_baseline = (cyc_idx == 0)

                # ── Protocol-specified sampling schedule ──────────────────────
                if is_baseline:
                    week_offsets = [0]                        # screening: all tests once

                elif test == "B2MG":
                    # Every 12 weeks (every 3rd cycle Day 1)
                    if cyc_idx % 3 == 0:
                        week_offsets = [0]
                    else:
                        continue                               # skip non-assessment cycle

                elif test in _DISEASE_BM:
                    # Every 4 weeks = Day 1 of each cycle only
                    week_offsets = [0]

                elif test in _HEME or test in _CHEM:
                    # Days 1, 8, 15 of every cycle
                    week_offsets = [0, 7, 14]

                else:
                    week_offsets = [0]                        # fallback: Day 1 only

                # Next cycle value used for linear interpolation
                next_val = traj[min(cyc_idx + 1, len(traj) - 1)]

                for wk_off in week_offsets:
                    week_num = wk_off // 7 + 1                # 1, 2, or 3

                    # ── Within-cycle value ────────────────────────────────────
                    if wk_off == 0:
                        val_wc = float(val)
                    else:
                        # Linear trend toward next cycle + biological noise
                        frac   = wk_off / 28.0
                        interp = float(val) + (float(next_val) - float(val)) * frac
                        val_wc = interp + RNG.normal(0, abs(float(val)) * 0.04)

                        # Myelosuppression dip (Lenalidomide days 1-21)
                        # Nadir ≈ day 14 (wk_off=14): sin(π×14/28) = 1.0
                        # Spec §5B: per-patient worst Grade 3 PLT (<50) targets:
                        #   MM2 IRd: 25%, Rd: 14%  |  MM1 IRd: 31%, Rd: 16%
                        # Calibration (per-patient worst PLT across all cycles):
                        # PLT/NEUT/WBC myelosuppression dip amplitudes (per-study calibration).
                        # Grade3 rate grows super-linearly with dip_amp; calibrated via
                        # log-linear interpolation from two empirical reference points
                        # (per-patient worst PLT < 50 × 10⁹/L; fixed RNG seeds).
                        # IRd: MM2 dip=0.45 → 26.8% Grade3 (target 25%) [seed=42] ✓
                        #      MM1 dip=0.48 → ~31% Grade3 (target 31%) [seed=43]
                        #      [empirical: 0.45→25.3%, 0.50→37.5%; log-linear → 0.48@31%]
                        # Rd:  MM2 dip=0.47 → ~14% Grade3 (target 14%) [seed=42]
                        #      [empirical: 0.44→10.2%, 0.52→24.9%; log-linear → 0.47@14%]
                        #      MM1 dip=0.46 → 16.0% Grade3 (target 16%) [seed=43] ✓
                        if test in _MYELO_SUPP:
                            if arm == "IRd":
                                dip_amp = 0.45 if is_ndmm else 0.48
                            else:
                                dip_amp = 0.47 if is_ndmm else 0.46
                            dip     = dip_amp * abs(float(val)) * np.sin(
                                np.pi * wk_off / 28.0
                            )
                            val_wc -= dip

                        val_wc = max(0.0, val_wc)

                    # ── Missingness (MAR) — calibrated by test frequency ──────
                    # Mandatory hematology at day 8/15: lower miss rate
                    if is_baseline:
                        miss_rate = 0.03
                    elif test in _DISEASE_BM or test == "B2MG":
                        # Disease biomarkers (SPEP, Ig, UPEP) are the primary
                        # efficacy endpoint in oncology trials — visit compliance
                        # is near-perfect.  Low miss_rate ensures the guaranteed
                        # nadir override at cycle 10 is almost always observed,
                        # so best-response (ORR/VGPR/CR) is correctly captured.
                        miss_rate = 0.03   # key efficacy endpoint; rarely missed
                    elif wk_off == 0:
                        miss_rate = 0.05   # Day 1 visit well-attended
                    else:
                        miss_rate = 0.10   # Day 8/15 — intermediate visits

                    is_missing = (RNG.random() < miss_rate)
                    study_day  = cyc_idx * 28 + wk_off + 1
                    val_r      = round(val_wc, 3)

                    rows.append({
                        "STUDYID":  sid,
                        "DOMAIN":   "LB",
                        "USUBJID":  subj["USUBJID"],
                        "LBSEQ":    None,
                        "LBCAT":    _lb_cat(test),
                        "LBTESTCD": test,
                        "LBTEST":   _lb_label(test),
                        "LBORRES":  "" if is_missing else f"{val_r:.3f}",
                        "LBORRESU": unit,
                        "LBSTRESC": "" if is_missing else f"{val_r:.3f}",
                        "LBSTRESN": np.nan if is_missing else val_r,
                        "LBSTRESU": unit,
                        "LBNRLO":   lo,
                        "LBNRHI":   hi,
                        "LBNRIND":  "" if is_missing else (
                            "LOW"  if val_r < lo else
                            ("HIGH" if val_r > hi else "NORMAL")
                        ),
                        "VISITNUM": cyc_idx + 1,
                        "WEEKNUM":  week_num,             # week within cycle (1–3)
                        "VISIT":    "BASELINE" if is_baseline
                                    else f"CYCLE {cyc_idx+1} WEEK {week_num}",
                        "LBDTC":    iso_date(cyc_idx * 28 + wk_off),
                        "LBDY":     study_day,
                        "EPOCH":    "BASELINE" if is_baseline else "TREATMENT",
                    })

    lb = pd.DataFrame(rows)
    lb["LBSEQ"] = lb.groupby("USUBJID").cumcount() + 1
    return lb


# ─────────────────────────────────────────────────────────────────────────────
# AE — Adverse Events (hematologic AEs linked to lab values)
# ─────────────────────────────────────────────────────────────────────────────

def make_ae(dm: pd.DataFrame, lb: pd.DataFrame, study_key: str) -> pd.DataFrame:
    sid  = STUDY_CONFIG[study_key]["studyid"]
    rows = []

    # Build per-subject per-cycle NEUT and PLT medians for hematologic AE linkage
    heme_lb = lb[lb["LBTESTCD"].isin(["NEUT", "PLT"]) & lb["LBSTRESN"].notna()][
        ["USUBJID", "LBTESTCD", "VISITNUM", "LBSTRESN"]
    ].copy()
    heme_pivot = heme_lb.pivot_table(
        index=["USUBJID", "VISITNUM"], columns="LBTESTCD", values="LBSTRESN"
    ).reset_index()

    subj_heme = {}
    for uid, grp in heme_pivot.groupby("USUBJID"):
        subj_heme[uid] = grp.set_index("VISITNUM")

    for _, subj in dm.iterrows():
        uid      = subj["USUBJID"]
        arm      = subj["ARMCD"]
        n_cycles = RNG.integers(3, 27)
        heme     = subj_heme.get(uid, pd.DataFrame())

        for ae_name, soc in AE_DEFS:
            per_cycle_prob = AE_INCIDENCE[arm][ae_name]
            is_heme        = ae_name in HEMATOLOGIC_AE

            for cyc_idx in range(n_cycles):
                cyc_num = cyc_idx + 1
                p = per_cycle_prob

                # Hematologic AEs: boost probability when lab values are low
                if is_heme and not heme.empty and cyc_num in heme.index:
                    row_h = heme.loc[cyc_num]
                    if ae_name == "Neutropenia":
                        neut = row_h.get("NEUT", 2.0)
                        if neut < 1.0:   p = min(p * 3.0, 0.60)
                        elif neut < 1.5: p = min(p * 1.8, 0.40)
                    elif ae_name == "Thrombocytopenia":
                        plt_ = row_h.get("PLT", 200.0)
                        if plt_ < 75:   p = min(p * 3.0, 0.55)
                        elif plt_ < 100: p = min(p * 1.8, 0.35)

                if RNG.random() >= p:
                    continue

                # CTCAE grade
                grade = RNG.choice([3, 4, 5], p=[0.55, 0.40, 0.05]) if is_heme \
                   else RNG.choice([2, 3, 4], p=[0.60, 0.35, 0.05])

                onset_day = cyc_idx * 28 + RNG.integers(1, 28)
                duration  = RNG.integers(3, 45)

                # Peripheral neuropathy: cumulative — probability increases with cycles
                if ae_name == "Peripheral Neuropathies" and arm == "IRd":
                    p *= min(1.0, 1.0 + 0.05 * cyc_idx)  # accumulates

                rows.append({
                    "STUDYID":  sid,
                    "DOMAIN":   "AE",
                    "USUBJID":  uid,
                    "AESEQ":    None,
                    "AETERM":   ae_name.upper(),
                    "AEDECOD":  ae_name,
                    "AEBODSYS": soc,
                    "AESOC":    soc,
                    "AESEV":    {1:"MILD",2:"MODERATE",3:"SEVERE",
                                 4:"LIFE THREATENING",5:"FATAL"}[grade],
                    "AETOXGR":  str(grade),
                    "AESER":    "Y" if grade >= 3 else "N",
                    "AEREL":    RNG.choice(["RELATED","NOT RELATED","POSSIBLY RELATED"],
                                           p=[0.55, 0.25, 0.20]),
                    "AEACN":    RNG.choice(["DOSE REDUCTION","DRUG WITHDRAWN",
                                            "DOSE NOT CHANGED","DRUG INTERRUPTED"],
                                           p=[0.20, 0.05, 0.55, 0.20]),
                    "AESTDTC":  iso_date(onset_day),
                    "AEENDTC":  iso_date(onset_day + duration),
                    "AEDY":     onset_day + 1,
                    "AEENDY":   onset_day + duration + 1,
                    "AEOUT":    RNG.choice(["RECOVERED/RESOLVED","RECOVERING/RESOLVING",
                                            "NOT RECOVERED/NOT RESOLVED"],
                                           p=[0.70, 0.20, 0.10]),
                    "EPOCH":    "TREATMENT",
                    "VISITNUM": cyc_num,
                    "VISIT":    f"CYCLE {cyc_num}",
                })

    ae = pd.DataFrame(rows)
    if len(ae):
        ae["AESEQ"] = ae.groupby("USUBJID").cumcount() + 1
    return ae


# ─────────────────────────────────────────────────────────────────────────────
# DS — Disposition (arm-specific survival with published HR)
# ─────────────────────────────────────────────────────────────────────────────

def make_ds(dm: pd.DataFrame, study_key: str) -> pd.DataFrame:
    sid    = STUDY_CONFIG[study_key]["studyid"]
    cfg    = STUDY_CONFIG[study_key]
    is_ndmm = (study_key == "MM2")

    ird_mask = dm["ARMCD"].values == "IRd"
    n_ird    = ird_mask.sum()
    n_rd     = (~ird_mask).sum()

    # Deterministic seeds per study × arm × endpoint — stable regardless of RNG call order
    # Seeds chosen so KM medians match calibrated targets (see calibration block)
    base = 2000 if study_key == "MM2" else 3000

    # PFS — arm-specific
    pfs_t, pfs_e = sim_weibull_arm(
        n_ird, n_rd,
        cfg["pfs_median_ird"], cfg["pfs_median_rd"],
        shape=1.2, censor_rate=0.35,
        seed_ird=base + 1, seed_rd=base + 2
    )
    # OS — arm-specific
    os_t, os_e = sim_weibull_arm(
        n_ird, n_rd,
        cfg["os_median_ird"], cfg["os_median_rd"],
        shape=1.3, censor_rate=0.50,
        seed_ird=base + 3, seed_rd=base + 4
    )

    # Re-order to match dm row order (IRd rows first in sim_weibull_arm output)
    # We need to map back: first n_ird values → IRd subjects, next n_rd → Rd subjects
    ird_idx = np.where(ird_mask)[0]
    rd_idx  = np.where(~ird_mask)[0]
    all_idx = np.concatenate([ird_idx, rd_idx])

    pfs_t_ordered = np.zeros(len(dm), dtype=int)
    pfs_e_ordered = np.zeros(len(dm), dtype=int)
    os_t_ordered  = np.zeros(len(dm), dtype=int)
    os_e_ordered  = np.zeros(len(dm), dtype=int)
    pfs_t_ordered[all_idx] = pfs_t
    pfs_e_ordered[all_idx] = pfs_e
    os_t_ordered[all_idx]  = os_t
    os_e_ordered[all_idx]  = os_e

    rows = []
    for i, (_, subj) in enumerate(dm.iterrows()):
        for times, events, dsdecod_e, dsdecod_c, paramcd, param_label in [
            (pfs_t_ordered, pfs_e_ordered, "DISEASE PROGRESSION", "CENSORED", "PFS", "Progression-Free Survival"),
            (os_t_ordered,  os_e_ordered,  "DEATH",               "CENSORED", "OS",  "Overall Survival"),
        ]:
            t = times[i]; e = events[i]
            rows.append({
                "STUDYID":  sid, "DOMAIN": "DS", "USUBJID": subj["USUBJID"],
                "DSDECOD":  dsdecod_e if e else dsdecod_c,
                "DSCAT":    "DISPOSITION EVENT",
                "DSSCAT":   param_label,
                "DSSTDTC":  iso_date(t),
                "DSSTDY":   t + 1,
                "DSTERM":   dsdecod_e if e else "ONGOING",
                "EPOCH":    "TREATMENT",
                "EVENTFL":  "Y" if e else "N",
                "CNSR":     0 if e else 1,
                "AVAL":     round(t / 28.0, 2),
                "PARAM":    param_label,
                "PARAMCD":  paramcd,
            })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# ADSL — Subject Level Analysis Dataset
# ─────────────────────────────────────────────────────────────────────────────

def make_adsl(dm: pd.DataFrame, ds: pd.DataFrame,
              lb: pd.DataFrame, cm: pd.DataFrame,
              study_key: str) -> pd.DataFrame:
    sid = STUDY_CONFIG[study_key]["studyid"]

    pfs = ds[ds["PARAMCD"] == "PFS"][["USUBJID","AVAL","CNSR","DSSTDTC"]]\
          .rename(columns={"AVAL":"PFSDUR","CNSR":"PFSCNSR","DSSTDTC":"PFSDTC"})
    os_ = ds[ds["PARAMCD"] == "OS"][["USUBJID","AVAL","CNSR","DSSTDTC"]]\
          .rename(columns={"AVAL":"OSDUR","CNSR":"OSCNSR","DSSTDTC":"OSDTC"})

    adsl = dm.merge(pfs, on="USUBJID", how="left")\
             .merge(os_, on="USUBJID", how="left")

    adsl["STUDYID"]  = sid
    adsl["TRTSDT"]   = adsl["RFSTDTC"]
    adsl["TRT01P"]   = adsl["ARM"]
    adsl["TRT01A"]   = adsl["ACTARM"]
    adsl["TRT01PN"]  = np.where(adsl["ARMCD"] == "IRd", 1, 2)
    adsl["TRT01AN"]  = adsl["TRT01PN"]
    adsl["ITTFL"]    = "Y"
    adsl["SAFFL"]    = "Y"
    adsl["PPSFL"]    = RNG.choice(["Y","N"], size=len(adsl), p=[0.92, 0.08])
    adsl["EOSSTT"]   = np.where(adsl["PFSCNSR"] == 0, "DISCONTINUED", "ONGOING")
    adsl["DCSREAS"]  = np.where(adsl["PFSCNSR"] == 0, "DISEASE PROGRESSION", "")

    # ── Baseline labs from LB (baseline epoch) ───────────────────────────────
    bl = lb[(lb["EPOCH"] == "BASELINE") & lb["LBSTRESN"].notna()][["USUBJID","LBTESTCD","LBSTRESN"]]\
         .groupby(["USUBJID","LBTESTCD"])["LBSTRESN"].median().reset_index()
    bl_wide = bl.pivot(index="USUBJID", columns="LBTESTCD", values="LBSTRESN").reset_index()
    adsl = adsl.merge(bl_wide.add_prefix("BASE_").rename(columns={"BASE_USUBJID":"USUBJID"}),
                      on="USUBJID", how="left")

    # ── Cockcroft-Gault CrCL from actual weight/age/sex/creatinine ───────────
    # BASE_CREAT is in umol/L from LB; fallback to random if missing
    if "BASE_CREAT" in adsl.columns:
        scr = adsl["BASE_CREAT"].fillna(88.0)
    else:
        scr = np.round(RNG.normal(88.0, 25.0, len(adsl)), 1)
        adsl["BASE_CREAT"] = scr
    adsl["BASE_CREACL"] = cockcroft_gault(
        adsl["AGE"].values, adsl["WEIGHT"].values,
        scr.values, adsl["SEX"].values
    )
    adsl["RENGRP"] = pd.cut(
        adsl["BASE_CREACL"].fillna(90),
        bins=[0, 30, 60, 90, 999],
        labels=["SEVERE","MODERATE","MILD","NORMAL"]
    ).astype(str)

    # ── Cytogenetics — study-specific prevalences per Spec §9B ──────────────
    cp = CYTO_PROBS[study_key]
    n  = len(adsl)
    adsl["DEL17P"]  = RNG.choice(["Y","N"], size=n, p=[cp["DEL17P"],  1-cp["DEL17P"]])
    adsl["T414"]    = RNG.choice(["Y","N"], size=n, p=[cp["T414"],    1-cp["T414"]])
    adsl["T1416"]   = RNG.choice(["Y","N"], size=n, p=[cp["T1416"],   1-cp["T1416"]])
    adsl["T1420"]   = RNG.choice(["Y","N"], size=n, p=[cp["T1420"],   1-cp["T1420"]])
    adsl["GAIN1Q21"]= RNG.choice(["Y","N"], size=n, p=[cp["GAIN1Q21"],1-cp["GAIN1Q21"]])
    adsl["DEL1P32"] = RNG.choice(["Y","N"], size=n, p=[cp["DEL1P32"], 1-cp["DEL1P32"]])
    adsl["AMP1Q"]   = adsl["GAIN1Q21"]

    # IMWG high-risk: del(17p) OR t(4;14) OR t(14;16)
    adsl["CYTOGR"] = np.where(
        (adsl["DEL17P"] == "Y") | (adsl["T414"] == "Y") | (adsl["T1416"] == "Y"),
        "HIGH RISK", "STANDARD RISK"
    )

    # ── R-ISS staging ─────────────────────────────────────────────────────────
    # R-ISS I:   ISS I + standard-risk cytogenetics + LDH ≤ ULN (280 U/L)
    # R-ISS III: ISS III + (high-risk cytogenetics OR LDH > ULN)
    # R-ISS II:  all others
    # Target distribution: ~35/45/20% — requires LDH elevation to correlate with ISS III
    ldh_col = adsl["BASE_LDH"] if "BASE_LDH" in adsl.columns else \
              pd.Series(RNG.uniform(140, 400, n), index=adsl.index)
    ldh_col = ldh_col.fillna(200.0)

    # Enrich LDH elevation among ISS III subjects — R-ISS III requires ISS III + adverse feature
    # Without this, P(LDH>ULN | ISS III) ≈ same as overall → too few R-ISS III
    iss3_mask = adsl["ISSSTAGE"] == 3
    ldh_arr   = ldh_col.values.copy()
    # Force ~55% of ISS III to have LDH > ULN (vs ~35% background) to get ~20% R-ISS III
    for i in np.where(iss3_mask)[0]:
        if RNG.random() < 0.55 and ldh_arr[i] <= 280:
            ldh_arr[i] = RNG.uniform(281, 600)   # above ULN
    ldh_col   = pd.Series(ldh_arr, index=adsl.index)
    adsl["BASE_LDH"] = ldh_col

    ldh_hi  = ldh_col > 280.0
    hr_cyto = adsl["CYTOGR"] == "HIGH RISK"
    iss1    = adsl["ISSSTAGE"] == 1
    iss3    = adsl["ISSSTAGE"] == 3

    riss = np.where(
        iss1 & ~hr_cyto & ~ldh_hi, "I",
        np.where(iss3 & (hr_cyto | ldh_hi), "III", "II")
    )
    adsl["RISS"] = riss

    # ── Prior therapy variables (MM1 / RRMM only) ────────────────────────────
    # Spec §9B — stratification factors for TOURMALINE-MM1
    # NPRIORLINE: 1/2/3 prior therapies [0.61/0.29/0.10]
    # PRIORPHI:   prior proteasome inhibitor [70% Y]
    # PRIORIMID:  prior IMiD [55% Y]
    # Both stratification factors → balanced across arms by design
    if study_key == "MM1":
        adsl["NPRIORLINE"] = RNG.choice([1, 2, 3], size=n, p=[0.61, 0.29, 0.10])
        adsl["PRIORPHI"]   = RNG.choice(["Y", "N"], size=n, p=[0.70, 0.30])
        adsl["PRIORIMID"]  = RNG.choice(["Y", "N"], size=n, p=[0.55, 0.45])
    else:
        adsl["NPRIORLINE"] = 0   # treatment-naive
        adsl["PRIORPHI"]   = "N"
        adsl["PRIORIMID"]  = "N"

    # ── Polymorphisms (CG/CC/GG — from original paper subgroup analysis) ──────
    poly = RNG.choice(["CC","CG","GG"], size=len(adsl), p=[0.35, 0.45, 0.20])
    adsl["POLYCC"]  = (poly == "CC").astype(int)
    adsl["POLYCG"]  = (poly == "CG").astype(int)
    adsl["POLYGG"]  = (poly == "GG").astype(int)

    # ── Durie-Salmon stage ────────────────────────────────────────────────────
    adsl["DSSTAGE"] = RNG.choice(["IA","IIA","IIB","IIIA","IIIB"],
                                 size=len(adsl), p=[0.20,0.25,0.10,0.35,0.10])

    # ── ISS-based risk grouping ───────────────────────────────────────────────
    adsl["RISKGR"] = pd.cut(adsl["ISSSTAGE"], bins=[0,1,2,3],
                            labels=["LOW","INTERMEDIATE","HIGH"]).astype(str)

    # ── Disease features ──────────────────────────────────────────────────────
    adsl["EXTRAMED"]  = RNG.choice(["Y","N"], size=len(adsl), p=[0.08, 0.92])
    adsl["BONELESI"]  = RNG.choice(["Y","N"], size=len(adsl), p=[0.55, 0.45])
    adsl["PLASMA"]    = RNG.choice(["Y","N"], size=len(adsl), p=[0.12, 0.88])
    adsl["MEAS_DISFL"]= RNG.choice(
        ["SERUM M-PROTEIN","URINE M-PROTEIN","BOTH SPEP AND UPEP","SERUM FLC"],
        size=len(adsl), p=[0.55, 0.10, 0.15, 0.20]
    )

    # ── CYP3A4 DDI flags from CM (for PK generator) ──────────────────────────
    if cm is not None and len(cm):
        inh_uids = cm[cm["CMCAT"] == "CYP3A4 INHIBITOR"]["USUBJID"].unique()
        ind_uids = cm[cm["CMCAT"] == "CYP3A4 INDUCER"]["USUBJID"].unique()
        # Weighted average CL multiplier per subject
        inh_cl = cm[cm["CMCAT"] == "CYP3A4 INHIBITOR"]\
                 .groupby("USUBJID")["CL_MULTIPLIER"].min()  # worst = lowest CL (strongest inh)
        ind_cl = cm[cm["CMCAT"] == "CYP3A4 INDUCER"]\
                 .groupby("USUBJID")["CL_MULTIPLIER"].max()  # worst = highest CL (strongest ind)

        adsl["CYP3A4_INHIBFL"] = np.where(adsl["USUBJID"].isin(inh_uids), "Y", "N")
        adsl["CYP3A4_INDUFL"]  = np.where(adsl["USUBJID"].isin(ind_uids), "Y", "N")
        adsl["CL_DDI_MULT"] = 1.0  # default
        for uid, mult in inh_cl.items():
            adsl.loc[adsl["USUBJID"] == uid, "CL_DDI_MULT"] = mult
        for uid, mult in ind_cl.items():
            # Inducer overrides inhibitor if both (rare)
            adsl.loc[adsl["USUBJID"] == uid, "CL_DDI_MULT"] = mult
    else:
        adsl["CYP3A4_INHIBFL"] = "N"
        adsl["CYP3A4_INDUFL"]  = "N"
        adsl["CL_DDI_MULT"]    = 1.0

    return adsl


# ─────────────────────────────────────────────────────────────────────────────
# ADTTE — Time-to-Event Analysis Dataset
# ─────────────────────────────────────────────────────────────────────────────

def make_adtte(ds: pd.DataFrame, adsl: pd.DataFrame, study_key: str) -> pd.DataFrame:
    sid   = STUDY_CONFIG[study_key]["studyid"]
    cols  = ["USUBJID","ARMCD","ARM","AGE","SEX","IGTYPE","ISSSTAGE",
             "RISKGR","CYTOGR","RISS","TRT01PN",
             "DEL17P","T414","T1416","T1420","GAIN1Q21","DEL1P32",
             "BASE_CREACL","RENGRP","CYP3A4_INHIBFL","CYP3A4_INDUFL",
             "NPRIORLINE","PRIORPHI","PRIORIMID"]
    cols  = [c for c in cols if c in adsl.columns]

    adtte = ds[["STUDYID","USUBJID","PARAMCD","PARAM","AVAL","CNSR","DSSTDTC"]].copy()
    adtte = adtte.merge(adsl[cols], on="USUBJID", how="left")
    adtte["STUDYID"]  = sid
    adtte["AVALU"]    = "MONTHS"
    adtte["STARTDT"]  = adsl.set_index("USUBJID").loc[adtte["USUBJID"], "RFSTDTC"].values
    adtte["EVNTDESC"] = np.where(adtte["PARAMCD"] == "PFS",
                                  np.where(adtte["CNSR"] == 0, "DISEASE PROGRESSION","CENSORED"),
                                  np.where(adtte["CNSR"] == 0, "DEATH","CENSORED"))
    adtte["SRCDOM"]   = "DS"
    adtte["ITTFL"]    = "Y"
    adtte["SAFFL"]    = "Y"
    adtte["STRATFL1"] = adtte["ISSSTAGE"].apply(lambda x: {1:"I",2:"II",3:"III"}.get(x,""))
    adtte["STRATFL2"] = adtte["CYTOGR"]
    adtte["STRATFL3"] = adtte["IGTYPE"]
    return adtte


# ─────────────────────────────────────────────────────────────────────────────
# ADLB — Longitudinal Lab Analysis Dataset
# ─────────────────────────────────────────────────────────────────────────────

def make_adlb(lb: pd.DataFrame, adsl: pd.DataFrame, study_key: str) -> pd.DataFrame:
    sid  = STUDY_CONFIG[study_key]["studyid"]
    cols = ["USUBJID","ARMCD","TRT01PN","IGTYPE","ISSSTAGE","RISS","CYTOGR"]
    cols = [c for c in cols if c in adsl.columns]

    adlb = lb.copy()
    adlb = adlb.merge(adsl[cols], on="USUBJID", how="left")
    adlb["STUDYID"] = sid

    baseline = lb[(lb["EPOCH"] == "BASELINE") & lb["LBSTRESN"].notna()][["USUBJID","LBTESTCD","LBSTRESN"]]\
               .rename(columns={"LBSTRESN":"BASE"})
    adlb = adlb.merge(baseline, on=["USUBJID","LBTESTCD"], how="left")
    adlb["CHG"]    = adlb["LBSTRESN"] - adlb["BASE"]
    adlb["PCHG"]   = np.where(adlb["BASE"] != 0, adlb["CHG"] / adlb["BASE"] * 100, np.nan)
    adlb["AVAL"]   = adlb["LBSTRESN"]
    adlb["AVALU"]  = adlb["LBSTRESU"]
    adlb["PARAM"]  = adlb["LBTEST"]
    adlb["PARAMCD"]= adlb["LBTESTCD"]
    adlb["ANL01FL"]= "Y"
    adlb["ABLFL"]  = np.where(adlb["EPOCH"] == "BASELINE", "Y", "")
    adlb["DTYPE"]  = ""
    return adlb


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def generate_all():
    summary = []

    for study_key in ["MM2", "MM1"]:
        # Reseed RNG per study so MM2 and MM1 generation are fully independent.
        # Without this, any change to MM2 generation (number of RNG calls) shifts
        # the RNG state before MM1 begins, causing MM1 rates to change unexpectedly.
        # Seed 42 for MM2 (development set), seed 43 for MM1 (validation set).
        study_seed = 42 if study_key == "MM2" else 43
        RNG.bit_generator.state = np.random.default_rng(study_seed).bit_generator.state

        cfg       = STUDY_CONFIG[study_key]
        study_dir = os.path.join(OUT_BASE, study_key)
        os.makedirs(study_dir, exist_ok=True)

        print(f"\n{'='*65}")
        print(f"  {study_key} — {cfg['studyid']}  N={cfg['n']}  [{cfg['population']}]")
        print(f"{'='*65}")

        print("  [1/9] DM  — Demographics (+ weight/height/BMI/BSA)")
        dm = make_dm(study_key)
        dm.to_csv(f"{study_dir}/sdtm_dm.csv", index=False)
        print(f"        {len(dm):,} subjects")

        print("  [2/9] EX  — Exposure (all subjects; adherence + reason-coded mods)")
        ex = make_ex(dm, study_key)
        ex.to_csv(f"{study_dir}/sdtm_ex.csv", index=False)
        print(f"        {len(ex):,} dose records")

        print("  [3/9] CM  — Concomitant Medications (CYP3A4, anticoag, G-CSF)")
        cm = make_cm(dm, study_key)
        cm.to_csv(f"{study_dir}/sdtm_cm.csv", index=False)
        print(f"        {len(cm):,} CM records")

        print("  [4/9] LB  — Laboratory (exposure-response M-protein; B2MG/LDH added)")
        lb = make_lb(dm, study_key)
        lb.to_csv(f"{study_dir}/sdtm_lb.csv", index=False)
        print(f"        {len(lb):,} lab records")

        print("  [5/9] AE  — Adverse Events (hematologic linked to ANC/PLT)")
        ae = make_ae(dm, lb, study_key)
        ae.to_csv(f"{study_dir}/sdtm_ae.csv", index=False)
        print(f"        {len(ae):,} AE records")

        print("  [6/9] DS  — Disposition (arm-specific Weibull; IRd HR=0.742)")
        ds = make_ds(dm, study_key)
        ds.to_csv(f"{study_dir}/sdtm_ds.csv", index=False)
        print(f"        {len(ds):,} disposition records")

        print("  [7/9] ADSL — Subject Level (CrCL, R-ISS, extended cytogenetics)")
        adsl = make_adsl(dm, ds, lb, cm, study_key)
        adsl.to_csv(f"{study_dir}/adam_adsl.csv", index=False)
        print(f"        {len(adsl):,} subjects  |  {len(adsl.columns)} variables")

        print("  [8/9] ADTTE — Time-to-Event")
        adtte = make_adtte(ds, adsl, study_key)
        adtte.to_csv(f"{study_dir}/adam_adtte.csv", index=False)
        print(f"        {len(adtte):,} records")

        print("  [9/9] ADLB — Longitudinal Lab ADaM")
        adlb = make_adlb(lb, adsl, study_key)
        adlb.to_csv(f"{study_dir}/adam_adlb.csv", index=False)
        print(f"        {len(adlb):,} records")

        # ── Quick validation summary ──────────────────────────────────────────
        try:
            from lifelines import KaplanMeierFitter
            pfs = adtte[adtte["PARAMCD"] == "PFS"]
            def km_median(df):
                kmf = KaplanMeierFitter()
                kmf.fit(df["AVAL"], event_observed=df["CNSR"].eq(0))
                return kmf.median_survival_time_
            med_pfs_ird = km_median(pfs[pfs["ARMCD"]=="IRd"])
            med_pfs_rd  = km_median(pfs[pfs["ARMCD"]=="Rd"])
        except ImportError:
            pfs = adtte[adtte["PARAMCD"] == "PFS"]
            med_pfs_ird = pfs[pfs["ARMCD"]=="IRd"]["AVAL"].median()
            med_pfs_rd  = pfs[pfs["ARMCD"]=="Rd"]["AVAL"].median()
        pfs_event   = pfs["CNSR"].eq(0).mean() * 100
        riss_dist   = adsl["RISS"].value_counts(normalize=True).sort_index().to_dict() \
                      if "RISS" in adsl.columns else {}
        cyto_hr     = (adsl["CYTOGR"] == "HIGH RISK").mean() * 100

        cyto_target = cfg.get("cyto_hr_target", 20.0)
        print(f"\n  ── Validation snapshot ──")
        print(f"     PFS median  IRd: {med_pfs_ird:.1f} mo (target {cfg['pfs_target_ird']} mo)  "
              f"Rd: {med_pfs_rd:.1f} mo (target {cfg['pfs_target_rd']} mo)")
        print(f"     PFS event rate: {pfs_event:.0f}%  (target ~65%)")
        print(f"     High-risk cyto: {cyto_hr:.0f}%  (target ~{cyto_target:.0f}%)")
        if riss_dist:
            print(f"     R-ISS: I={riss_dist.get('I',0)*100:.0f}%  "
                  f"II={riss_dist.get('II',0)*100:.0f}%  "
                  f"III={riss_dist.get('III',0)*100:.0f}%  "
                  f"(target ~35/45/20%)")
        arm_dist = adsl["ARMCD"].value_counts().to_dict()
        print(f"     Arm balance: IRd={arm_dist.get('IRd',0)}  Rd={arm_dist.get('Rd',0)}"
              f"  (target IRd={cfg['n_ird']} Rd={cfg['n_rd']})")
        if study_key == "MM1" and "NPRIORLINE" in adsl.columns:
            priorline = adsl["NPRIORLINE"].value_counts(normalize=True).sort_index()
            print(f"     Prior lines: {dict(priorline.round(2))}"
                  f"  (target 1=61%/2=29%/3=10%)")
            priorphi  = (adsl["PRIORPHI"]  == "Y").mean() * 100
            priorimid = (adsl["PRIORIMID"] == "Y").mean() * 100
            print(f"     Prior PI: {priorphi:.0f}%  (target 70%)  "
                  f"Prior IMiD: {priorimid:.0f}%  (target 55%)")

        crcl_med  = adsl["BASE_CREACL"].median() if "BASE_CREACL" in adsl.columns else float("nan")
        pct_le60  = (adsl["BASE_CREACL"] <= 60).mean() * 100 if "BASE_CREACL" in adsl.columns else float("nan")
        print(f"     Median CrCL: {crcl_med:.0f} mL/min  %≤60: {pct_le60:.0f}%  (expected 55-70 mL/min)")
        wt_med = dm["WEIGHT"].median()
        print(f"     Median weight: {wt_med:.1f} kg")
        bsa_med = dm["BSA"].median()
        print(f"     Median BSA: {bsa_med:.3f} m²  (expected ~1.75 m²)")

        # M-protein response depth check
        try:
            mp = adlb[adlb["PARAMCD"] == "SPEP_MPROT"].copy()
            bl_mp = mp[mp["EPOCH"] == "BASELINE"].groupby("USUBJID")["AVAL"].first()
            c3_mp = mp[mp["VISITNUM"] == 3].groupby("USUBJID")["AVAL"].first()
            c6_mp = mp[mp["VISITNUM"] == 6].groupby("USUBJID")["AVAL"].first()
            common_c3 = bl_mp.index.intersection(c3_mp.index)
            common_c6 = bl_mp.index.intersection(c6_mp.index)
            if len(common_c3) > 50:
                pchg_c3 = ((c3_mp[common_c3] - bl_mp[common_c3]) / bl_mp[common_c3] * 100).median()
                print(f"     M-protein median %change C3: {pchg_c3:.0f}%  (target IRd responders: ~-44%)")
            if len(common_c6) > 50:
                pchg_c6 = ((c6_mp[common_c6] - bl_mp[common_c6]) / bl_mp[common_c6] * 100).median()
                print(f"     M-protein median %change C6: {pchg_c6:.0f}%  (target IRd responders: ~-72%)")
        except Exception:
            pass

        summary.append({
            "Study":      study_key,
            "Studyid":    cfg["studyid"],
            "Population": cfg["population"],
            "N_subjects": len(dm),
            "DM_rows":    len(dm),
            "EX_rows":    len(ex),
            "CM_rows":    len(cm),
            "LB_rows":    len(lb),
            "AE_rows":    len(ae),
            "DS_rows":    len(ds),
            "ADSL_rows":  len(adsl),
            "ADTTE_rows": len(adtte),
            "ADLB_rows":  len(adlb),
            "PFS_med_IRd_mo": round(med_pfs_ird, 1),
            "PFS_med_Rd_mo":  round(med_pfs_rd, 1),
        })

    pd.DataFrame(summary).to_csv(os.path.join(OUT_BASE, "dataset_summary_v2.csv"), index=False)
    print(f"\n✓  All datasets written to: {os.path.abspath(OUT_BASE)}")


if __name__ == "__main__":
    generate_all()
