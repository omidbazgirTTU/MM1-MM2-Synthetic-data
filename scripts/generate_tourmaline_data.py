"""
Synthetic TOURMALINE MM1 / MM2 Trial Data Generator
=====================================================
Generates SDTM and ADaM datasets compatible with CDISC standards, based on
the paper: "Joint AI-driven event prediction and longitudinal modeling in
newly diagnosed and relapsed multiple myeloma" (Hussain et al., 2024)

SDTM Domains: DM, EX, LB, AE, DS
ADaM Datasets: ADSL, ADTTE, ADLB
"""

import numpy as np
import pandas as pd
from scipy.stats import weibull_min
import warnings
warnings.filterwarnings("ignore")

RNG = np.random.default_rng(42)

# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTS & CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

STUDY_CONFIG = {
    "MM2": {
        "studyid": "TOURMALINE-MM2",
        "protocol": "NCT01850524",
        "population": "NDMM",
        "n": 703,
        "median_age": 73,
        "age_range": (48, 90),
        "pct_female": 0.499,
        "median_dx_months": 1.11,
    },
    "MM1": {
        "studyid": "TOURMALINE-MM1",
        "protocol": "NCT01564537",
        "population": "RRMM",
        "n": 720,
        "median_age": 66,
        "age_range": (30, 91),
        "pct_female": 0.433,
        "median_dx_months": 42.8,
    },
}

# Ig-type distributions (from Table 1)
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
            "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER": 0.006, "OTHER": 0.010,
            "NOT REPORTED": 0.024},
}

ISS_DIST = {
    "MM2": {1: 0.461, 2: 0.374, 3: 0.164},
    "MM1": {1: 0.636, 2: 0.244, 3: 0.120},
}

# Drug dosing: standard doses (mg), 28-day cycles
DRUG_DOSES = {
    "IXAZOMIB":     {"dose": 4.0,   "unit": "mg",  "route": "ORAL", "freq": "WEEKLY (DAYS 1,8,15)"},
    "LENALIDOMIDE": {"dose": 25.0,  "unit": "mg",  "route": "ORAL", "freq": "DAILY (DAYS 1-21)"},
    "DEXAMETHASONE":{"dose": 40.0,  "unit": "mg",  "route": "ORAL", "freq": "WEEKLY (DAYS 1,8,15,22)"},
    "PLACEBO":      {"dose": 0.0,   "unit": "mg",  "route": "ORAL", "freq": "WEEKLY (DAYS 1,8,15)"},
}

# Lab normal ranges (for normalization per paper eq. 1)
LAB_NORMAL_RANGES = {
    "ALBUMIN":          (35,   50,   "g/L"),
    "ALP":              (40,  130,   "U/L"),
    "ALT":              (7,    56,   "U/L"),
    "AST":              (10,   40,   "U/L"),
    "BILI":             (3,    21,   "umol/L"),
    "BUN":              (2.5,  7.1,  "mmol/L"),
    "CA":               (2.2,  2.6,  "mmol/L"),
    "CL":               (98,  106,   "mmol/L"),
    "CA_CORR":          (2.2,  2.6,  "mmol/L"),
    "GLUC":             (3.9,  5.6,  "mmol/L"),
    "K":                (3.5,  5.0,  "mmol/L"),
    "MG":               (0.7,  1.0,  "mmol/L"),
    "PHOS":             (0.8,  1.5,  "mmol/L"),
    "CO2":              (22,   29,   "mmol/L"),
    "CREAT":            (44,  106,   "umol/L"),
    "GFR":              (60,  120,   "mL/min/1.73m2"),
    "HCT":              (36,   46,   "%"),
    "HGB":              (120, 160,   "g/L"),
    "LYMPH":            (1.0,  4.8,  "10^9/L"),
    "MONO":             (0.2,  1.0,  "10^9/L"),
    "NEUT":             (1.8,  7.7,  "10^9/L"),
    "PLT":              (150, 400,   "10^9/L"),
    "WBC":              (4.5, 11.0,  "10^9/L"),
    "LDH":              (140, 280,   "U/L"),
    "PROT":             (60,   80,   "g/L"),
    "GLOB":             (20,   35,   "g/L"),
    "NA":               (136, 145,   "mmol/L"),
    "SPEP_GAMMA":       (7,    16,   "g/L"),
    "SPEP_KAPPA":       (3.3, 19.4,  "mg/L"),
    "KAPPA_LAMBDA":     (0.26, 1.65, "ratio"),
    "SPEP_LAMBDA":      (5.7, 26.3,  "mg/L"),
    "SPEP_MPROT":       (0,    5,    "g/L"),
    "IGA":              (0.7,  4.0,  "g/L"),
    "IGG":              (7.0, 16.0,  "g/L"),
    "IGM":              (0.4,  2.3,  "g/L"),
    "URINE_ALB":        (0,    30,   "%"),
    "UPEP_MPROT":       (0,    80,   "mg/day"),
    "URATE":            (155, 428,   "umol/L"),
}

# Adverse event definitions (from paper + CTCAE grades)
AE_DEFS = [
    ("Acute Renal Failure",     "RENAL AND URINARY DISORDERS",      "C-index=0.62"),
    ("Cardiac Arrhythmias",     "CARDIAC DISORDERS",                ""),
    ("Diarrhea",                "GASTROINTESTINAL DISORDERS",       ""),
    ("Heart Failure",           "CARDIAC DISORDERS",                ""),
    ("Hypotension",             "VASCULAR DISORDERS",               "C-index=0.66"),
    ("Liver Impairment",        "HEPATOBILIARY DISORDERS",          ""),
    ("Nausea",                  "GASTROINTESTINAL DISORDERS",       ""),
    ("Neutropenia",             "BLOOD AND LYMPHATIC SYSTEM DISORDERS", "C-index=0.59"),
    ("Peripheral Neuropathies", "NERVOUS SYSTEM DISORDERS",         ""),
    ("Rash",                    "SKIN AND SUBCUTANEOUS TISSUE DISORDERS", ""),
    ("Thrombocytopenia",        "BLOOD AND LYMPHATIC SYSTEM DISORDERS", "C-index=0.85"),
    ("Vomiting",                "GASTROINTESTINAL DISORDERS",       ""),
]

# AE grade thresholds (hematologic ≥3, non-hematologic ≥2)
HEMATOLOGIC_AE = {"Neutropenia", "Thrombocytopenia"}

# AE incidence rates per arm (approximate, per cycle)
AE_INCIDENCE = {
    "IRd": {"Acute Renal Failure": 0.04, "Cardiac Arrhythmias": 0.02, "Diarrhea": 0.15,
            "Heart Failure": 0.02, "Hypotension": 0.03, "Liver Impairment": 0.02,
            "Nausea": 0.12, "Neutropenia": 0.10, "Peripheral Neuropathies": 0.06,
            "Rash": 0.08, "Thrombocytopenia": 0.08, "Vomiting": 0.07},
    "Rd":  {"Acute Renal Failure": 0.03, "Cardiac Arrhythmias": 0.02, "Diarrhea": 0.12,
            "Heart Failure": 0.02, "Hypotension": 0.025, "Liver Impairment": 0.015,
            "Nausea": 0.09, "Neutropenia": 0.12, "Peripheral Neuropathies": 0.04,
            "Rash": 0.06, "Thrombocytopenia": 0.06, "Vomiting": 0.05},
}


# ─────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def weighted_choice(options: dict, n: int) -> np.ndarray:
    keys = list(options.keys())
    probs = np.array(list(options.values()))
    probs /= probs.sum()
    return RNG.choice(keys, size=n, p=probs)


def iso_date(days_from_ref: int, ref_year: int = 2015) -> str:
    import datetime
    ref = datetime.date(ref_year, 1, 1)
    d = ref + datetime.timedelta(days=int(days_from_ref))
    return d.strftime("%Y-%m-%d")


def days_to_dtc(days: int) -> str:
    return iso_date(days)


def simulate_survival(n: int, median_months: float, censor_rate: float = 0.35):
    """Weibull survival with right censoring."""
    scale = median_months * 28  # convert to days
    shape = 1.2
    event_times = weibull_min.rvs(shape, scale=scale, size=n, random_state=RNG.integers(9999))
    censor_times = RNG.exponential(scale * 2, size=n)
    observed_time = np.minimum(event_times, censor_times)
    event_flag = (event_times <= censor_times).astype(int)
    return observed_time.astype(int), event_flag


# ─────────────────────────────────────────────────────────────────────────────
# DM — Demographics
# ─────────────────────────────────────────────────────────────────────────────

def make_dm(study_key: str) -> pd.DataFrame:
    cfg = STUDY_CONFIG[study_key]
    n = cfg["n"]
    study_id = cfg["studyid"]

    usubjid = [f"{study_id}-{i+1:04d}" for i in range(n)]
    subjid  = [f"{i+1:04d}" for i in range(n)]

    ages = RNG.integers(cfg["age_range"][0], cfg["age_range"][1]+1, size=n)
    # bias toward median
    ages = np.clip(
        np.round(RNG.normal(cfg["median_age"], 8, n)).astype(int),
        cfg["age_range"][0], cfg["age_range"][1]
    )

    sex = RNG.choice(["F","M"], size=n, p=[cfg["pct_female"], 1-cfg["pct_female"]])
    race = weighted_choice(RACE_DIST[study_key], n)
    ig_type = weighted_choice(IG_TYPE_DIST[study_key], n)
    lchain = RNG.choice(["KAPPA","LAMBDA","BICLONAL"], size=n, p=[0.60, 0.35, 0.05])
    iss_probs = np.array(list(ISS_DIST[study_key].values()), dtype=float)
    iss_probs /= iss_probs.sum()
    iss = np.array([
        RNG.choice([1,2,3], p=iss_probs) for _ in range(n)
    ])
    ecog = RNG.choice([0,1,2], size=n, p=[0.35, 0.50, 0.15])

    dx_months = np.abs(RNG.exponential(cfg["median_dx_months"], size=n))

    # Treatment assignment: 50/50 randomization
    arm = RNG.choice(["IRd","Rd"], size=n, p=[0.504, 0.496])

    # Random study entry dates (2015–2018)
    rfstdtc = np.array([iso_date(RNG.integers(0, 3*365)) for _ in range(n)])

    dm = pd.DataFrame({
        "STUDYID":   study_id,
        "DOMAIN":    "DM",
        "USUBJID":   usubjid,
        "SUBJID":    subjid,
        "RFSTDTC":   rfstdtc,
        "RFENDTC":   "",            # filled later
        "SITEID":    RNG.choice([f"SITE{i:02d}" for i in range(1,16)], size=n),
        "AGE":       ages,
        "AGEU":      "YEARS",
        "SEX":       sex,
        "RACE":      race,
        "ETHNIC":    RNG.choice(["HISPANIC OR LATINO","NOT HISPANIC OR LATINO","NOT REPORTED"],
                                size=n, p=[0.05, 0.90, 0.05]),
        "COUNTRY":   RNG.choice(["USA","FRA","DEU","GBR","JPN","ITA","ESP","AUS"],
                                size=n, p=[0.35,0.12,0.10,0.10,0.09,0.08,0.08,0.08]),
        "ARMCD":     arm,
        "ARM":       np.where(arm=="IRd",
                              "Ixazomib + Lenalidomide + Dexamethasone",
                              "Placebo + Lenalidomide + Dexamethasone"),
        "ACTARMCD":  arm,
        "ACTARM":    np.where(arm=="IRd",
                              "Ixazomib + Lenalidomide + Dexamethasone",
                              "Placebo + Lenalidomide + Dexamethasone"),
        "DTHDTC":    "",
        "DTHFL":     "",
        # Supplemental
        "IGTYPE":    ig_type,
        "LCTYPE":    lchain,
        "ISSSTAGE":  iss,
        "ECOG":      ecog,
        "DXMONTHS":  np.round(dx_months, 2),
        "POPULATION": cfg["population"],
    })
    return dm


# ─────────────────────────────────────────────────────────────────────────────
# EX — Exposure (Dosing)
# ─────────────────────────────────────────────────────────────────────────────

def make_ex(dm: pd.DataFrame, study_key: str) -> pd.DataFrame:
    rows = []
    study_id = STUDY_CONFIG[study_key]["studyid"]
    max_cycles = 26  # ~2 years

    for _, subj in dm.iterrows():
        arm = subj["ARMCD"]
        rfst = subj["RFSTDTC"]

        # Determine actual number of cycles per patient (censored/event)
        n_cycles = RNG.integers(1, max_cycles + 1)

        drugs = ["LENALIDOMIDE", "DEXAMETHASONE"]
        drugs += ["IXAZOMIB"] if arm == "IRd" else ["PLACEBO"]

        for drug in drugs:
            info = DRUG_DOSES[drug]
            for cyc in range(1, n_cycles + 1):
                # Compute dosing days within cycle
                cycle_start_day = (cyc - 1) * 28
                if "WEEKLY" in info["freq"]:
                    dose_days = [cycle_start_day + d for d in [0, 7, 14]]
                    if "22" in info["freq"]:
                        dose_days.append(cycle_start_day + 21)
                else:  # DAILY days 1-21
                    dose_days = list(range(cycle_start_day, cycle_start_day + 21))

                for dd in dose_days:
                    # Occasional dose modifications
                    dose_mod = RNG.choice([1.0, 0.75, 0.5], p=[0.85, 0.10, 0.05])
                    actual_dose = round(info["dose"] * dose_mod, 1) if drug != "PLACEBO" else 0.0

                    rows.append({
                        "STUDYID":   study_id,
                        "DOMAIN":    "EX",
                        "USUBJID":   subj["USUBJID"],
                        "EXSEQ":     None,
                        "EXTRT":     drug,
                        "EXDOSE":    actual_dose,
                        "EXDOSU":    info["unit"],
                        "EXDOSFRM":  "TABLET",
                        "EXDOSFRQ":  info["freq"],
                        "EXROUTE":   info["route"],
                        "EXSTDTC":   iso_date(dd),
                        "EXENDTC":   iso_date(dd),
                        "VISITNUM":  cyc,
                        "VISIT":     f"CYCLE {cyc}",
                        "EPOCH":     "TREATMENT",
                        "EXDOSMOD":  "DOSE REDUCTION" if dose_mod < 1.0 else "",
                    })

    ex = pd.DataFrame(rows)
    ex["EXSEQ"] = ex.groupby("USUBJID").cumcount() + 1
    return ex


# ─────────────────────────────────────────────────────────────────────────────
# LB — Laboratory Tests (Longitudinal Biomarkers)
# ─────────────────────────────────────────────────────────────────────────────

def _simulate_lab_trajectory(test: str, ig_type: str, arm: str,
                              n_cycles: int, is_ndmm: bool) -> np.ndarray:
    """
    Simulate physiologically plausible longitudinal lab values.
    Disease biomarkers (M-protein, IgG, IgA) decrease with treatment;
    hemoglobin and GFR show treatment-dependent trends.
    """
    lo, hi, _ = LAB_NORMAL_RANGES[test]
    mid = (lo + hi) / 2
    rng_scale = (hi - lo) / 4

    # Baseline: abnormal for disease markers, near normal for others
    disease_markers = {"SPEP_MPROT", "IGA", "IGG", "IGM", "SPEP_KAPPA", "SPEP_LAMBDA",
                       "SPEP_GAMMA", "UPEP_MPROT"}
    progression_markers = {"CREAT", "BUN", "LDH", "CA", "CA_CORR"}
    anemia_markers = {"HGB", "HCT"}
    protective_markers = {"GFR", "ALBUMIN"}

    vals = np.zeros(n_cycles)

    if test in disease_markers:
        # High at baseline, decreasing with treatment (response)
        # IRd has better response than Rd (lower residual)
        base = hi * RNG.uniform(2.0, 5.0) if ig_type in {"IgG","IgA"} else hi * RNG.uniform(0.5, 2.0)
        resp_rate = 0.75 if arm == "IRd" else 0.65
        nadir = base * (1 - resp_rate)
        for t in range(n_cycles):
            decay = np.exp(-0.3 * t)
            relapse_noise = RNG.normal(0, base * 0.05)
            vals[t] = max(0, nadir + (base - nadir) * decay + relapse_noise)

    elif test in anemia_markers:
        # Low at baseline, mild recovery
        base = lo * RNG.uniform(0.70, 0.90)
        for t in range(n_cycles):
            vals[t] = base + (mid - base) * min(1, t / 12) + RNG.normal(0, rng_scale * 0.1)

    elif test in protective_markers:
        base = lo * RNG.uniform(0.85, 1.0)
        for t in range(n_cycles):
            vals[t] = base + RNG.normal(0, rng_scale * 0.1)

    elif test in progression_markers:
        base = mid * RNG.uniform(0.9, 1.5)
        for t in range(n_cycles):
            trend = RNG.normal(0, rng_scale * 0.05)
            vals[t] = max(lo * 0.5, base + trend * t)

    else:
        base = mid + RNG.normal(0, rng_scale * 0.5)
        for t in range(n_cycles):
            vals[t] = max(lo * 0.3, base + RNG.normal(0, rng_scale * 0.08))

    return np.clip(vals, lo * 0.1, hi * 8)


def make_lb(dm: pd.DataFrame, study_key: str, max_cycles: int = 24) -> pd.DataFrame:
    rows = []
    study_id = STUDY_CONFIG[study_key]["studyid"]
    is_ndmm = study_key == "MM2"

    tests = list(LAB_NORMAL_RANGES.keys())

    for _, subj in dm.iterrows():
        n_cycles = RNG.integers(2, max_cycles + 1)
        arm = subj["ARMCD"]
        ig = subj["IGTYPE"]

        for test in tests:
            _, _, unit = LAB_NORMAL_RANGES[test]
            traj = _simulate_lab_trajectory(test, ig, arm, n_cycles, is_ndmm)

            for cyc_idx, val in enumerate(traj):
                # Occasional missing (LOCF in original paper)
                if RNG.random() < 0.15:
                    continue

                rows.append({
                    "STUDYID":  study_id,
                    "DOMAIN":   "LB",
                    "USUBJID":  subj["USUBJID"],
                    "LBSEQ":    None,
                    "LBCAT":    _lb_cat(test),
                    "LBTESTCD": test,
                    "LBTEST":   _lb_label(test),
                    "LBORRES":  f"{val:.3f}",
                    "LBORRESU": unit,
                    "LBSTRESC": f"{val:.3f}",
                    "LBSTRESN": round(float(val), 3),
                    "LBSTRESU": unit,
                    "LBNRLO":   LAB_NORMAL_RANGES[test][0],
                    "LBNRHI":   LAB_NORMAL_RANGES[test][1],
                    "LBNRIND":  _nrind(val, test),
                    "VISITNUM": cyc_idx + 1,
                    "VISIT":    f"CYCLE {cyc_idx+1}" if cyc_idx > 0 else "BASELINE",
                    "LBDTC":    iso_date(cyc_idx * 28),
                    "LBDY":     cyc_idx * 28 + 1,
                    "EPOCH":    "TREATMENT" if cyc_idx > 0 else "BASELINE",
                })

    lb = pd.DataFrame(rows)
    lb["LBSEQ"] = lb.groupby("USUBJID").cumcount() + 1
    return lb


def _lb_cat(test: str) -> str:
    heme = {"HGB","HCT","LYMPH","MONO","NEUT","PLT","WBC"}
    immuno = {"IGA","IGG","IGM","SPEP_KAPPA","SPEP_LAMBDA","SPEP_GAMMA",
              "SPEP_MPROT","UPEP_MPROT","KAPPA_LAMBDA"}
    chem = {"ALBUMIN","ALP","ALT","AST","BILI","BUN","CA","CL","CA_CORR",
            "GLUC","K","MG","PHOS","CO2","CREAT","GFR","LDH","PROT","GLOB",
            "NA","URATE","URINE_ALB"}
    if test in heme:   return "HEMATOLOGY"
    if test in immuno: return "SERUM IMMUNOGLOBULINS"
    return "CHEMISTRY"


def _lb_label(test: str) -> str:
    labels = {
        "ALBUMIN": "Albumin", "ALP": "Alkaline Phosphatase",
        "ALT": "Alanine Aminotransferase", "AST": "Aspartate Aminotransferase",
        "BILI": "Bilirubin", "BUN": "Blood Urea Nitrogen",
        "CA": "Calcium", "CL": "Chloride", "CA_CORR": "Corrected Calcium",
        "GLUC": "Glucose", "K": "Potassium", "MG": "Magnesium",
        "PHOS": "Phosphate", "CO2": "Carbon Dioxide",
        "CREAT": "Creatinine", "GFR": "Glomerular Filtration Rate",
        "HCT": "Hematocrit", "HGB": "Hemoglobin",
        "LYMPH": "Lymphocytes", "MONO": "Monocytes",
        "NEUT": "Neutrophils", "PLT": "Platelets", "WBC": "Leukocytes",
        "LDH": "Lactate Dehydrogenase",
        "PROT": "Protein", "GLOB": "Serum Globulin", "NA": "Sodium",
        "SPEP_GAMMA": "SPEP Gamma Globulin",
        "SPEP_KAPPA": "Free SPEP Kappa Light Chain",
        "KAPPA_LAMBDA": "Free SPEP Kappa/Lambda Ratio",
        "SPEP_LAMBDA": "Free SPEP Lambda Light Chain",
        "SPEP_MPROT": "SPEP Monoclonal Protein",
        "IGA": "Immunoglobulin A", "IGG": "Immunoglobulin G",
        "IGM": "Immunoglobulin M",
        "URINE_ALB": "Urine Albumin", "UPEP_MPROT": "UPEP Monoclonal Protein",
        "URATE": "Urate",
    }
    return labels.get(test, test)


def _nrind(val: float, test: str) -> str:
    lo, hi, _ = LAB_NORMAL_RANGES[test]
    if val < lo:  return "LOW"
    if val > hi:  return "HIGH"
    return "NORMAL"


# ─────────────────────────────────────────────────────────────────────────────
# AE — Adverse Events
# ─────────────────────────────────────────────────────────────────────────────

def make_ae(dm: pd.DataFrame, study_key: str) -> pd.DataFrame:
    rows = []
    study_id = STUDY_CONFIG[study_key]["studyid"]

    for _, subj in dm.iterrows():
        arm = subj["ARMCD"]
        n_cycles = RNG.integers(2, 27)

        for ae_name, soc, _ in AE_DEFS:
            # Incidence probability over follow-up
            per_cycle_prob = AE_INCIDENCE[arm][ae_name]
            occurred = RNG.random(n_cycles) < per_cycle_prob

            for cyc_idx in np.where(occurred)[0]:
                # CTCAE grade
                is_heme = ae_name in HEMATOLOGIC_AE
                if is_heme:
                    grade = RNG.choice([3, 4, 5], p=[0.55, 0.40, 0.05])
                else:
                    grade = RNG.choice([2, 3, 4], p=[0.60, 0.35, 0.05])

                onset_day = cyc_idx * 28 + RNG.integers(1, 28)
                duration  = RNG.integers(3, 45)

                rows.append({
                    "STUDYID":   study_id,
                    "DOMAIN":    "AE",
                    "USUBJID":   subj["USUBJID"],
                    "AESEQ":     None,
                    "AETERM":    ae_name.upper(),
                    "AEDECOD":   ae_name,
                    "AEBODSYS":  soc,
                    "AESOC":     soc,
                    "AESEV":     {1:"MILD",2:"MODERATE",3:"SEVERE",
                                  4:"LIFE THREATENING",5:"FATAL"}[grade],
                    "AETOXGR":   str(grade),
                    "AESER":     "Y" if grade >= 3 else "N",
                    "AEREL":     RNG.choice(["RELATED","NOT RELATED","POSSIBLY RELATED"],
                                            p=[0.55, 0.25, 0.20]),
                    "AEACN":     RNG.choice(["DOSE REDUCTION","DRUG WITHDRAWN",
                                             "DOSE NOT CHANGED","DRUG INTERRUPTED"],
                                            p=[0.20, 0.05, 0.55, 0.20]),
                    "AESTDTC":   iso_date(onset_day),
                    "AEENDTC":   iso_date(onset_day + duration),
                    "AEDY":      onset_day + 1,
                    "AEENDY":    onset_day + duration + 1,
                    "AEOUT":     RNG.choice(["RECOVERED/RESOLVED",
                                             "RECOVERING/RESOLVING",
                                             "NOT RECOVERED/NOT RESOLVED"],
                                            p=[0.70, 0.20, 0.10]),
                    "EPOCH":     "TREATMENT",
                    "VISITNUM":  cyc_idx + 1,
                    "VISIT":     f"CYCLE {cyc_idx+1}",
                })

    ae = pd.DataFrame(rows)
    ae["AESEQ"] = ae.groupby("USUBJID").cumcount() + 1
    return ae


# ─────────────────────────────────────────────────────────────────────────────
# DS — Disposition
# ─────────────────────────────────────────────────────────────────────────────

def make_ds(dm: pd.DataFrame, study_key: str) -> pd.DataFrame:
    rows = []
    study_id = STUDY_CONFIG[study_key]["studyid"]
    is_ndmm = (study_key == "MM2")

    # Survival parameters (median PFS ~18mo NDMM, ~12mo RRMM)
    median_pfs = 18 if is_ndmm else 12
    median_os  = 40 if is_ndmm else 30

    pfs_times, pfs_events = simulate_survival(len(dm), median_pfs)
    os_times,  os_events  = simulate_survival(len(dm), median_os)

    for i, (_, subj) in enumerate(dm.iterrows()):
        # PFS event
        rows.append({
            "STUDYID":  study_id,
            "DOMAIN":   "DS",
            "USUBJID":  subj["USUBJID"],
            "DSSEQ":    1,
            "DSDECOD":  "DISEASE PROGRESSION" if pfs_events[i] else "CENSORED",
            "DSCAT":    "DISPOSITION EVENT",
            "DSSCAT":   "PROGRESSION FREE SURVIVAL",
            "DSSTDTC":  iso_date(int(pfs_times[i])),
            "DSSTDY":   int(pfs_times[i]) + 1,
            "DSTERM":   "DISEASE PROGRESSION" if pfs_events[i] else "ONGOING",
            "EPOCH":    "TREATMENT",
            "EVENTFL":  "Y" if pfs_events[i] else "N",
            "CNSR":     0 if pfs_events[i] else 1,
            "AVAL":     round(pfs_times[i] / 28, 2),  # months
            "PARAM":    "Progression-Free Survival",
            "PARAMCD":  "PFS",
        })
        # OS event
        rows.append({
            "STUDYID":  study_id,
            "DOMAIN":   "DS",
            "USUBJID":  subj["USUBJID"],
            "DSSEQ":    2,
            "DSDECOD":  "DEATH" if os_events[i] else "CENSORED",
            "DSCAT":    "DISPOSITION EVENT",
            "DSSCAT":   "OVERALL SURVIVAL",
            "DSSTDTC":  iso_date(int(os_times[i])),
            "DSSTDY":   int(os_times[i]) + 1,
            "DSTERM":   "DEATH" if os_events[i] else "ONGOING",
            "EPOCH":    "TREATMENT",
            "EVENTFL":  "Y" if os_events[i] else "N",
            "CNSR":     0 if os_events[i] else 1,
            "AVAL":     round(os_times[i] / 28, 2),
            "PARAM":    "Overall Survival",
            "PARAMCD":  "OS",
        })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# ADaM: ADSL — Subject Level Analysis Dataset
# ─────────────────────────────────────────────────────────────────────────────

def make_adsl(dm: pd.DataFrame, ds: pd.DataFrame, study_key: str) -> pd.DataFrame:
    study_id = STUDY_CONFIG[study_key]["studyid"]

    pfs = ds[ds["PARAMCD"] == "PFS"][["USUBJID","AVAL","CNSR","DSSTDTC"]].rename(
        columns={"AVAL":"PFSDUR","CNSR":"PFSCNSR","DSSTDTC":"PFSDTC"})
    os_ = ds[ds["PARAMCD"] == "OS"][["USUBJID","AVAL","CNSR","DSSTDTC"]].rename(
        columns={"AVAL":"OSDUR","CNSR":"OSCNSR","DSSTDTC":"OSDTC"})

    adsl = dm.merge(pfs, on="USUBJID", how="left")\
             .merge(os_, on="USUBJID", how="left")

    adsl["STUDYID"]   = study_id
    adsl["TRTSDT"]    = adsl["RFSTDTC"]
    adsl["TRT01P"]    = adsl["ARM"]
    adsl["TRT01A"]    = adsl["ACTARM"]
    adsl["TRT01PN"]   = np.where(adsl["ARMCD"] == "IRd", 1, 2)
    adsl["TRT01AN"]   = adsl["TRT01PN"]
    adsl["ITTFL"]     = "Y"
    adsl["SAFFL"]     = "Y"
    adsl["PPSFL"]     = RNG.choice(["Y","N"], size=len(adsl), p=[0.92,0.08])
    adsl["EOSSTT"]    = np.where(adsl["PFSCNSR"] == 0, "DISCONTINUED", "ONGOING")
    adsl["DCSREAS"]   = np.where(adsl["PFSCNSR"] == 0, "DISEASE PROGRESSION", "")

    # ISS risk grouping
    adsl["RISKGR"]    = pd.cut(adsl["ISSSTAGE"], bins=[0,1,2,3],
                               labels=["LOW","INTERMEDIATE","HIGH"])

    # Cytogenetics (simplified)
    adsl["CYTOGR"]    = RNG.choice(["HIGH RISK","STANDARD RISK"],
                                   size=len(adsl), p=[0.22, 0.78])
    adsl["DEL17P"]    = RNG.choice(["Y","N"], size=len(adsl), p=[0.10,0.90])
    adsl["T414"]      = RNG.choice(["Y","N"], size=len(adsl), p=[0.15,0.85])
    adsl["AMP1Q"]     = RNG.choice(["Y","N"], size=len(adsl), p=[0.25,0.75])

    # Polymorphisms (CC, CG, GG)
    poly = RNG.choice(["CC","CG","GG"], size=len(adsl), p=[0.35,0.45,0.20])
    adsl["POLYCC"]    = (poly == "CC").astype(int)
    adsl["POLYCG"]    = (poly == "CG").astype(int)
    adsl["POLYGG"]    = (poly == "GG").astype(int)

    # Durie-Salmon stage
    adsl["DSSTAGE"]   = RNG.choice(["IA","IIA","IIB","IIIA","IIIB"],
                                   size=len(adsl), p=[0.20,0.25,0.10,0.35,0.10])

    # Baseline labs
    adsl["BASE_CREAT"] = np.round(RNG.normal(88, 25, len(adsl)), 1)   # umol/L
    adsl["BASE_CREACL"] = np.round(RNG.normal(65, 20, len(adsl)), 1)  # mL/min
    adsl["BASE_ALB"]   = np.round(RNG.normal(38, 5,  len(adsl)), 1)   # g/L
    adsl["BASE_B2MG"]  = RNG.choice(["<3.5","3.5-5.5",">5.5"],
                                    size=len(adsl), p=[0.46,0.37,0.17])
    adsl["BASE_FLC"]   = np.round(np.abs(RNG.exponential(120, len(adsl))), 1)  # mg/L
    adsl["BASE_BMPC"]  = np.round(RNG.uniform(5, 90, len(adsl)), 1)   # %
    adsl["BASE_PC"]    = np.round(RNG.uniform(5, 90, len(adsl)), 1)   # %
    adsl["EXTRAMED"]   = RNG.choice(["Y","N"], size=len(adsl), p=[0.08,0.92])
    adsl["BONELESI"]   = RNG.choice(["Y","N"], size=len(adsl), p=[0.55,0.45])
    adsl["PLASMA"]     = RNG.choice(["Y","N"], size=len(adsl), p=[0.12,0.88])
    adsl["MEAS_DISFL"] = RNG.choice(["SERUM M-PROTEIN","URINE M-PROTEIN",
                                      "BOTH SPEP AND UPEP","SERUM FLC"],
                                    size=len(adsl), p=[0.55,0.10,0.15,0.20])

    return adsl


# ─────────────────────────────────────────────────────────────────────────────
# ADaM: ADTTE — Time-to-Event Analysis Dataset
# ─────────────────────────────────────────────────────────────────────────────

def make_adtte(ds: pd.DataFrame, adsl: pd.DataFrame, study_key: str) -> pd.DataFrame:
    study_id = STUDY_CONFIG[study_key]["studyid"]

    adtte = ds[["STUDYID","USUBJID","PARAMCD","PARAM","AVAL","CNSR","DSSTDTC"]].copy()
    adtte = adtte.merge(adsl[["USUBJID","ARMCD","ARM","AGE","SEX","IGTYPE",
                               "ISSSTAGE","RISKGR","CYTOGR","TRT01PN"]], on="USUBJID")
    adtte["STUDYID"]  = study_id
    adtte["AVALU"]    = "MONTHS"
    adtte["STARTDT"]  = adsl.set_index("USUBJID").loc[adtte["USUBJID"], "RFSTDTC"].values
    adtte["EVNTDESC"] = np.where(adtte["PARAMCD"] == "PFS",
                                 np.where(adtte["CNSR"] == 0, "DISEASE PROGRESSION", "CENSORED"),
                                 np.where(adtte["CNSR"] == 0, "DEATH", "CENSORED"))
    adtte["SRCDOM"]   = "DS"
    adtte["ITTFL"]    = "Y"
    adtte["SAFFL"]    = "Y"

    # Stratification variables (used in stratified Cox in paper)
    adtte["STRATFL1"] = adtte["ISSSTAGE"].apply(lambda x: "I" if x==1 else ("II" if x==2 else "III"))
    adtte["STRATFL2"] = adtte["CYTOGR"]
    adtte["STRATFL3"] = adtte["IGTYPE"]

    return adtte


# ─────────────────────────────────────────────────────────────────────────────
# ADaM: ADLB — Longitudinal Lab Analysis Dataset
# ─────────────────────────────────────────────────────────────────────────────

def make_adlb(lb: pd.DataFrame, adsl: pd.DataFrame, study_key: str) -> pd.DataFrame:
    study_id = STUDY_CONFIG[study_key]["studyid"]

    adlb = lb.copy()
    adlb = adlb.merge(adsl[["USUBJID","ARMCD","TRT01PN","IGTYPE","ISSSTAGE"]],
                      on="USUBJID", how="left")
    adlb["STUDYID"] = study_id

    # Baseline value
    baseline = lb[lb["EPOCH"] == "BASELINE"][["USUBJID","LBTESTCD","LBSTRESN"]]\
               .rename(columns={"LBSTRESN": "BASE"})
    adlb = adlb.merge(baseline, on=["USUBJID","LBTESTCD"], how="left")

    # Change from baseline
    adlb["CHG"]    = adlb["LBSTRESN"] - adlb["BASE"]
    adlb["PCHG"]   = np.where(adlb["BASE"] != 0,
                               (adlb["CHG"] / adlb["BASE"]) * 100, np.nan)
    adlb["AVAL"]   = adlb["LBSTRESN"]
    adlb["DTYPE"]  = ""
    adlb["ANL01FL"] = "Y"
    adlb["PARAM"]   = adlb["LBTEST"]
    adlb["PARAMCD"] = adlb["LBTESTCD"]
    adlb["AVALU"]   = adlb["LBSTRESU"]
    adlb["ABLFL"]   = np.where(adlb["EPOCH"] == "BASELINE", "Y", "")

    return adlb


# ─────────────────────────────────────────────────────────────────────────────
# MAIN — Generate All Datasets
# ─────────────────────────────────────────────────────────────────────────────

def generate_all():
    import os
    out_dir = "/mnt/user-data/outputs/tourmaline_synthetic"
    os.makedirs(out_dir, exist_ok=True)

    summary = []

    for study_key in ["MM2", "MM1"]:
        print(f"\n{'='*60}")
        print(f"  Generating {study_key} ({STUDY_CONFIG[study_key]['studyid']})")
        print(f"  N = {STUDY_CONFIG[study_key]['n']}, Population = {STUDY_CONFIG[study_key]['population']}")
        print(f"{'='*60}")

        study_dir = os.path.join(out_dir, study_key)
        os.makedirs(study_dir, exist_ok=True)

        # SDTM
        print("  [1/7] DM  — Demographics")
        dm = make_dm(study_key)
        dm.to_csv(f"{study_dir}/sdtm_dm.csv", index=False)

        print("  [2/7] EX  — Exposure/Dosing  (may take a moment...)")
        # Use subset for EX to avoid huge file
        dm_ex_sample = dm.sample(min(200, len(dm)), random_state=42)
        ex = make_ex(dm_ex_sample, study_key)
        ex.to_csv(f"{study_dir}/sdtm_ex.csv", index=False)

        print("  [3/7] LB  — Laboratory Results")
        dm_lb_sample = dm.sample(min(300, len(dm)), random_state=42)
        lb = make_lb(dm_lb_sample, study_key)
        lb.to_csv(f"{study_dir}/sdtm_lb.csv", index=False)

        print("  [4/7] AE  — Adverse Events")
        ae = make_ae(dm, study_key)
        ae.to_csv(f"{study_dir}/sdtm_ae.csv", index=False)

        print("  [5/7] DS  — Disposition / Survival")
        ds = make_ds(dm, study_key)
        ds.to_csv(f"{study_dir}/sdtm_ds.csv", index=False)

        # ADaM
        print("  [6/7] ADSL — Subject Level")
        adsl = make_adsl(dm, ds, study_key)
        adsl.to_csv(f"{study_dir}/adam_adsl.csv", index=False)

        print("  [7/7] ADTTE / ADLB — TTE & Longitudinal")
        adtte = make_adtte(ds, adsl, study_key)
        adtte.to_csv(f"{study_dir}/adam_adtte.csv", index=False)

        adlb = make_adlb(lb, adsl, study_key)
        adlb.to_csv(f"{study_dir}/adam_adlb.csv", index=False)

        # Collect summary
        summary.append({
            "Study": study_key,
            "Studyid": STUDY_CONFIG[study_key]["studyid"],
            "Population": STUDY_CONFIG[study_key]["population"],
            "N_subjects": len(dm),
            "DM_rows": len(dm),
            "EX_rows": len(ex),
            "LB_rows": len(lb),
            "AE_rows": len(ae),
            "DS_rows": len(ds),
            "ADSL_rows": len(adsl),
            "ADTTE_rows": len(adtte),
            "ADLB_rows": len(adlb),
        })

        print(f"\n  Summary for {study_key}:")
        print(f"    DM: {len(dm):,} subjects")
        print(f"    EX: {len(ex):,} exposure records")
        print(f"    LB: {len(lb):,} lab records")
        print(f"    AE: {len(ae):,} adverse events")
        print(f"    DS: {len(ds):,} disposition records")
        print(f"    ADSL: {len(adsl):,} | ADTTE: {len(adtte):,} | ADLB: {len(adlb):,}")

    # Write combined summary
    pd.DataFrame(summary).to_csv(f"{out_dir}/dataset_summary.csv", index=False)
    print(f"\n✓ All datasets written to: {out_dir}")
    return out_dir


if __name__ == "__main__":
    generate_all()
