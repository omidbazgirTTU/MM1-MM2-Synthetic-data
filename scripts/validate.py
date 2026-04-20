"""
TOURMALINE Synthetic Data Validator
=====================================
Compares generated dataset statistics against published TOURMALINE trial values.
Prints a pass/warn/fail report for each check.

Usage:
    python scripts/validate.py

Reference targets:
  - Moreau et al. NEJM 2016 (TOURMALINE-MM2): Table 1, Fig 1, Table S4
  - Moreau et al. Lancet Oncol 2016 (TOURMALINE-MM1): Table 1, Fig 1
  - Gupta et al. Clin Pharmacokinet 2017 (Ixazomib popPK)
"""

import os
import sys
import numpy as np
import pandas as pd
try:
    from lifelines import KaplanMeierFitter
    HAS_LIFELINES = True
except ImportError:
    HAS_LIFELINES = False

def km_median(aval, cnsr):
    """Return KM median survival time. Falls back to raw median if lifelines unavailable."""
    if HAS_LIFELINES:
        kmf = KaplanMeierFitter()
        kmf.fit(aval, event_observed=cnsr.eq(0))
        return kmf.median_survival_time_
    return aval.median()  # biased fallback

BASE = os.path.join(os.path.dirname(__file__), "..")

PASS  = "\033[92m  PASS\033[0m"
WARN  = "\033[93m  WARN\033[0m"
FAIL  = "\033[91m  FAIL\033[0m"

def check(label, value, target, tol_pct=15.0, fmt=".1f"):
    """Print pass/warn/fail for a numeric check."""
    if target == 0:
        err_pct = abs(value - target)
    else:
        err_pct = abs(value - target) / abs(target) * 100
    icon = PASS if err_pct <= tol_pct else (WARN if err_pct <= tol_pct * 2 else FAIL)
    print(f"{icon}  {label:<52} got={value:{fmt}}  target={target:{fmt}}  "
          f"err={err_pct:.0f}%")

def check_range(label, value, lo, hi, fmt=".1f"):
    icon = PASS if lo <= value <= hi else FAIL
    print(f"{icon}  {label:<52} got={value:{fmt}}  target=[{lo:{fmt}}, {hi:{fmt}}]")

def section(title):
    print(f"\n{'─'*70}")
    print(f"  {title}")
    print(f"{'─'*70}")

results = {"pass": 0, "warn": 0, "fail": 0}


def validate_study(study_key):
    sdir = os.path.join(BASE, study_key)
    is_ndmm = (study_key == "MM2")

    # Published targets
    if is_ndmm:
        tgt = {
            "median_age": 73, "pct_female": 49.9,
            "pfs_med_ird": 20.6, "pfs_med_rd": 14.7,
            "os_med_ird": 53.6, "os_med_rd": 51.6,
            "pfs_event_pct": 62.0,
            "iss1_pct": 46.1, "iss3_pct": 16.4,
            "igg_pct": 57.3,
            "hr_cyto_pct": 22.0,
            "median_weight": 72.0, "bsa_range": (1.55, 2.10),
            "crcl_range": (30.0, 120.0), "crcl_median": 62.0,
            # Palumbo 2015 NEJM (R-ISS development cohort, N=4445): I=28%, II=62%, III=10%
            "riss_i_pct": 28.0, "riss_iii_pct": 10.0,
            # AE incidences (any grade, published approximate %)
            "ae_diarrhea_ird": 42.0, "ae_neut_ird": 26.0,
            "ae_thrombo_ird": 19.0, "ae_pn_ird": 19.0,
        }
    else:
        tgt = {
            "median_age": 66, "pct_female": 43.3,
            "pfs_med_ird": 20.6, "pfs_med_rd": 14.7,
            "os_med_ird": 39.0, "os_med_rd": 33.0,
            "pfs_event_pct": 62.0,
            "iss1_pct": 63.6, "iss3_pct": 12.0,
            "igg_pct": 54.0,
            "hr_cyto_pct": 22.0,
            "median_weight": 75.0, "bsa_range": (1.55, 2.15),
            "crcl_range": (30.0, 130.0), "crcl_median": 68.0,
            "riss_i_pct": 28.0, "riss_iii_pct": 10.0,
            "ae_diarrhea_ird": 45.0, "ae_neut_ird": 26.0,
            "ae_thrombo_ird": 20.0, "ae_pn_ird": 21.0,
        }

    print(f"\n{'='*70}")
    print(f"  VALIDATING {study_key}  ({('NDMM / MM2' if is_ndmm else 'RRMM / MM1')})")
    print(f"{'='*70}")

    # ── Demographics ──────────────────────────────────────────────────────────
    section("Demographics (sdtm_dm.csv)")
    dm = pd.read_csv(f"{sdir}/sdtm_dm.csv")

    check("Median age (years)",        dm["AGE"].median(),              tgt["median_age"],    tol_pct=5)
    check("% Female",                  (dm["SEX"]=="F").mean()*100,     tgt["pct_female"],    tol_pct=10)
    check("% ISS Stage I",             (dm["ISSSTAGE"]==1).mean()*100,  tgt["iss1_pct"],      tol_pct=15)
    check("% ISS Stage III",           (dm["ISSSTAGE"]==3).mean()*100,  tgt["iss3_pct"],      tol_pct=20)
    check("% IgG subtype",             (dm["IGTYPE"]=="IgG").mean()*100,tgt["igg_pct"],       tol_pct=10)
    check("Median weight (kg)",        dm["WEIGHT"].median(),           tgt["median_weight"], tol_pct=15)
    check_range("BSA range (m²)",      dm["BSA"].median(),              *tgt["bsa_range"])

    # ── Survival ──────────────────────────────────────────────────────────────
    section("Survival (adam_adtte.csv)")
    adtte = pd.read_csv(f"{sdir}/adam_adtte.csv")
    pfs   = adtte[adtte["PARAMCD"] == "PFS"]
    os_   = adtte[adtte["PARAMCD"] == "OS"]

    # Use KM median (accounts for censoring) — raw median is biased low
    pfs_ird = pfs[pfs["ARMCD"]=="IRd"]
    pfs_rd  = pfs[pfs["ARMCD"]=="Rd"]
    os_ird  = os_[os_["ARMCD"]=="IRd"]
    os_rd   = os_[os_["ARMCD"]=="Rd"]
    med_pfs_ird = km_median(pfs_ird["AVAL"], pfs_ird["CNSR"])
    med_pfs_rd  = km_median(pfs_rd["AVAL"],  pfs_rd["CNSR"])
    med_os_ird  = km_median(os_ird["AVAL"],  os_ird["CNSR"])
    med_os_rd   = km_median(os_rd["AVAL"],   os_rd["CNSR"])
    pfs_ev_pct  = pfs["CNSR"].eq(0).mean() * 100

    check("PFS median IRd (months)",   med_pfs_ird,  tgt["pfs_med_ird"],  tol_pct=15)
    check("PFS median Rd  (months)",   med_pfs_rd,   tgt["pfs_med_rd"],   tol_pct=15)
    check("OS  median IRd (months)",   med_os_ird,   tgt["os_med_ird"],   tol_pct=20)
    check("OS  median Rd  (months)",   med_os_rd,    tgt["os_med_rd"],    tol_pct=20)
    check("PFS event rate (%)",        pfs_ev_pct,   tgt["pfs_event_pct"],tol_pct=15)

    # Crude HR estimate (log of ratio of medians — approximate for Weibull)
    hr_approx = med_pfs_rd / med_pfs_ird
    check("PFS HR approx (Rd/IRd medians)",  hr_approx,  0.742, tol_pct=20, fmt=".3f")

    # ── ADSL ─────────────────────────────────────────────────────────────────
    section("ADSL derived variables (adam_adsl.csv)")
    adsl = pd.read_csv(f"{sdir}/adam_adsl.csv")

    check_range("Median CrCL (mL/min)", adsl["BASE_CREACL"].median(), *tgt["crcl_range"])
    check("Median CrCL (mL/min)",       adsl["BASE_CREACL"].median(), tgt["crcl_median"], tol_pct=25)
    check("% High-risk cytogenetics",   (adsl["CYTOGR"]=="HIGH RISK").mean()*100, tgt["hr_cyto_pct"], tol_pct=20)

    if "RISS" in adsl.columns:
        riss = adsl["RISS"].value_counts(normalize=True) * 100
        check("% R-ISS I",   riss.get("I",  0), tgt["riss_i_pct"],   tol_pct=30)
        check("% R-ISS III", riss.get("III", 0), tgt["riss_iii_pct"], tol_pct=30)
    else:
        print(f"{WARN}  R-ISS column not found in ADSL")

    # Extended cytogenetics present?
    for col, exp_pct in [("T1416", 3.0), ("T1420", 1.0), ("GAIN1Q21", 35.0), ("DEL1P32", 10.0)]:
        if col in adsl.columns:
            check(f"% {col}",  (adsl[col]=="Y").mean()*100, exp_pct, tol_pct=50)
        else:
            print(f"{WARN}  {col} column missing in ADSL")

    # CYP3A4 DDI flags
    if "CYP3A4_INHIBFL" in adsl.columns:
        pct_inh = (adsl["CYP3A4_INHIBFL"]=="Y").mean() * 100
        check_range("% CYP3A4 inhibitor use", pct_inh, 5.0, 22.0)
    else:
        print(f"{WARN}  CYP3A4_INHIBFL missing in ADSL")

    # ── Concomitant Meds ──────────────────────────────────────────────────────
    section("Concomitant Medications (sdtm_cm.csv)")
    cm_path = f"{sdir}/sdtm_cm.csv"
    if os.path.exists(cm_path):
        cm  = pd.read_csv(cm_path)
        n   = len(dm)
        inh_pct = len(cm[cm["CMCAT"]=="CYP3A4 INHIBITOR"]["USUBJID"].unique()) / n * 100
        acg_pct = len(cm[cm["CMCAT"]=="ANTICOAGULANT / ANTIPLATELET"]["USUBJID"].unique()) / n * 100
        sc_pct  = len(cm[cm["CMCAT"]=="SUPPORTIVE CARE"]["USUBJID"].unique()) / n * 100
        check_range("% subjects with CYP3A4 inhibitor", inh_pct,  5.0, 22.0)
        check_range("% subjects with anticoagulant",    acg_pct, 50.0, 95.0)
        check_range("% subjects with supportive care",  sc_pct,  50.0, 99.0)
    else:
        print(f"{FAIL}  sdtm_cm.csv not found")

    # ── Adverse Events ────────────────────────────────────────────────────────
    section("Adverse Events (sdtm_ae.csv)")
    ae  = pd.read_csv(f"{sdir}/sdtm_ae.csv")
    n   = len(dm)
    ird_uids = set(dm[dm["ARMCD"]=="IRd"]["USUBJID"])

    def ae_pct(ae_name, arm_uids):
        return len(ae[(ae["AEDECOD"]==ae_name) & (ae["USUBJID"].isin(arm_uids))]["USUBJID"].unique()) \
               / len(arm_uids) * 100

    check("% Diarrhea    (IRd)",         ae_pct("Diarrhea", ird_uids),              tgt["ae_diarrhea_ird"], tol_pct=30)
    check("% Neutropenia (IRd)",         ae_pct("Neutropenia", ird_uids),            tgt["ae_neut_ird"],     tol_pct=30)
    check("% Thrombocytopenia (IRd)",    ae_pct("Thrombocytopenia", ird_uids),       tgt["ae_thrombo_ird"],  tol_pct=30)
    check("% Periph. Neuropathy (IRd)",  ae_pct("Peripheral Neuropathies", ird_uids),tgt["ae_pn_ird"],       tol_pct=30)

    # ── Lab Biomarkers ────────────────────────────────────────────────────────
    section("Longitudinal Biomarkers (sdtm_lb.csv)")
    lb = pd.read_csv(f"{sdir}/sdtm_lb.csv")

    # Baseline HGB (expected ~10.5 g/dL = 105 g/L in MM)
    bl_hgb = lb[(lb["LBTESTCD"]=="HGB") & (lb["EPOCH"]=="BASELINE")]["LBSTRESN"]
    if len(bl_hgb):
        check("Baseline median HGB (g/L)", bl_hgb.median(), 105.0, tol_pct=20)
        check_range("Baseline HGB range",  bl_hgb.median(), 80.0, 130.0)
    else:
        print(f"{WARN}  No baseline HGB records found")

    # M-protein response: best % change from baseline (IRd vs Rd)
    mp = lb[lb["LBTESTCD"]=="SPEP_MPROT"].copy()
    if len(mp):
        bl_mp = mp[mp["EPOCH"]=="BASELINE"][["USUBJID","LBSTRESN"]]\
                .rename(columns={"LBSTRESN":"BASE"})
        mp = mp.merge(bl_mp, on="USUBJID", how="left")
        mp = mp.merge(dm[["USUBJID","ARMCD"]], on="USUBJID", how="left")
        mp["PCHG"] = (mp["LBSTRESN"] - mp["BASE"]) / mp["BASE"].replace(0, np.nan) * 100
        best_pchg = mp.groupby(["USUBJID","ARMCD"])["PCHG"].min().reset_index()
        best_ird  = best_pchg[best_pchg["ARMCD"]=="IRd"]["PCHG"].median()
        best_rd   = best_pchg[best_pchg["ARMCD"]=="Rd"]["PCHG"].median()
        # Bimodal: 18-28% non-resp @ ~0%, rest deep responders → overall median ~−88/−78%
        # Consistent with TOURMALINE VGPR+ 63%/55% (MM2) and 48%/39% (MM1)
        check("Median best M-protein %chg IRd", best_ird, -90.0, tol_pct=15)
        check("Median best M-protein %chg Rd",  best_rd,  -78.0, tol_pct=15)
    else:
        print(f"{WARN}  No SPEP_MPROT records found")

    # New biomarkers present
    for test in ["B2MG", "LDH", "UPEP_MPROT"]:
        n_records = (lb["LBTESTCD"] == test).sum()
        if n_records > 0:
            print(f"{PASS}  {test:<20} present  ({n_records:,} records)")
        else:
            print(f"{FAIL}  {test:<20} MISSING")

    # Missing rate
    miss_rate = lb.groupby("LBTESTCD").size()
    expected_complete = len(dm) * 12  # rough: 12 cycles × all subjects
    overall_miss = 1 - lb["LBSTRESN"].notna().mean()
    check_range("Overall LB missing rate", overall_miss * 100, 10.0, 20.0)

    # ── EX ────────────────────────────────────────────────────────────────────
    section("Exposure / Dosing (sdtm_ex.csv)")
    ex = pd.read_csv(f"{sdir}/sdtm_ex.csv")
    all_subj = dm["USUBJID"].nunique()
    ex_subj  = ex["USUBJID"].nunique()
    check("EX covers all DM subjects (%)", ex_subj / all_subj * 100, 100.0, tol_pct=1)

    if "EXDOSMOD" in ex.columns:
        dose_mod_pct = (ex["EXDOSMOD"].ne("") & ex["EXDOSMOD"].notna()).mean() * 100
        check_range("% dose records with reduction", dose_mod_pct, 10.0, 25.0)

    if "EXADH" in ex.columns:
        print(f"{PASS}  Adherence (EXADH) column present")
    else:
        print(f"{WARN}  EXADH (adherence) column missing")

    if "EXDOSMOD_REASON" in ex.columns:
        print(f"{PASS}  Dose-mod reason (EXDOSMOD_REASON) column present")
    else:
        print(f"{WARN}  EXDOSMOD_REASON column missing")

    print()


def main():
    for sk in ["MM2", "MM1"]:
        path = os.path.join(BASE, sk)
        if not os.path.exists(path):
            print(f"\n[SKIP] {sk} directory not found: {path}")
            continue
        validate_study(sk)

    print(f"\n{'='*70}")
    print("  Validation complete. Review WARN/FAIL items above.")
    print("  Tolerance: PASS ≤15%, WARN ≤30%, FAIL >30% from target.")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
