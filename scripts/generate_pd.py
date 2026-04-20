"""
TOURMALINE Synthetic PD Data Generator
=======================================
Adds the missing pharmacodynamic (PD) domains to the existing synthetic dataset.

New outputs
-----------
  sdtm_lb.csv      — BMPC (bone marrow plasma cell %) records appended
  adam_adlb.csv    — BMPC rows appended
  sdtm_rs.csv      — NEW: IMWG response assessments per cycle per subject
  adam_adrs.csv    — NEW: Best overall response, time-to-response, duration-of-response

IMWG 2016 Response Criteria  (Kumar et al. Lancet Oncol 2016, PMID 27511158)
-----------------------------------------------------------------------------
  sCR  : CR + normal FLC ratio + BM <5% PC (no clonal cells by IHC)
  CR   : Negative serum+urine immunofixation (IF) + <5% BM plasma cells
  VGPR : ≥90% reduction in serum M-protein (OR IF positive but SPEP undetectable)
  PR   : ≥50% reduction in serum M-protein
  MR   : ≥25% but <49% reduction in serum M-protein
  SD   : Does not meet criteria for PR or PD
  PD   : ≥25% increase from lowest confirmed value (nadir) with abs increase ≥0.5 g/dL

Immunofixation (IF) negativity is modelled as:
  IF_NEG = (serum M-protein < 0.5 g/L) AND (urine M-protein < 20 mg/day)

FLC ratio normalisation (required for sCR) is modelled as:
  FLC_NORMAL = (KAPPA_LAMBDA ratio between 0.26 and 1.65) after treatment response

Response assessments are performed at every cycle (per TOURMALINE trial protocol),
with the caveat that BMPC (bone marrow biopsy) is assessed at:
  Baseline · Cycle 3 · Cycle 6 · Cycle 12 · End of Treatment

Seed: 42 (reproducible)
"""

import os
import sys
import numpy as np
import pandas as pd

BASE = os.path.join(os.path.dirname(__file__), "..")
RNG  = np.random.default_rng(seed=42)

# ─── Response ordinal ranking (higher = better) ───────────────────────────────
RESP_RANK = {"NE": 0, "PD": 1, "SD": 2, "MR": 3, "PR": 4, "VGPR": 5, "CR": 6, "sCR": 7}
RESP_LABEL = {
    "sCR":  "Stringent Complete Response",
    "CR":   "Complete Response",
    "VGPR": "Very Good Partial Response",
    "PR":   "Partial Response",
    "MR":   "Minimal Response",
    "SD":   "Stable Disease",
    "PD":   "Progressive Disease",
    "NE":   "Not Evaluable",
}

# IMWG threshold for PD: ≥25% increase from nadir, absolute increase ≥ 5 g/L
PD_REL_THRESHOLD  = 0.25   # 25% increase from nadir
PD_ABS_THRESHOLD  = 5.0    # g/L  (≈ 0.5 g/dL)

# BMPC target distributions (% plasma cells in bone marrow)
# Baseline in MM: typically 20–80% (median ~40%)
# Normal: <5%
BMPC_ASSESSMENT_VISITS = {1, 3, 6, 12, 24}   # cycles with bone marrow biopsy
BMPC_MISSING_RATE = 0.22                       # invasive — more often missed

# Published response rates (TOURMALINE) — used to document targets, not hard-coded
RESP_TARGETS = {
    "MM2": {   # NDMM — Moreau 2019 NEJM
        "IRd": {"ORR": 82.0, "VGPR_plus": 63.0, "CR_plus": 28.0, "sCR": 15.0},
        "Rd":  {"ORR": 75.0, "VGPR_plus": 55.0, "CR_plus": 14.0, "sCR":  8.0},
    },
    "MM1": {   # RRMM — Moreau 2016 Lancet Oncol
        "IRd": {"ORR": 78.3, "VGPR_plus": 48.1, "CR_plus": 11.7, "sCR":  5.0},
        "Rd":  {"ORR": 71.5, "VGPR_plus": 38.8, "CR_plus":  6.6, "sCR":  3.0},
    },
}


# ─── Helpers ──────────────────────────────────────────────────────────────────

def iso_date(day_offset: int, ref: str = "2016-01-01") -> str:
    d = pd.Timestamp(ref) + pd.Timedelta(days=int(day_offset))
    return d.strftime("%Y-%m-%d")


def _imwg_response(pchg_from_base: float, pchg_from_nadir: float,
                   abs_from_nadir: float, if_neg: bool,
                   bmpc_lt5: bool, flc_normal: bool) -> str:
    """
    Apply IMWG 2016 criteria to assign a response category.

    pchg_from_base  : % change in M-protein from baseline (negative = reduction)
    pchg_from_nadir : % change from nadir to current (positive = increase from nadir)
    abs_from_nadir  : absolute increase in M-protein from nadir (g/L)
    if_neg          : True if serum+urine immunofixation negative
    bmpc_lt5        : True if bone marrow plasma cells < 5%
    flc_normal      : True if serum FLC ratio normal (0.26-1.65)
    """
    # PD — ≥25% increase from nadir AND absolute increase ≥ 5 g/L
    if pchg_from_nadir >= PD_REL_THRESHOLD * 100 and abs_from_nadir >= PD_ABS_THRESHOLD:
        return "PD"

    if pd.isna(pchg_from_base):
        return "NE"

    if pchg_from_base > -25.0:
        return "SD"

    if pchg_from_base > -50.0:
        return "MR"

    if pchg_from_base > -90.0:
        return "PR"

    # Deep response (≥90% reduction)
    if if_neg and bmpc_lt5 and flc_normal:
        return "sCR"
    if if_neg and bmpc_lt5:
        return "CR"
    return "VGPR"   # ≥90% but not meeting CR criteria


# ─── BMPC simulation ──────────────────────────────────────────────────────────

def _sim_bmpc(baseline_mp: float, mp_vals: np.ndarray,
              baseline_mp_pop: float, arm: str) -> np.ndarray:
    """
    Simulate bone marrow plasma cell % at each cycle.

    Correlation with M-protein:  r ~ 0.70 (published; Hillengass 2019)
    Baseline BMPC:  log-normal, correlated with serum M-protein
    """
    n = len(mp_vals)
    if n == 0:
        return np.array([])

    # Baseline BMPC correlated with baseline M-protein level
    # Typical MM baseline BMPC: 20–80% (lognormal; median ~40%)
    mp_frac = np.clip(baseline_mp / (baseline_mp_pop + 1e-6), 0.3, 3.0)
    bmpc_base = np.clip(
        RNG.lognormal(mean=np.log(35.0 * mp_frac), sigma=0.45),
        5.0, 95.0
    )

    # Response: BMPC mirrors M-protein reduction with some lag
    bmpc_vals = np.zeros(n)
    mp_pchg   = (mp_vals - baseline_mp) / (baseline_mp + 1e-6) * 100

    for i, pchg in enumerate(mp_pchg):
        lag_factor = 0.80 if arm == "IRd" else 0.65   # IRd clears BM faster
        raw_reduction = np.clip(pchg * lag_factor, -100, 50)   # % of baseline
        bmpc_i = bmpc_base * (1 + raw_reduction / 100.0)
        # Add noise + floor at 0.1
        bmpc_vals[i] = np.clip(
            bmpc_i + RNG.normal(0, bmpc_base * 0.08),
            0.1, 95.0
        )

    return np.round(bmpc_vals, 1)


# ─── IMWG response derivation per subject ────────────────────────────────────

def _derive_response_series(subj_mp: pd.DataFrame,
                             subj_upep: pd.DataFrame,
                             bmpc_map: dict,
                             flc_vals: dict,
                             arm: str) -> pd.DataFrame:
    """
    For one subject: derive per-cycle IMWG response from M-protein, UPEP, BMPC, FLC.

    Returns DataFrame with columns:
        VISITNUM, RSC_IMWG, SPEP_VAL, PCHG_BASE, PCHG_NADIR, IF_NEG, BMPC_LT5, FLC_NORM
    """
    if subj_mp.empty:
        return pd.DataFrame()

    subj_mp = subj_mp.sort_values("VISITNUM").dropna(subset=["LBSTRESN"])
    if subj_mp.empty:
        return pd.DataFrame()

    bl_row = subj_mp[subj_mp["EPOCH"] == "BASELINE"]
    if bl_row.empty:
        return pd.DataFrame()

    bl_mp = bl_row["LBSTRESN"].iloc[0]

    # Merge UPEP baseline
    if not subj_upep.empty:
        bl_upep_row = subj_upep[subj_upep["EPOCH"] == "BASELINE"]
        bl_upep = bl_upep_row["LBSTRESN"].iloc[0] if not bl_upep_row.empty else 100.0
    else:
        bl_upep = 100.0

    rows = []
    nadir_mp  = bl_mp      # track nadir for PD detection

    for _, row in subj_mp.iterrows():
        vis   = int(row["VISITNUM"])
        mp_v  = row["LBSTRESN"]
        epoch = row["EPOCH"]

        # Update nadir
        if mp_v < nadir_mp:
            nadir_mp = mp_v

        # % change from baseline
        pchg_base  = (mp_v - bl_mp) / (bl_mp + 1e-6) * 100

        # % change from nadir
        pchg_nadir = (mp_v - nadir_mp) / (nadir_mp + 1e-6) * 100
        abs_nadir  = mp_v - nadir_mp     # g/L absolute increase from nadir

        # Immunofixation negative: very low absolute M-protein level
        # IF can detect down to ~0.02 g/L; we model IF_NEG when M-protein < 0.5 g/L
        # and UPEP < 20 mg/day (if available)
        upep_row = subj_upep[subj_upep["VISITNUM"] == vis]
        upep_v   = upep_row["LBSTRESN"].iloc[0] if (not upep_row.empty and
                   upep_row["LBSTRESN"].notna().any()) else np.nan
        if_neg   = (mp_v < 0.5) and (pd.isna(upep_v) or upep_v < 20.0)

        # BMPC < 5%
        bmpc_val  = bmpc_map.get(vis, np.nan)
        bmpc_lt5  = (not pd.isna(bmpc_val)) and (bmpc_val < 5.0)

        # FLC ratio normal (0.26–1.65) after treatment
        flc_v     = flc_vals.get(vis, np.nan)
        flc_norm  = (not pd.isna(flc_v)) and (0.26 <= flc_v <= 1.65)

        resp = _imwg_response(
            pchg_base, pchg_nadir, abs_nadir,
            if_neg, bmpc_lt5, flc_norm
        )

        rows.append({
            "VISITNUM":   vis,
            "EPOCH":      epoch,
            "RSC_IMWG":   resp,
            "SPEP_VAL":   round(mp_v, 3),
            "PCHG_BASE":  round(pchg_base, 1),
            "PCHG_NADIR": round(pchg_nadir, 1),
            "ABS_NADIR":  round(abs_nadir, 3),
            "IF_NEG":     if_neg,
            "BMPC_VAL":   bmpc_val,
            "BMPC_LT5":   bmpc_lt5,
            "FLC_VAL":    flc_v,
            "FLC_NORM":   flc_norm,
        })

    return pd.DataFrame(rows)


# ─── SDTM RS domain ───────────────────────────────────────────────────────────

def make_rs(dm: pd.DataFrame, lb: pd.DataFrame,
            bmpc_df: pd.DataFrame, study_key: str) -> pd.DataFrame:
    """
    Generate SDTM RS (Disease Response) domain.
    One row per subject per assessment visit.
    """
    sid  = study_key
    rows = []
    seq  = 1

    mp_all   = lb[lb["LBTESTCD"] == "SPEP_MPROT"].copy()
    upep_all = lb[lb["LBTESTCD"] == "UPEP_MPROT"].copy()
    flc_all  = lb[lb["LBTESTCD"] == "KAPPA_LAMBDA"].copy()

    # Build BMPC lookup per subject per visit
    bmpc_pivot = {}
    if not bmpc_df.empty:
        for uid, grp in bmpc_df.groupby("USUBJID"):
            bmpc_pivot[uid] = grp.dropna(subset=["LBSTRESN"])\
                                 .set_index("VISITNUM")["LBSTRESN"].to_dict()

    for _, subj in dm.iterrows():
        uid  = subj["USUBJID"]
        arm  = subj["ARMCD"]

        subj_mp   = mp_all[mp_all["USUBJID"] == uid]
        subj_upep = upep_all[upep_all["USUBJID"] == uid]
        subj_flc  = flc_all[flc_all["USUBJID"] == uid]

        bmpc_map = bmpc_pivot.get(uid, {})
        flc_map  = {}
        for _, fr in subj_flc.dropna(subset=["LBSTRESN"]).iterrows():
            flc_map[int(fr["VISITNUM"])] = fr["LBSTRESN"]

        resp_df = _derive_response_series(subj_mp, subj_upep, bmpc_map, flc_map, arm)
        if resp_df.empty:
            continue

        for _, r in resp_df.iterrows():
            rows.append({
                "STUDYID":   f"TOURMALINE-{study_key}",
                "DOMAIN":    "RS",
                "USUBJID":   uid,
                "RSSEQ":     seq,
                "RSTESTCD":  "OVRLRESP",
                "RSTEST":    "Overall Response",
                "RSCAT":     "EFFICACY",
                "RSORRES":   r["RSC_IMWG"],
                "RSSTRESC":  r["RSC_IMWG"],
                "RSEVAL":    "INVESTIGATOR",
                "RSEVALID":  "INVESTIGATOR",
                "VISITNUM":  r["VISITNUM"],
                "VISIT":     "BASELINE" if r["VISITNUM"] == 1 else f"CYCLE {int(r['VISITNUM'])}",
                "EPOCH":     r["EPOCH"],
                "RSDTC":     iso_date((r["VISITNUM"] - 1) * 28),
                "RSDY":      (r["VISITNUM"] - 1) * 28 + 1,
                # Analysis variables (denormalised for ML convenience)
                "SPEP_VAL":  r["SPEP_VAL"],
                "PCHG_BASE": r["PCHG_BASE"],
                "PCHG_NADIR":r["PCHG_NADIR"],
                "IF_NEG":    "Y" if r["IF_NEG"] else "N",
                "BMPC_VAL":  r["BMPC_VAL"],
                "BMPC_LT5":  "Y" if r["BMPC_LT5"] else "N",
                "FLC_NORM":  "Y" if r["FLC_NORM"] else "N",
                "RESP_RANK": RESP_RANK.get(r["RSC_IMWG"], 0),
            })
            seq += 1

    return pd.DataFrame(rows)


# ─── ADAM ADRS domain ─────────────────────────────────────────────────────────

def make_adrs(rs: pd.DataFrame, dm: pd.DataFrame, study_key: str) -> pd.DataFrame:
    """
    Generate ADAM ADRS (Analysis Dataset for Response).
    One row per subject: BOR, TTR (time to response), DOR (duration of response).

    Response hierarchy: sCR > CR > VGPR > PR > MR > SD > PD > NE
    ORR (Overall Response Rate) = proportion with BOR ≥ PR
    """
    sid  = f"TOURMALINE-{study_key}"
    rows = []

    for _, subj in dm.iterrows():
        uid  = subj["USUBJID"]
        arm  = subj["ARMCD"]

        subj_rs = rs[rs["USUBJID"] == uid].copy()
        if subj_rs.empty:
            rows.append({
                "STUDYID": sid, "USUBJID": uid, "ARMCD": arm,
                "PARAMCD": "BOR", "PARAM": "Best Overall Response",
                "AVALC": "NE", "AVAL": 0, "TRSPFL": "N", "BORSFL": "Y",
                "TTR": np.nan, "DOR": np.nan, "DOR_CNSR": 1,
            })
            continue

        subj_rs = subj_rs.sort_values("VISITNUM")

        # Best overall response (BOR)
        best_resp  = "NE"
        best_rank  = 0
        best_visit = np.nan

        for _, r in subj_rs.iterrows():
            rk = RESP_RANK.get(r["RSORRES"], 0)
            if rk > best_rank:
                best_rank  = rk
                best_resp  = r["RSORRES"]
                best_visit = r["VISITNUM"]

        # Time to response (TTR): days from first dose to first ≥PR
        ttr = np.nan
        first_pr_visit = np.nan
        for _, r in subj_rs.iterrows():
            if RESP_RANK.get(r["RSORRES"], 0) >= RESP_RANK["PR"]:
                first_pr_visit = r["VISITNUM"]
                ttr = (r["VISITNUM"] - 1) * 28   # days from start
                break

        # Duration of response (DOR): from first ≥PR to PD or censoring
        dor      = np.nan
        dor_cnsr = 1    # 1 = censored, 0 = event (PD or death)
        if not pd.isna(ttr):
            # Find first PD after first PR
            post_pr = subj_rs[subj_rs["VISITNUM"] >= first_pr_visit]
            pd_rows = post_pr[post_pr["RSORRES"] == "PD"]
            if not pd_rows.empty:
                pd_visit = pd_rows["VISITNUM"].iloc[0]
                dor      = (pd_visit - first_pr_visit) * 28
                dor_cnsr = 0
            else:
                # Censored at last assessment
                dor      = (subj_rs["VISITNUM"].max() - first_pr_visit) * 28
                dor_cnsr = 1

        trsp = "Y" if best_rank >= RESP_RANK["PR"] else "N"

        rows.append({
            "STUDYID":    sid,
            "USUBJID":    uid,
            "ARMCD":      arm,
            "PARAMCD":    "BOR",
            "PARAM":      "Best Overall Response",
            "AVALC":      best_resp,
            "AVAL":       best_rank,
            "TRSPFL":     trsp,      # treatment responder (≥PR)
            "BORSFL":     "Y",       # flag for BOR record
            "TTR":        ttr,       # days from first dose to first ≥PR
            "DOR":        dor,       # days of maintained response
            "DOR_CNSR":   dor_cnsr,  # 0=PD event, 1=censored
            "BOR_VISIT":  best_visit,
        })

    adrs = pd.DataFrame(rows)

    # Add response-category flags
    for cat, code in [("SCR","sCR"), ("CR","CR"), ("VGPR","VGPR"),
                      ("PR","PR"), ("MR","MR"), ("SD","SD"), ("PD","PD")]:
        adrs[f"RESP_{cat}FL"] = np.where(adrs["AVALC"] == code, "Y", "N")

    adrs["VGPR_PLUS_FL"] = np.where(adrs["AVAL"] >= RESP_RANK["VGPR"], "Y", "N")
    adrs["CR_PLUS_FL"]   = np.where(adrs["AVAL"] >= RESP_RANK["CR"],   "Y", "N")

    # Convert TTR to months for readability
    adrs["TTR_MO"]  = (adrs["TTR"] / 30.4375).round(1)
    adrs["DOR_MO"]  = (adrs["DOR"] / 30.4375).round(1)

    return adrs


# ─── BMPC LB records ──────────────────────────────────────────────────────────

def make_bmpc_lb(dm: pd.DataFrame, lb: pd.DataFrame, study_key: str) -> pd.DataFrame:
    """
    Simulate BMPC (Bone Marrow Plasma Cell %) records for SDTM LB.
    Assessed at: Baseline (C1), Cycle 3, Cycle 6, Cycle 12, Cycle 24 (EOT).
    ~22% missing rate per assessment.
    Correlated with individual M-protein trajectory.
    """
    sid  = f"TOURMALINE-{study_key}"
    rows = []

    mp_all = lb[lb["LBTESTCD"] == "SPEP_MPROT"].copy()
    mp_all = mp_all.dropna(subset=["LBSTRESN"])

    # Population baseline M-protein for normalisation
    bl_mp_pop = mp_all[mp_all["EPOCH"] == "BASELINE"]["LBSTRESN"].median()
    if pd.isna(bl_mp_pop):
        bl_mp_pop = 15.0

    seq = 1
    for _, subj in dm.iterrows():
        uid = subj["USUBJID"]
        arm = subj["ARMCD"]

        subj_mp = mp_all[mp_all["USUBJID"] == uid].sort_values("VISITNUM")
        if subj_mp.empty:
            continue

        max_vis   = int(subj_mp["VISITNUM"].max())
        bl_mp_row = subj_mp[subj_mp["EPOCH"] == "BASELINE"]
        if bl_mp_row.empty:
            continue
        bl_mp = bl_mp_row["LBSTRESN"].iloc[0]

        # Get M-protein values at BMPC assessment visits
        assessment_visits = sorted([v for v in BMPC_ASSESSMENT_VISITS if v <= max_vis])
        if not assessment_visits:
            continue

        mp_at_assessments = []
        for v in assessment_visits:
            mp_row = subj_mp[subj_mp["VISITNUM"] == v]
            mp_v   = mp_row["LBSTRESN"].iloc[0] if not mp_row.empty else np.nan
            mp_at_assessments.append(mp_v)

        mp_arr  = np.array(mp_at_assessments, dtype=float)
        valid   = ~np.isnan(mp_arr)
        mp_arr[~valid] = bl_mp   # impute missing M-protein with baseline for BMPC sim

        bmpc_arr = _sim_bmpc(bl_mp, mp_arr, bl_mp_pop, arm)

        for i, vis in enumerate(assessment_visits):
            # Simulate missing (bone marrow biopsy not always done)
            if RNG.random() < BMPC_MISSING_RATE:
                continue

            bmpc_v  = bmpc_arr[i]
            epoch   = "BASELINE" if vis == 1 else "TREATMENT"

            rows.append({
                "STUDYID":  sid,
                "DOMAIN":   "LB",
                "USUBJID":  uid,
                "LBSEQ":    seq,
                "LBCAT":    "BONE MARROW",
                "LBTESTCD": "BMPC",
                "LBTEST":   "Bone Marrow Plasma Cells",
                "LBORRES":  f"{bmpc_v:.1f}",
                "LBORRESU": "%",
                "LBSTRESC": f"{bmpc_v:.1f}",
                "LBSTRESN": round(float(bmpc_v), 1),
                "LBSTRESU": "%",
                "LBNRLO":   0.0,
                "LBNRHI":   5.0,
                "LBNRIND":  "LOW" if bmpc_v < 0 else (
                            "HIGH" if bmpc_v > 5.0 else "NORMAL"),
                "VISITNUM": vis,
                "VISIT":    "BASELINE" if vis == 1 else f"CYCLE {vis}",
                "LBDTC":    iso_date((vis - 1) * 28),
                "LBDY":     (vis - 1) * 28 + 1,
                "EPOCH":    epoch,
            })
            seq += 1

    return pd.DataFrame(rows)


# ─── Validation snapshot ──────────────────────────────────────────────────────

def print_response_summary(adrs: pd.DataFrame, study_key: str):
    targets = RESP_TARGETS.get(study_key, {})
    print(f"\n  ── Response summary ({study_key}) ──")
    for arm in ["IRd", "Rd"]:
        sub = adrs[adrs["ARMCD"] == arm]
        n   = len(sub)
        if n == 0:
            continue

        orr     = (sub["TRSPFL"] == "Y").mean() * 100
        vgpr_p  = (sub["VGPR_PLUS_FL"] == "Y").mean() * 100
        cr_p    = (sub["CR_PLUS_FL"] == "Y").mean() * 100
        scr_p   = (sub["RESP_SCRFL"] == "Y").mean() * 100
        med_ttr = sub.loc[sub["TRSPFL"] == "Y", "TTR_MO"].median()
        med_dor = sub.loc[sub["DOR_CNSR"] == 0, "DOR_MO"].median()

        tgt     = targets.get(arm, {})
        print(f"\n     {arm} (N={n}):")
        print(f"       ORR  (≥PR): {orr:.1f}%   "
              f"[published {tgt.get('ORR','-')}%]")
        print(f"       VGPR+    : {vgpr_p:.1f}%   "
              f"[published {tgt.get('VGPR_plus','-')}%]")
        print(f"       CR+      : {cr_p:.1f}%   "
              f"[published {tgt.get('CR_plus','-')}%]")
        print(f"       sCR      : {scr_p:.1f}%   "
              f"[published {tgt.get('sCR','-')}%]")
        print(f"       Median TTR: {med_ttr:.1f} mo  "
              f"Median DOR (events): {med_dor:.1f} mo")


# ─── MAIN ─────────────────────────────────────────────────────────────────────

def generate_pd(study_key: str):
    sdir = os.path.join(BASE, study_key)
    print(f"\n{'='*60}")
    print(f"  Generating PD domains — {study_key}")
    print(f"{'='*60}")

    dm   = pd.read_csv(f"{sdir}/sdtm_dm.csv")
    lb   = pd.read_csv(f"{sdir}/sdtm_lb.csv")

    # ── 1. BMPC records ───────────────────────────────────────────────────────
    print("  [1/4] BMPC — Bone Marrow Plasma Cell % ...")
    bmpc_lb = make_bmpc_lb(dm, lb, study_key)
    print(f"        {len(bmpc_lb):,} BMPC records "
          f"({bmpc_lb['USUBJID'].nunique()} subjects)")

    # Append to sdtm_lb.csv
    lb_new = pd.concat([lb, bmpc_lb], ignore_index=True)
    lb_new.to_csv(f"{sdir}/sdtm_lb.csv", index=False)

    # Rebuild adam_adlb
    adsl = pd.read_csv(f"{sdir}/adam_adsl.csv")
    cols_adsl = ["USUBJID","ARMCD","TRT01PN","IGTYPE","ISSSTAGE","RISS","CYTOGR"]
    cols_adsl = [c for c in cols_adsl if c in adsl.columns]

    adlb = lb_new.copy()
    adlb = adlb.merge(adsl[cols_adsl], on="USUBJID", how="left")
    baseline_lb = lb_new[
        (lb_new["EPOCH"] == "BASELINE") & lb_new["LBSTRESN"].notna()
    ][["USUBJID","LBTESTCD","LBSTRESN"]].rename(columns={"LBSTRESN": "BASE"})
    adlb = adlb.merge(baseline_lb, on=["USUBJID","LBTESTCD"], how="left")
    adlb["CHG"]     = adlb["LBSTRESN"] - adlb["BASE"]
    adlb["PCHG"]    = np.where(adlb["BASE"].fillna(0) != 0,
                               adlb["CHG"] / adlb["BASE"] * 100, np.nan)
    adlb["AVAL"]    = adlb["LBSTRESN"]
    adlb["AVALU"]   = adlb["LBSTRESU"]
    adlb["PARAM"]   = adlb["LBTEST"]
    adlb["PARAMCD"] = adlb["LBTESTCD"]
    adlb["ANL01FL"] = "Y"
    adlb["ABLFL"]   = np.where(adlb["EPOCH"] == "BASELINE", "Y", "")
    adlb["DTYPE"]   = ""
    adlb.to_csv(f"{sdir}/adam_adlb.csv", index=False)
    print(f"        sdtm_lb.csv + adam_adlb.csv updated")

    # ── 2. RS — IMWG response per cycle ──────────────────────────────────────
    print("  [2/4] RS  — IMWG response per cycle ...")
    rs = make_rs(dm, lb_new, bmpc_lb, study_key)
    rs.to_csv(f"{sdir}/sdtm_rs.csv", index=False)
    print(f"        {len(rs):,} RS records ({rs['USUBJID'].nunique()} subjects)")

    # ── 3. ADRS — Best overall response summary ───────────────────────────────
    print("  [3/4] ADRS — Best Overall Response + TTR + DOR ...")
    adrs = make_adrs(rs, dm, study_key)
    adrs.to_csv(f"{sdir}/adam_adrs.csv", index=False)
    print(f"        {len(adrs):,} ADRS records")

    # ── 4. Summary ────────────────────────────────────────────────────────────
    print("  [4/4] Validation snapshot ...")
    print_response_summary(adrs, study_key)

    # BOR distribution
    bor_dist = adrs.groupby(["ARMCD","AVALC"]).size().unstack(fill_value=0)
    print(f"\n  BOR distribution:\n{bor_dist.to_string()}")

    print(f"\n  Files written to {sdir}/")
    print(f"    sdtm_lb.csv   (BMPC appended)")
    print(f"    adam_adlb.csv (BMPC appended)")
    print(f"    sdtm_rs.csv   (NEW)")
    print(f"    adam_adrs.csv (NEW)")


def main():
    for sk in ["MM2", "MM1"]:
        path = os.path.join(BASE, sk)
        if not os.path.exists(path):
            print(f"[SKIP] {sk} not found")
            continue
        generate_pd(sk)

    print(f"\n{'='*60}")
    print("  PD generation complete.")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
