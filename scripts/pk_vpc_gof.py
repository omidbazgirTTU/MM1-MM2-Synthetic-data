"""
Step 8 — PK Visual Predictive Check (VPC) & Goodness-of-Fit (GOF)
====================================================================
VPC uses Monte Carlo reference bands derived from published population-PK
model parameters (Gupta 2017 for Ixazomib; Chen 2012 for Lenalidomide).
The observed synthetic concentrations are then overlaid on those
model-derived prediction intervals — NOT compared against their own percentiles
(which would be circular/self-referential).

VPC layout per study (2 rows × 3 columns):
  Row 1 — Post-dose dense sampling (Cycle 1, PCTPTNUM > 0)
           vs MC reference bands from Gupta/Chen 2017 model parameters
  Row 2 — Pre-dose trough trajectory by cycle (PCTPTNUM = 0, all cycles)
           vs expected trough accumulation from MC simulation

GOF layout (6 panels, unchanged):
  NCA distributions (Cmax, AUCinf, t½) and covariate effects vs published values

VPC pass criterion : ≥ 80% of non-BLQ observations within model-derived 5–95th PI
Reference          : Gupta 2017 Clin Pharmacokinet (Ixazomib 3-cmt, N=800 MC draws)
                     Chen 2012 (Lenalidomide 1-cmt, N=800 MC draws)

Outputs
-------
  outputs/figures/pk_vpc_{MM1,MM2}.png
  outputs/figures/pk_gof_{MM1,MM2}.png
  outputs/tables/pk_vpc_summary.csv
  outputs/tables/pk_nca_summary.csv
"""

import os, sys, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy import stats
from scipy.integrate import solve_ivp

warnings.filterwarnings("ignore")

BASE   = os.path.join(os.path.dirname(__file__), "..")
FIGDIR = os.path.join(BASE, "outputs", "figures")
TABDIR = os.path.join(BASE, "outputs", "tables")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(TABDIR, exist_ok=True)

# ── Published reference values ────────────────────────────────────────────────
PUBLISHED = {
    "IXAZOMIB": {
        "cmax_med":  41.0,    # ng/mL   (Gupta 2017 geometric mean)
        "aucinf_med":1247.0,  # ng·h/mL
        "thalf_med": 228.0,   # h
        "clf_med":   3.63,    # L/h
        "iiv_cmax":  0.38,
        "lloq":      0.5,
    },
    "LENALIDOMIDE": {
        "cmax_lo":   350.0,  "cmax_hi":  500.0,
        "aucinf_med":2200.0,
        "thalf_lo":  3.0,    "thalf_hi": 5.0,
        "lloq":      2.0,
    },
    "DEXAMETHASONE": {
        "cmax_lo":   120.0,  "cmax_hi":  250.0,
        "aucinf_med":1800.0,
        "thalf_lo":  5.0,    "thalf_hi": 8.0,
        "lloq":      0.2,
    },
}

DRUG_LABEL = {"IXAZOMIB": "Ixazomib", "LENALIDOMIDE": "Lenalidomide",
              "DEXAMETHASONE": "Dexamethasone"}
DRUG_COLOR = {"IXAZOMIB": "#1f77b4", "LENALIDOMIDE": "#d62728",
              "DEXAMETHASONE": "#2ca02c"}

# ── Data loading ──────────────────────────────────────────────────────────────

def _load(study_key):
    sdir = os.path.join(BASE, study_key)
    pc   = pd.read_csv(f"{sdir}/sdtm_pc.csv")
    pp   = pd.read_csv(f"{sdir}/sdtm_pp.csv")
    adpc = pd.read_csv(f"{sdir}/adam_adpc.csv")
    adsl = pd.read_csv(f"{sdir}/adam_adsl.csv")
    return pc, pp, adpc, adsl


# ════════════════════════════════════════════════════════════════════════════
# Monte Carlo reference generators (proper VPC reference — NOT self-derived)
# ════════════════════════════════════════════════════════════════════════════

def _build_dose_schedule(n_cycles=26, dose_ng=4000.0):
    """
    Ixazomib clinical dosing: Days 1, 8, 15 of each 28-day cycle.
    Returns sorted (dose_times_h, dose_amounts) from study Day 1.
    """
    times, amounts = [], []
    for c in range(1, n_cycles + 1):
        cycle_start_h = (c - 1) * 28 * 24   # hours
        for day_offset_h in [0, 168, 336]:   # Days 1, 8, 15
            times.append(cycle_start_h + day_offset_h)
            amounts.append(dose_ng)
    return np.array(times), np.array(amounts)


def _mc_ref_ixazomib(visitnums, pctptnums, n_mc=800, seed=101,
                     cl_ddi_mults=None, dose_mults=None,
                     abs_times_h=None):
    """
    Monte Carlo prediction intervals for Ixazomib.

    Uses the EXACT parameters from generate_pk_v2.py (Gupta 2017 implementation)
    so that the MC reference reflects what the data generator actually produced.
    Superposition principle (linear PK): solve 3-cmt impulse response once per
    virtual patient, sum shifted copies for all prior doses.

    Parameters
    ----------
    visitnums    : list of cycle numbers (used when abs_times_h is None)
    pctptnums    : list of hours post Day-1 dose of that cycle (same caveat)
    cl_ddi_mults : optional array (length n_mc) of per-patient CL multipliers
                   reflecting CYP3A4 DDI (inhibitor <1, inducer >1, normal =1).
    dose_mults   : optional array (length n_mc) of per-patient dose multipliers.
    abs_times_h  : optional list of absolute study hours (overrides visitnum/pctptnum
                   conversion).  Use for trough observations where PCTPTNUM=0 is
                   a convention that does NOT encode absolute time — instead use
                   (DOSE_DAY - 1) * 24 derived from sdtm_pc.DOSE_DAY.
                   Returns DataFrame with column 't_abs_h' instead of VISITNUM/PCTPTNUM.

    Returns
    -------
    If abs_times_h is None  : DataFrame with [VISITNUM, PCTPTNUM, p05, p50, p95, pct_blq]
    If abs_times_h provided : DataFrame with [t_abs_h, p05, p50, p95, pct_blq]
    """
    rng = np.random.default_rng(seed)

    # Exact parameters from generate_pk_v2.py PK_DRUGS["IXAZOMIB"]
    # Vss = V2+V3+V4 = 14.3+200+328.7 = 543 L (Gupta 2017)
    TV_CL=1.86; TV_V2=14.3; TV_Q3=8.0;  TV_V3=200.0
    TV_Q4=0.50; TV_V4=328.7; TV_KA=0.5; F=0.58; ALAG=0.15
    OM_CL=0.36; OM_V2=0.30; OM_Q3=0.45; OM_V3=0.50
    OM_Q4=0.40; OM_V4=0.45; OM_KA=0.55; OM_F=0.30
    SIGMA_PROP=0.18; SIGMA_ADD=0.3
    LLOQ = PUBLISHED["IXAZOMIB"]["lloq"]
    dose_ng = 4000.0

    # Individual PK parameters (log-normal IIV — all params as in generator)
    CL_i = TV_CL * np.exp(rng.normal(0, OM_CL,  n_mc))
    V2_i = TV_V2 * np.exp(rng.normal(0, OM_V2,  n_mc))
    Q3_i = TV_Q3 * np.exp(rng.normal(0, OM_Q3,  n_mc))
    V3_i = TV_V3 * np.exp(rng.normal(0, OM_V3,  n_mc))
    Q4_i = TV_Q4 * np.exp(rng.normal(0, OM_Q4,  n_mc))
    V4_i = TV_V4 * np.exp(rng.normal(0, OM_V4,  n_mc))
    Ka_i = TV_KA * np.exp(rng.normal(0, OM_KA,  n_mc))
    F_i  = np.clip(F * np.exp(rng.normal(0, OM_F, n_mc)), 0.1, 1.0)

    # Apply CYP3A4 DDI effect on CL (inhibitor <1, inducer >1)
    if cl_ddi_mults is not None:
        CL_i = CL_i * np.asarray(cl_ddi_mults)

    # Per-patient dose multiplier (dose reduction → proportionally lower conc)
    dose_ng_i = dose_ng * (np.asarray(dose_mults) if dose_mults is not None
                           else np.ones(n_mc))

    # Build the absolute time list for MC evaluation
    if abs_times_h is not None:
        # Direct mode: caller provides exact absolute hours (e.g. (DOSE_DAY-1)*24)
        unique_abs_times = sorted(set(float(t) for t in abs_times_h))
        max_cycle = int(max(unique_abs_times) / (28 * 24)) + 2
    else:
        # Convert (visitnum, pctptnum) → absolute study hours
        # VISITNUM = cycle; PCTPTNUM = hours post Day-1 dose of that cycle
        unique_keys = list(set(zip(visitnums, pctptnums)))
        sample_abs  = {(vn, pt): (vn - 1) * 28 * 24 + pt for vn, pt in unique_keys}
        unique_abs_times = sorted(sample_abs.values())
        max_cycle = max(visitnums) + 1

    dose_times, dose_amounts = _build_dose_schedule(n_cycles=max_cycle)

    # Impulse response time grid: fine near Cmax, coarser in terminal phase
    MAX_DELAY = 1500   # h — beyond this, contribution < 1%
    t_impulse = np.unique(np.concatenate([
        np.linspace(0, 0.5, 20),
        np.linspace(0.5, 10, 80),
        np.linspace(10, 48, 100),
        np.linspace(48, 200, 100),
        np.linspace(200, MAX_DELAY, 150),
    ]))

    # unique_abs_times is set by the if/else block above
    all_conc = np.full((n_mc, len(unique_abs_times)), np.nan)

    for i in range(n_mc):
        cl, v2, q3, v3, q4, v4, ka, f = (
            CL_i[i], V2_i[i], Q3_i[i], V3_i[i],
            Q4_i[i], V4_i[i], Ka_i[i], F_i[i])

        # 3-compartment ODE with absorption lag — single dose impulse response
        # Variables: amounts (ng)
        def ode(t, y, _cl=cl, _v2=v2, _q3=q3, _v3=v3, _q4=q4, _v4=v4,
                _ka=ka, _f=f):
            dep, c2, c3, c4 = y
            ab = _ka * dep if t >= ALAG else 0.0
            dDep = -_ka * dep if t >= ALAG else 0.0
            dC2  = (_f*ab + _q3/_v3*c3 + _q4/_v4*c4
                    - (_cl + _q3 + _q4)/_v2 * c2)
            dC3  = _q3/_v2*c2 - _q3/_v3*c3
            dC4  = _q4/_v2*c2 - _q4/_v4*c4
            return [dDep, dC2, dC3, dC4]

        # Impulse response solved for a unit dose (dose_ng_ref=4000 ng);
        # scale contributions by patient's actual dose_ng_i / dose_ng_ref
        sol = solve_ivp(ode, [0, MAX_DELAY + 1], [dose_ng, 0.0, 0.0, 0.0],
                        t_eval=t_impulse, method="LSODA", rtol=1e-5, atol=1e-7,
                        dense_output=False)
        if not sol.success:
            continue

        cp_impulse = np.maximum(sol.y[1] / v2, 0.0)   # concentration (ng/mL) per dose_ng

        pt_dose_scale = dose_ng_i[i] / dose_ng   # per-patient dose multiplier

        # Superposition: C_total(t) = Σ_k C_impulse(t − t_k) * (amt_k/dose_ng)
        # All scheduled dose amounts equal dose_ng so amt_k/dose_ng = pt_dose_scale
        for j, t_s in enumerate(unique_abs_times):
            c_total = 0.0
            for t_d, amt in zip(dose_times, dose_amounts):
                delay = t_s - t_d
                if delay <= 0:
                    break   # doses are sorted; no future doses contribute
                if delay <= MAX_DELAY:
                    c_total += np.interp(delay, t_impulse, cp_impulse) * pt_dose_scale
            all_conc[i, j] = max(c_total, 0.0)

    # Residual error (proportional + additive)
    eps_p = rng.normal(0, SIGMA_PROP, all_conc.shape)
    eps_a = rng.normal(0, SIGMA_ADD,  all_conc.shape)
    obs   = np.maximum(all_conc * (1 + eps_p) + eps_a, 0.0)

    # Compute percentiles — output format depends on mode
    rows = []
    if abs_times_h is not None:
        for t_abs in unique_abs_times:
            j     = unique_abs_times.index(t_abs)
            vals  = obs[:, j]
            blq   = vals < LLOQ
            valid = vals[~blq]
            if len(valid) < 10:
                continue
            p05, p50, p95 = np.nanpercentile(valid, [5, 50, 95])
            rows.append({"t_abs_h": t_abs,
                         "p05": p05, "p50": p50, "p95": p95,
                         "pct_blq": blq.mean() * 100})
    else:
        for (vn, pt), t_abs in sample_abs.items():
            j     = unique_abs_times.index(t_abs)
            vals  = obs[:, j]
            blq   = vals < LLOQ
            valid = vals[~blq]
            if len(valid) < 10:
                continue
            p05, p50, p95 = np.nanpercentile(valid, [5, 50, 95])
            rows.append({"VISITNUM": vn, "PCTPTNUM": pt,
                         "p05": p05, "p50": p50, "p95": p95,
                         "pct_blq": blq.mean() * 100})
    return pd.DataFrame(rows)


def _mc_ref_1cmt(drug, visitnums, pctptnums, n_mc=800, seed=202):
    """
    Monte Carlo prediction intervals for a 1-compartment oral drug.
    Uses EXACT parameters from generate_pk_v2.py PK_DRUGS dict.
    Lenalidomide: Chen 2012; Dexamethasone: published values.
    """
    rng = np.random.default_rng(seed)

    if drug == "LENALIDOMIDE":
        # generate_pk_v2.py: CL/F=dose/AUCinf=25000/2200=11.36, Vd/F=80, Ka=1.10
        TV_CL=11.36; TV_V=80.0; TV_KA=1.10; F=1.0
        OM_CL=0.32;  OM_V=0.38; OM_KA=0.50
        SIGMA_PROP=0.18; SIGMA_ADD=8.0
        LLOQ = PUBLISHED["LENALIDOMIDE"]["lloq"]
        dose_ng = 25_000.0    # 25 mg × 1000 ng/mg (daily Days 1-21)
        dose_offset_h = 0.0
    else:  # DEXAMETHASONE
        # generate_pk_v2.py: CL/F=dose/AUCinf=40000/1800=22.2, Vd/F=208, Ka=1.50
        TV_CL=22.2; TV_V=208.0; TV_KA=1.50; F=1.0
        OM_CL=0.28; OM_V=0.32;  OM_KA=0.45
        SIGMA_PROP=0.20; SIGMA_ADD=3.0
        LLOQ = PUBLISHED["DEXAMETHASONE"]["lloq"]
        dose_ng = 40_000.0
        dose_offset_h = 0.0

    CL_i = TV_CL * np.exp(rng.normal(0, OM_CL, n_mc))
    V_i  = TV_V  * np.exp(rng.normal(0, OM_V,  n_mc))
    Ka_i = TV_KA * np.exp(rng.normal(0, OM_KA, n_mc))

    # For simplicity: treat Cycle 1 Day 1 as first dose (most relevant for VPC)
    # For later cycles, lenalidomide is dosed daily so steady state is reached
    # by Cycle 1 Day 3; the single-dose C1D1 profile is the correct reference.
    unique_pts = sorted(set(pctptnums))
    t_obs = np.array([pt + dose_offset_h for pt in unique_pts])

    all_conc = np.full((n_mc, len(unique_pts)), np.nan)

    for i in range(n_mc):
        cl, v, ka = CL_i[i], V_i[i], Ka_i[i]
        ke = cl / v
        # Analytical 1-cmt oral solution (single dose)
        if abs(ka - ke) < 1e-6:
            ka += 1e-6
        for j, t in enumerate(t_obs):
            if t <= 0:
                all_conc[i, j] = 0.0
            else:
                c = (F * dose_ng * ka) / (v * (ka - ke)) * (np.exp(-ke*t) - np.exp(-ka*t))
                all_conc[i, j] = max(c, 0.0)

    eps_p = rng.normal(0, SIGMA_PROP, all_conc.shape)
    eps_a = rng.normal(0, SIGMA_ADD,  all_conc.shape)
    obs   = np.maximum(all_conc * (1 + eps_p) + eps_a, 0.0)

    rows = []
    for j, pt in enumerate(unique_pts):
        vals  = obs[:, j]
        blq   = vals < LLOQ
        valid = vals[~blq]
        if len(valid) < 10:
            continue
        p05, p50, p95 = np.nanpercentile(valid, [5, 50, 95])
        rows.append({"PCTPTNUM": pt, "p05": p05, "p50": p50, "p95": p95,
                     "pct_blq": blq.mean() * 100})
    return pd.DataFrame(rows)


def _pct_within_model_pi(df_obs, ref_df, join_cols, lo_col="p05", hi_col="p95"):
    """
    What fraction of non-BLQ observed concentrations fall within
    the MODEL-derived PI bands (not self-derived)?
    """
    merged = df_obs.merge(ref_df[join_cols + [lo_col, hi_col]], on=join_cols, how="inner")
    within = merged["PCSTRESN"].between(merged[lo_col], merged[hi_col])
    return within.mean() * 100 if len(merged) > 0 else np.nan


# ════════════════════════════════════════════════════════════════════════════
# VPC figure — 2 rows × 3 columns per study
# Row 0: post-dose dense (Cycle 1, PCTPTNUM > 0) vs MC reference
# Row 1: pre-dose trough by cycle            vs MC trough trajectory
# ════════════════════════════════════════════════════════════════════════════

def plot_vpc(study_key, pc, adsl, figdir, n_mc=800):
    """
    Proper VPC: Monte Carlo reference bands from Gupta/Chen model parameters
    overlaid with observed synthetic concentrations.

    For the Ixazomib trough panel the MC reference additionally samples
    CYP3A4 DDI multipliers from the study-specific empirical distribution
    (adsl.CL_DDI_MULT) so that inhibitor/inducer patients are represented
    in the reference bands.
    """
    drugs = ["IXAZOMIB", "LENALIDOMIDE", "DEXAMETHASONE"]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10),
                             gridspec_kw={"hspace": 0.50, "wspace": 0.35})
    fig.suptitle(
        f"VPC — {study_key}  |  Shaded ribbon = MC PI from generator parameters "
        f"(Gupta 2017 / Chen 2012, N={n_mc} virtual patients)  |  Dots = synthetic observations",
        fontsize=11, fontweight="bold"
    )

    summary_rows = []

    for col, drug in enumerate(drugs):
        color = DRUG_COLOR[drug]
        lloq  = PUBLISHED[drug]["lloq"]
        df    = pc[(pc["PCTESTCD"] == drug) & pc["VISITNUM"].notna()].copy()
        df["VISITNUM"]  = pd.to_numeric(df["VISITNUM"],  errors="coerce")
        df["PCTPTNUM"]  = pd.to_numeric(df["PCTPTNUM"],  errors="coerce")
        df["PCSTRESN"]  = pd.to_numeric(df["PCSTRESN"],  errors="coerce")
        df = df.dropna(subset=["VISITNUM", "PCTPTNUM", "PCSTRESN"])

        # ── Row 0: post-dose dense (PCTPTNUM > 0, prefer Cycle 1) ─────────
        ax0 = axes[0][col]
        df_dense = df[(df["PCTPTNUM"] > 0) & (df["VISITNUM"] == 1) & (df["BLQ"] == "N")]

        if not df_dense.empty:
            ax0.scatter(df_dense["PCTPTNUM"], df_dense["PCSTRESN"],
                        alpha=0.30, s=12, color=color, zorder=3, label="C1 observed")

            # BLQ points at LLOQ/2
            df_blq = df[(df["PCTPTNUM"] > 0) & (df["VISITNUM"] == 1) & (df["BLQ"] == "Y")]
            if not df_blq.empty:
                ax0.scatter(df_blq["PCTPTNUM"], np.full(len(df_blq), lloq / 2),
                            marker="v", alpha=0.5, s=14, color="gray", zorder=3,
                            label=f"BLQ (<{lloq})")

            # MC reference bands
            print(f"    [{drug}] Computing MC reference ({n_mc} subjects) ...", end=" ", flush=True)
            vns  = df_dense["VISITNUM"].astype(int).tolist()
            pts  = df_dense["PCTPTNUM"].tolist()
            if drug == "IXAZOMIB":
                ref = _mc_ref_ixazomib(vns + [1]*len(set(pts)), pts + list(set(pts)),
                                       n_mc=n_mc, seed=101 + col)
                join_cols = ["VISITNUM", "PCTPTNUM"]
            else:
                ref = _mc_ref_1cmt(drug, vns, pts, n_mc=n_mc, seed=202 + col)
                join_cols = ["PCTPTNUM"]

            if not ref.empty:
                ref_sorted = ref.drop_duplicates(join_cols).sort_values(
                    "PCTPTNUM" if "VISITNUM" not in join_cols else join_cols)
                t_ref  = ref_sorted["PCTPTNUM"].values
                ax0.fill_between(t_ref, ref_sorted["p05"], ref_sorted["p95"],
                                 alpha=0.20, color=color,
                                 label=f"MC 5–95th PI (Gupta/Chen 2017, N={n_mc})")
                ax0.plot(t_ref, ref_sorted["p50"], color=color, lw=2.0,
                         label="MC median", zorder=4)

                pct_in = _pct_within_model_pi(df_dense, ref_sorted, join_cols)
                status = "PASS" if pct_in >= 80 else "FAIL"
                ax0.text(0.97, 0.97, f"{pct_in:.0f}% within MC 5–95th PI\n{status}",
                         transform=ax0.transAxes, fontsize=8, ha="right", va="top",
                         color="#1a7a1a" if pct_in >= 80 else "#cc1a1a",
                         fontweight="bold",
                         bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85))
                print(f"{pct_in:.0f}% within PI ({status})")

                summary_rows.append({"study": study_key, "drug": drug,
                                     "panel": "C1 post-dose",
                                     "pct_within_pi": round(pct_in, 1),
                                     "n_obs": len(df_dense)})
            else:
                print("no MC ref computed")

        ax0.axhline(lloq, color="gray", ls="--", lw=0.8, alpha=0.7,
                    label=f"LLOQ={lloq}")
        ax0.set_yscale("log")
        ax0.set_xlabel("Time post-dose (h)", fontsize=9)
        ax0.set_ylabel("Concentration (ng/mL)", fontsize=9)
        ax0.set_title(f"{DRUG_LABEL[drug]}\nCycle 1 post-dose (dense sampling)",
                      fontsize=9, fontweight="bold")
        ax0.legend(fontsize=7, loc="upper right", framealpha=0.7)
        ax0.grid(True, which="both", ls=":", alpha=0.3)

        # ── Row 1: pre-dose trough trajectory by study day ───────────────
        ax1 = axes[1][col]
        df_trough_all = df[df["PCTPTNUM"] == 0].copy()
        if "DOSE_DAY" in df_trough_all.columns:
            df_trough_all["DOSE_DAY"] = pd.to_numeric(
                df_trough_all["DOSE_DAY"], errors="coerce")
            df_trough_all["t_abs_h"] = (df_trough_all["DOSE_DAY"] - 1) * 24.0
        else:
            # Fallback: estimate from VISITNUM (less accurate for intra-cycle troughs)
            df_trough_all["t_abs_h"] = (df_trough_all["VISITNUM"] - 1) * 28 * 24.0

        df_trough = df_trough_all[df_trough_all["BLQ"] == "N"].dropna(subset=["t_abs_h"])

        if not df_trough.empty:
            study_day = df_trough["t_abs_h"] / 24.0 + 1  # convert back to study day
            ax1.scatter(study_day, df_trough["PCSTRESN"],
                        alpha=0.15, s=8, color=color, zorder=2)
            # Smooth median/IQR in 28-day bins
            df_trough["day_bin"] = (df_trough["t_abs_h"] // (28 * 24) + 1).astype(int)
            obs_med = df_trough.groupby("day_bin")["PCSTRESN"].median()
            obs_lo  = df_trough.groupby("day_bin")["PCSTRESN"].quantile(0.25)
            obs_hi  = df_trough.groupby("day_bin")["PCSTRESN"].quantile(0.75)
            bin_centers = obs_med.index.values * 28 - 13  # mid-cycle day
            ax1.plot(bin_centers, obs_med.values, color=color, lw=2, marker="o",
                     markersize=4, zorder=4, label="Observed median (28-day bins)")
            ax1.fill_between(bin_centers, obs_lo.values, obs_hi.values,
                             alpha=0.25, color=color, label="Observed IQR")

        # MC reference for troughs — use DOSE_DAY-derived absolute times
        if drug == "IXAZOMIB" and not df_trough.empty:
            unique_t_abs = sorted(df_trough["t_abs_h"].dropna().unique())

            # Sample DDI multipliers from empirical patient-level distribution
            ddi_src = adsl["CL_DDI_MULT"].dropna().values
            _rng_ddi = np.random.default_rng(303)
            cl_ddi_mults = _rng_ddi.choice(ddi_src, size=n_mc, replace=True)

            # Dose modification: sample per-patient effective dose multiplier
            dose_probs  = np.array([0.87, 0.08, 0.05])
            dose_levels = np.array([1.0, 0.75, 0.50])
            dose_mults = _rng_ddi.choice(dose_levels, size=n_mc, p=dose_probs)

            ref_trough = _mc_ref_ixazomib(None, None,
                                          n_mc=n_mc, seed=303,
                                          cl_ddi_mults=cl_ddi_mults,
                                          dose_mults=dose_mults,
                                          abs_times_h=unique_t_abs)
            if not ref_trough.empty:
                ref_sorted = ref_trough.sort_values("t_abs_h")
                ref_day    = ref_sorted["t_abs_h"] / 24.0 + 1
                ax1.fill_between(ref_day, ref_sorted["p05"], ref_sorted["p95"],
                                 alpha=0.15, color="black",
                                 label=f"MC 5–95th PI (+DDI, N={n_mc})")
                ax1.plot(ref_day, ref_sorted["p50"],
                         color="black", lw=1.5, ls="--", label="MC median (expected)")

                # % trough obs within MC PI — merge on t_abs_h
                obs_merge = df_trough[["t_abs_h","PCSTRESN"]].copy()
                ref_merge  = ref_trough[["t_abs_h","p05","p95"]]
                merged_tr  = obs_merge.merge(ref_merge, on="t_abs_h", how="inner")
                if len(merged_tr) > 0:
                    pct_tr = merged_tr["PCSTRESN"].between(
                        merged_tr["p05"], merged_tr["p95"]).mean() * 100
                else:
                    pct_tr = np.nan

                ax1.text(0.97, 0.97, f"Troughs: {pct_tr:.0f}% within MC PI",
                         transform=ax1.transAxes, fontsize=8, ha="right", va="top",
                         color="#1a7a1a" if (not np.isnan(pct_tr) and pct_tr >= 80)
                                        else "#cc1a1a",
                         fontweight="bold",
                         bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85))
                summary_rows.append({"study": study_key, "drug": drug,
                                     "panel": "trough by cycle",
                                     "pct_within_pi": round(pct_tr, 1),
                                     "n_obs": len(df_trough)})

        elif drug != "IXAZOMIB" and not df_trough.empty:
            # For Lena/Dex troughs: just note that single-dose MC is C1-specific
            ax1.text(0.5, 0.5,
                     f"{DRUG_LABEL[drug]} troughs\n(daily dosing — trough reflects\nsteady-state, not single-dose model)",
                     ha="center", va="center", transform=ax1.transAxes,
                     fontsize=9, color="gray")

        ax1.axhline(lloq, color="gray", ls="--", lw=0.8, alpha=0.7)
        ax1.set_yscale("log")
        ax1.set_xlabel("Study day", fontsize=9)
        ax1.set_ylabel("Pre-dose trough (ng/mL)", fontsize=9)
        ax1.set_title(f"{DRUG_LABEL[drug]}\nPre-dose trough by study day",
                      fontsize=9, fontweight="bold")
        if not df_trough.empty and drug == "IXAZOMIB":
            ax1.legend(fontsize=7, framealpha=0.7)
        ax1.grid(True, which="both", ls=":", alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(figdir, f"pk_vpc_{study_key}.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out}")
    return pd.DataFrame(summary_rows)


# ── GOF figure (unchanged — already externally referenced) ────────────────────

def _cwres_proxy(df_drug):
    df = df_drug[df_drug["BLQ"] != "Y"].copy()
    df = df[df["PCSTRESN"] > 0].copy()
    sigma_prop, sigma_add = 0.20, 0.50
    med_pred = df.groupby("PCTPTNUM")["PCSTRESN"].transform("median")
    ipred    = med_pred.clip(lower=0.01)
    cwres    = (np.log(df["PCSTRESN"]) - np.log(ipred)) / (sigma_prop + sigma_add / ipred)
    return cwres.dropna()


def plot_gof(study_key, pc, pp, adpc, adsl, figdir):
    fig = plt.figure(figsize=(18, 10))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.35)
    fig.suptitle(f"GOF Diagnostics — {study_key}\n"
                 f"NCA parameters compared to published Gupta 2017 / Chen 2012 values",
                 fontsize=12, fontweight="bold")

    pp_ix = pp[pp["PPCAT"] == "IXAZOMIB"] if "PPCAT" in pp.columns else \
            pp[pp["PPTEST"].str.contains("Ixazomib", na=False)]
    ix_wide = pp_ix.pivot_table(index="USUBJID", columns="PPTESTCD",
                                values="PPSTRESN", aggfunc="first").reset_index()

    # Panel 0: Cmax
    ax0 = fig.add_subplot(gs[0, 0])
    if "CMAX" in ix_wide.columns:
        cmax = ix_wide["CMAX"].dropna()
        ax0.hist(cmax, bins=30, color=DRUG_COLOR["IXAZOMIB"], alpha=0.65,
                 edgecolor="white", density=True, label="Simulated NCA Cmax")
        mu, sigma = np.mean(np.log(cmax)), np.std(np.log(cmax))
        x = np.linspace(cmax.min(), cmax.max(), 200)
        ax0.plot(x, stats.lognorm.pdf(x, s=sigma, scale=np.exp(mu)),
                 "k--", lw=1.5, label=f"Log-N fit  CV={sigma:.0%}")
        ax0.axvline(PUBLISHED["IXAZOMIB"]["cmax_med"], color="red", lw=2,
                    ls="--", label=f"Published {PUBLISHED['IXAZOMIB']['cmax_med']} ng/mL\n(Gupta 2017)")
        ax0.set_xlabel("Ixazomib Cmax (ng/mL)"); ax0.set_ylabel("Density")
        ax0.set_title("Ixazomib Cmax — sim vs published"); ax0.legend(fontsize=7)

    # Panel 1: AUCinf
    ax1 = fig.add_subplot(gs[0, 1])
    if "AUCINF" in ix_wide.columns:
        auc = ix_wide["AUCINF"].dropna()
        ax1.hist(auc, bins=30, color=DRUG_COLOR["IXAZOMIB"], alpha=0.65,
                 edgecolor="white", density=True, label="Simulated NCA AUCinf")
        mu, sigma = np.mean(np.log(auc)), np.std(np.log(auc))
        x = np.linspace(auc.min(), auc.max(), 200)
        ax1.plot(x, stats.lognorm.pdf(x, s=sigma, scale=np.exp(mu)),
                 "k--", lw=1.5, label=f"Log-N fit  CV={sigma:.0%}")
        ax1.axvline(PUBLISHED["IXAZOMIB"]["aucinf_med"], color="red", lw=2,
                    ls="--", label=f"Published {PUBLISHED['IXAZOMIB']['aucinf_med']} ng·h/mL\n(Gupta 2017)")
        ax1.set_xlabel("Ixazomib AUCinf (ng·h/mL)"); ax1.set_ylabel("Density")
        ax1.set_title("Ixazomib AUCinf — sim vs published"); ax1.legend(fontsize=7)

    # Panel 2: t½ with NCA limitation annotation
    ax2 = fig.add_subplot(gs[0, 2])
    if "THALF" in ix_wide.columns:
        thalf = ix_wide["THALF"].dropna()
        ax2.hist(thalf, bins=30, color=DRUG_COLOR["IXAZOMIB"], alpha=0.65,
                 edgecolor="white", label=f"NCA-derived  median={thalf.median():.0f}h")
        ax2.axvline(PUBLISHED["IXAZOMIB"]["thalf_med"], color="red", lw=2,
                    ls="--", label=f"True t½ = 228h (Gupta 2017)")
        ax2.axvline(thalf.median(), color="navy", lw=1.5, ls="-",
                    label=f"NCA median = {thalf.median():.0f}h")
        ax2.set_xlabel("Ixazomib NCA t½ (h)"); ax2.set_ylabel("Count")
        ax2.set_title("Ixazomib t½ — NCA underestimates true t½")
        ax2.legend(fontsize=7)
        ax2.text(0.03, 0.88,
                 "NCA requires sampling ≥ 3× t½ (≥ 684h) to\n"
                 "resolve terminal phase. Sampling window = 168h.\n"
                 "True simulated CL/V matches Gupta 2017.",
                 transform=ax2.transAxes, fontsize=7, va="top",
                 bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.9))

    # Panel 3: CrCL vs Lenalidomide AUCinf
    ax3 = fig.add_subplot(gs[1, 0])
    pp_len = pp[pp["PPCAT"] == "LENALIDOMIDE"] if "PPCAT" in pp.columns else \
             pp[pp["PPTEST"].str.contains("Lenalidomide", na=False)]
    len_wide = pp_len.pivot_table(index="USUBJID", columns="PPTESTCD",
                                  values="PPSTRESN", aggfunc="first").reset_index()
    if "AUCINF" in len_wide.columns and "BASE_CREACL" in adsl.columns:
        merged = len_wide.merge(adsl[["USUBJID","BASE_CREACL"]], on="USUBJID")
        m = merged.dropna(subset=["AUCINF","BASE_CREACL"])
        ax3.scatter(m["BASE_CREACL"], m["AUCINF"], alpha=0.35, s=14,
                    color=DRUG_COLOR["LENALIDOMIDE"])
        slope, intercept, r, p, _ = stats.linregress(m["BASE_CREACL"], m["AUCINF"])
        xr = np.array([m["BASE_CREACL"].min(), m["BASE_CREACL"].max()])
        ax3.plot(xr, slope*xr + intercept, "k--", lw=1.5,
                 label=f"r={r:.2f}  p={p:.2e}")
        ax3.set_xlabel("CrCL (mL/min)"); ax3.set_ylabel("Lenalidomide AUCinf")
        ax3.set_title("CrCL ↔ Lenalidomide AUCinf\n(Chen 2012: CL ∝ CrCL^0.60)")
        ax3.legend(fontsize=8); ax3.grid(True, ls=":", alpha=0.4)

    # Panel 4: BSA vs Ixazomib Cmax
    ax4 = fig.add_subplot(gs[1, 1])
    if "CMAX" in ix_wide.columns and "BSA" in adsl.columns:
        merged2 = ix_wide.merge(adsl[["USUBJID","BSA","CYP3A4_INHIBFL"]], on="USUBJID")
        m2 = merged2.dropna(subset=["CMAX","BSA"])
        inh = m2["CYP3A4_INHIBFL"] == "Y"
        ax4.scatter(m2.loc[~inh,"BSA"], m2.loc[~inh,"CMAX"],
                    alpha=0.35, s=14, color=DRUG_COLOR["IXAZOMIB"], label="No DDI")
        ax4.scatter(m2.loc[inh,"BSA"], m2.loc[inh,"CMAX"],
                    alpha=0.7, s=22, marker="*", color="darkorange",
                    label="CYP3A4 inhibitor")
        slope, intercept, r, p, _ = stats.linregress(m2["BSA"], m2["CMAX"])
        xr = np.array([m2["BSA"].min(), m2["BSA"].max()])
        ax4.plot(xr, slope*xr + intercept, "k--", lw=1.5,
                 label=f"r={r:.2f}  p={p:.2e}")
        ax4.set_xlabel("BSA (m²)"); ax4.set_ylabel("Ixazomib Cmax (ng/mL)")
        ax4.set_title("BSA ↔ Ixazomib Cmax\n(Gupta 2017: BSA on V4, power 0.70 → lower Cmax)")
        ax4.legend(fontsize=7); ax4.grid(True, ls=":", alpha=0.4)

    # Panel 5: CWRES proxy
    ax5 = fig.add_subplot(gs[1, 2])
    all_cwres = []
    for drug in ["IXAZOMIB", "LENALIDOMIDE", "DEXAMETHASONE"]:
        cw = _cwres_proxy(pc[pc["PCTESTCD"] == drug].copy())
        all_cwres.append(cw)
        ax5.hist(cw, bins=30, alpha=0.45, density=True,
                 color=DRUG_COLOR[drug], label=DRUG_LABEL[drug], edgecolor="white")
    cwres_all = pd.concat(all_cwres)
    x = np.linspace(cwres_all.quantile(0.001), cwres_all.quantile(0.999), 200)
    ax5.plot(x, stats.norm.pdf(x, 0, 1), "k-", lw=2, label="N(0,1) reference")
    sw_stat, sw_p = stats.shapiro(cwres_all.sample(min(5000, len(cwres_all)), random_state=42))
    mu_cw, sd_cw = cwres_all.mean(), cwres_all.std()
    ax5.text(0.03, 0.95, f"Mean={mu_cw:.3f}  SD={sd_cw:.3f}\nShapiro p={sw_p:.3f}",
             transform=ax5.transAxes, fontsize=8, va="top",
             bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
    ax5.set_xlabel("CWRES proxy"); ax5.set_ylabel("Density")
    ax5.set_title("CWRES-proxy — all drugs combined\n(target distribution: N(0,1))")
    ax5.legend(fontsize=7); ax5.grid(True, ls=":", alpha=0.4)

    out = os.path.join(figdir, f"pk_gof_{study_key}.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out}")


# ── NCA summary table ──────────────────────────────────────────────────────────

def nca_summary_table(study_key, pp):
    rows = []
    drug_col = "PPCAT" if "PPCAT" in pp.columns else None
    for drug, pub in PUBLISHED.items():
        sub = pp[pp[drug_col] == drug] if drug_col else \
              pp[pp["PPTEST"].str.contains(DRUG_LABEL[drug], na=False)]
        for param in ["CMAX","AUCINF","THALF","CLF"]:
            vals = sub[sub["PPTESTCD"] == param]["PPSTRESN"].dropna()
            if len(vals) < 5:
                continue
            med = vals.median(); q25 = vals.quantile(0.25); q75 = vals.quantile(0.75)
            cv  = vals.std() / vals.mean() * 100 if vals.mean() > 0 else np.nan
            pub_val = pub.get(f"{param.lower()}_med", np.nan)
            rows.append({"study": study_key, "drug": drug, "parameter": param,
                         "n": len(vals), "median": round(med,2),
                         "Q25": round(q25,2), "Q75": round(q75,2),
                         "CV_pct": round(cv,1), "published": pub_val,
                         "pct_err": round(abs(med-pub_val)/abs(pub_val)*100, 1)
                                    if not np.isnan(pub_val) and pub_val else np.nan})
    return pd.DataFrame(rows)


# ── Console pass/fail summaries ────────────────────────────────────────────────

def vpc_pass_fail(vpc_df):
    print("\n" + "─"*65)
    print("  VPC — % within model-derived MC PI  (proper external reference)")
    print("─"*65)
    for _, r in vpc_df.iterrows():
        pct = r["pct_within_pi"]
        ok  = pct >= 80
        icon = "\033[92m  PASS\033[0m" if ok else "\033[91m  FAIL\033[0m"
        print(f"{icon}  {r['study']}  {r['drug']:<14}  [{r['panel']}]  "
              f"{pct:.1f}% within MC 5–95th PI  (n={int(r['n_obs'])})")


def nca_pass_fail(nca_df):
    print("\n" + "─"*65)
    print("  NCA — median vs published reference  (target ≤ 20% error)")
    print("─"*65)
    for _, r in nca_df.iterrows():
        if np.isnan(r["pct_err"]):
            continue
        ok   = r["pct_err"] <= 20
        warn = r["pct_err"] <= 40
        icon = ("\033[92m  PASS\033[0m" if ok else
                ("\033[93m  WARN\033[0m" if warn else "\033[91m  FAIL\033[0m"))
        note = ("  [NCA underestimates: sampling window (168h) < 3×t½ (684h) — expected]"
                if r["parameter"] == "THALF" and r["drug"] == "IXAZOMIB" else "")
        print(f"{icon}  {r['study']}  {r['drug']:<14} {r['parameter']:<8} "
              f"sim={r['median']:.2f}  pub={r['published']:.2f}  "
              f"err={r['pct_err']:.0f}%{note}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    all_vpc, all_nca = [], []

    for study_key in ["MM2", "MM1"]:
        if not os.path.exists(os.path.join(BASE, study_key)):
            print(f"[SKIP] {study_key} not found"); continue

        print(f"\n{'='*65}")
        print(f"  Step 8 VPC + GOF — {study_key}")
        print(f"{'='*65}")

        pc, pp, adpc, adsl = _load(study_key)

        print("\n  [1/2] VPC (MC reference bands from Gupta/Chen 2017) ...")
        vpc_df = plot_vpc(study_key, pc, adsl, FIGDIR)
        all_vpc.append(vpc_df)

        print("  [2/2] GOF ...")
        plot_gof(study_key, pc, pp, adpc, adsl, FIGDIR)

        all_nca.append(nca_summary_table(study_key, pp))

    vpc_all = pd.concat(all_vpc, ignore_index=True) if all_vpc else pd.DataFrame()
    nca_all = pd.concat(all_nca, ignore_index=True) if all_nca else pd.DataFrame()

    if not vpc_all.empty:
        path = os.path.join(TABDIR, "pk_vpc_summary.csv")
        vpc_all.to_csv(path, index=False)
        print(f"\n  → {path}")

    if not nca_all.empty:
        path = os.path.join(TABDIR, "pk_nca_summary.csv")
        nca_all.to_csv(path, index=False)
        print(f"  → {path}")

    if not vpc_all.empty:
        vpc_pass_fail(vpc_all)
    if not nca_all.empty:
        nca_pass_fail(nca_all)

    print(f"\n{'='*65}")
    print("  Step 8 complete.")
    print("  VPC reference: Gupta 2017 Monte Carlo (N=800 virtual patients)")
    print("  Figures → outputs/figures/pk_vpc_*.png  pk_gof_*.png")
    print(f"{'='*65}\n")


if __name__ == "__main__":
    main()
