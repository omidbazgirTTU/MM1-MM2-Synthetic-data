"""
TOURMALINE Synthetic Data Validation Script
============================================
Produces a full QC report against published TOURMALINE-MM1/MM2 targets.

Outputs
-------
outputs/validation/
    km_pfs.png               KM PFS curves (both studies, both arms)
    km_os.png                KM OS curves
    mprotein_waterfall.png   Best M-protein % change waterfall
    mprotein_spider.png      M-protein spider plot (sample of 40 patients/arm)
    plt_within_cycle.png     PLT nadir profile: Day 1 / 8 / 15 by arm
    pk_vpc_ixazomib.png      PK VPC — Ixazomib (Cycles 1 & 3)
    pk_vpc_lenalidomide.png  PK VPC — Lenalidomide
    covariate_dist.png       Key covariate distributions vs published targets
outputs/tables/
    validation_summary.csv   Pass/fail table for all key metrics
    response_rates.csv       ORR / VGPR / CR rates by study × arm
    covariate_table1.csv     Table 1 equivalent

Usage
-----
    python3 scripts/validate_data.py [MM2] [MM1]   (default: both)
"""

import os, sys, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

BASE     = os.path.join(os.path.dirname(__file__), "..")
OUT_FIG  = os.path.join(BASE, "outputs", "validation")
OUT_TAB  = os.path.join(BASE, "outputs", "tables")
os.makedirs(OUT_FIG, exist_ok=True)
os.makedirs(OUT_TAB, exist_ok=True)

# ── Published targets (Spec §2, §5, §6, §7) ──────────────────────────────────
TARGETS = {
    "MM2": {
        "pfs_ird":        35.3,    # months
        "pfs_rd":         21.8,
        "os_ird":         60.0,    # estimated (not mature)
        "os_rd":          48.0,
        "pfs_hr":         0.83,
        "n_ird":          351,
        "n_rd":           354,
        "median_age":     73,
        "pct_female":     45.0,
        "iss_i":          35.0,    # %
        "iss_ii":         35.0,
        "iss_iii":        30.0,
        "cyto_hr":        40.0,
        "del17p":         23.0,    # TOURMALINE-MM2 actual (calibrated for 40% high-risk)
        "t414":           18.0,    # TOURMALINE-MM2 actual (calibrated for 40% high-risk)
        "crcl_le60":      42.0,    # % patients CrCL ≤60
        "orr_ird":        82.0,    # Spec §5A MM2
        "orr_rd":         75.0,
        "vgpr_ird":       63.0,
        "vgpr_rd":        55.0,
        "cr_ird":         28.0,
        "cr_rd":          14.0,
        "g3_plt_ird":     25.0,    # % Grade 3 thrombocytopenia — Spec §5B
        "g3_plt_rd":      14.0,
    },
    "MM1": {
        "pfs_ird":        20.6,
        "pfs_rd":         14.7,
        "os_ird":         53.6,
        "os_rd":          51.6,
        "pfs_hr":         0.74,
        "n_ird":          360,
        "n_rd":           362,
        "median_age":     66,
        "pct_female":     43.0,
        "iss_i":          63.0,
        "iss_ii":         25.0,
        "iss_iii":        12.0,
        "cyto_hr":        20.0,
        "del17p":         10.0,
        "t414":           8.0,
        "crcl_le60":      30.0,
        "orr_ird":        78.0,    # Spec §5A MM1 / Moreau 2016
        "orr_rd":         72.0,
        "vgpr_ird":       48.0,
        "vgpr_rd":        39.0,
        "cr_ird":         12.0,
        "cr_rd":          7.0,
        "g3_plt_ird":     31.0,
        "g3_plt_rd":      16.0,
    },
}

ARM_COLORS = {"IRd": "#E05C3A", "Rd": "#3A7EBB"}
STUDY_COLORS = {"MM2": "#2E8B57", "MM1": "#8B2E5A"}

# ── KM estimator ──────────────────────────────────────────────────────────────

def km_curve(times_mo, events):
    """Returns (t_mo, surv) step-function arrays."""
    order  = np.argsort(times_mo)
    t      = np.array(times_mo)[order]
    e      = np.array(events)[order]
    n_risk = len(t)
    surv   = 1.0
    t_out  = [0.0];  s_out = [1.0]
    for i in range(len(t)):
        if e[i]:
            surv *= (1.0 - 1.0 / n_risk)
            t_out.append(t[i]);  s_out.append(surv)
        n_risk -= 1
    return np.array(t_out), np.array(s_out)

def km_median(times_mo, events):
    t, s = km_curve(times_mo, events)
    idx  = np.searchsorted(-s, -0.5)
    return t[idx] if idx < len(t) else float("inf")

def km_ci(times_mo, events, alpha=0.95):
    """Greenwood 95% CI at each event time."""
    order  = np.argsort(times_mo)
    t      = np.array(times_mo)[order]
    e      = np.array(events)[order]
    n      = len(t)
    n_risk = n; surv = 1.0; greenwood = 0.0
    t_out  = [0.0]; lo_out = [1.0]; hi_out = [1.0]
    for i in range(len(t)):
        if e[i]:
            d        = 1
            surv    *= (1.0 - d / n_risk)
            greenwood += d / (n_risk * (n_risk - d)) if n_risk > d else 0
            z  = 1.96
            se = surv * np.sqrt(greenwood)
            t_out.append(t[i])
            lo_out.append(max(0, surv - z * se))
            hi_out.append(min(1, surv + z * se))
        n_risk -= 1
    return np.array(t_out), np.array(lo_out), np.array(hi_out)

# ── Load data ─────────────────────────────────────────────────────────────────

def load_study(study_key):
    d = os.path.join(BASE, study_key)
    adsl  = pd.read_csv(f"{d}/adam_adsl.csv")
    adtte = pd.read_csv(f"{d}/adam_adtte.csv")
    adlb  = pd.read_csv(f"{d}/adam_adlb.csv")
    pc    = pd.read_csv(f"{d}/sdtm_pc.csv")
    return adsl, adtte, adlb, pc

# ── Figure 1: KM PFS + OS ─────────────────────────────────────────────────────

def plot_km(studies):
    fig, axes = plt.subplots(2, len(studies), figsize=(7 * len(studies), 10),
                             gridspec_kw={"hspace": 0.40, "wspace": 0.30})
    if len(studies) == 1:
        axes = axes.reshape(2, 1)

    results = {}   # store medians for validation table

    for col, sk in enumerate(studies):
        adsl, adtte, _, _ = load_study(sk)
        tgt = TARGETS[sk]

        for row, (param, ylabel, title) in enumerate([
            ("PFS", "Progression-Free Survival", "PFS"),
            ("OS",  "Overall Survival",           "OS"),
        ]):
            ax = axes[row][col]
            df = adtte[adtte["PARAMCD"] == param]

            med = {}
            for arm in ["IRd", "Rd"]:
                sub  = df[df["ARMCD"] == arm]
                t    = sub["AVAL"].values
                e    = (sub["CNSR"] == 0).astype(int).values
                t_km, s_km = km_curve(t, e)
                t_ci, lo, hi = km_ci(t, e)
                color = ARM_COLORS[arm]
                label = f"{arm}  (n={len(sub)})"
                ax.step(t_km, s_km, where="post", color=color, lw=2.0, label=label)
                ax.fill_between(t_ci, lo, hi, step="post", color=color, alpha=0.12)
                med[arm] = km_median(t, e)

            results[(sk, param)] = med

            # Published target lines
            for arm, key in [("IRd", f"{param.lower()}_ird"), ("Rd", f"{param.lower()}_rd")]:
                pub = tgt.get(key)
                if pub:
                    ax.axvline(pub, color=ARM_COLORS[arm], lw=1.0, ls="--", alpha=0.55)

            ax.set_xlim(0, df["AVAL"].max() * 1.05)
            ax.set_ylim(0, 1.05)
            ax.set_xlabel("Time (months)", fontsize=9)
            ax.set_ylabel("Survival Probability", fontsize=9)
            ax.set_title(f"{sk} — {title}", fontsize=10, fontweight="bold")
            ax.legend(fontsize=8, framealpha=0.7)
            ax.spines[["top", "right"]].set_visible(False)
            ax.tick_params(labelsize=8)

            # Median annotation
            for arm in ["IRd", "Rd"]:
                pub_key = f"{param.lower()}_{arm.lower()}"
                pub = tgt.get(pub_key)
                sim = med[arm]
                txt = (f"{arm}: sim={sim:.1f} mo\n"
                       f"      pub={pub:.1f} mo" if pub else f"{arm}: {sim:.1f} mo")
                ya  = 0.55 if arm == "IRd" else 0.40
                ax.text(0.97, ya, txt, transform=ax.transAxes,
                        fontsize=7, ha="right", va="top",
                        color=ARM_COLORS[arm], fontweight="bold")

    fig.suptitle("Kaplan–Meier Survival Curves\n(dashed lines = published targets)",
                 fontsize=11, y=0.98)
    path = os.path.join(OUT_FIG, "km_pfs_os.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {path}")
    return results

# ── Figure 2: M-protein waterfall + spider ────────────────────────────────────

def plot_mprotein(studies):
    fig, axes = plt.subplots(2, len(studies),
                             figsize=(8 * len(studies), 11),
                             gridspec_kw={"hspace": 0.50, "wspace": 0.35})
    if len(studies) == 1:
        axes = axes.reshape(2, 1)

    response_results = {}

    for col, sk in enumerate(studies):
        adsl, _, adlb, _ = load_study(sk)
        tgt = TARGETS[sk]

        mp = adlb[(adlb["PARAMCD"] == "SPEP_MPROT") & adlb["AVAL"].notna()].copy()
        bl = mp[mp["EPOCH"] == "BASELINE"].groupby("USUBJID")["AVAL"].first()
        best = mp[mp["EPOCH"] != "BASELINE"].groupby("USUBJID")["AVAL"].min()
        common = bl.index.intersection(best.index)
        pchg = ((best[common] - bl[common]) / bl[common] * 100).reset_index()
        pchg.columns = ["USUBJID", "PCHG"]
        pchg = pchg.merge(adsl[["USUBJID", "ARMCD"]], on="USUBJID")
        pchg = pchg.sort_values("PCHG")

        # ── Waterfall ──────────────────────────────────────────────────────
        ax_w = axes[0][col]
        colors = pchg["ARMCD"].map(ARM_COLORS).values
        ax_w.bar(range(len(pchg)), pchg["PCHG"].values, color=colors, width=1.0,
                 edgecolor="none", alpha=0.85)
        ax_w.axhline(-50,  color="black",  lw=0.8, ls="--", alpha=0.5, label="PR (−50%)")
        ax_w.axhline(-90,  color="#333333",lw=0.8, ls=":",  alpha=0.5, label="VGPR (−90%)")
        ax_w.axhline(  0,  color="black",  lw=0.6, alpha=0.4)
        ax_w.axhline( 25,  color="#CC0000",lw=0.8, ls="--", alpha=0.5, label="PD (+25%)")

        # Response rate from waterfall
        orr  = {}; vgpr = {}; cr = {}
        for arm in ["IRd", "Rd"]:
            sub = pchg[pchg["ARMCD"] == arm]["PCHG"]
            orr[arm]  = (sub <= -50).mean() * 100
            vgpr[arm] = (sub <= -90).mean() * 100
            cr[arm]   = (sub <= -99).mean() * 100

        response_results[sk] = {"orr": orr, "vgpr": vgpr, "cr": cr}

        legend_handles = [
            Line2D([0],[0], color=ARM_COLORS["IRd"], lw=4, label=f"IRd (n={len(pchg[pchg['ARMCD']=='IRd'])})"),
            Line2D([0],[0], color=ARM_COLORS["Rd"],  lw=4, label=f"Rd  (n={len(pchg[pchg['ARMCD']=='Rd'])})"),
            Line2D([0],[0], color="black", lw=1, ls="--", label="PR (−50%)"),
            Line2D([0],[0], color="#333333",lw=1, ls=":", label="VGPR (−90%)"),
            Line2D([0],[0], color="#CC0000",lw=1, ls="--", label="PD (+25%)"),
        ]
        ax_w.legend(handles=legend_handles, fontsize=7, loc="lower right", framealpha=0.7)
        ax_w.set_xlabel("Patients (ranked)", fontsize=9)
        ax_w.set_ylabel("Best M-protein % change", fontsize=9)
        ax_w.set_title(f"{sk} — M-protein Waterfall", fontsize=10, fontweight="bold")
        ax_w.tick_params(labelsize=8)
        ax_w.spines[["top", "right"]].set_visible(False)
        ax_w.set_xlim(-1, len(pchg))

        # ORR annotation
        for i, arm in enumerate(["IRd", "Rd"]):
            ax_w.text(0.02 + i*0.30, 0.97,
                      f"{arm}: ORR={orr[arm]:.0f}% (pub={tgt['orr_'+arm.lower()[:3].replace('rd','rd').replace('ird','ird')]:.0f}%)\n"
                      f"     VGPR+={vgpr[arm]:.0f}% (pub={tgt['vgpr_'+arm.lower()[:3].replace('rd','rd').replace('ird','ird')]:.0f}%)\n"
                      f"     CR  ={cr[arm]:.0f}% (pub={tgt['cr_'+arm.lower()[:3].replace('rd','rd').replace('ird','ird')]:.0f}%)",
                      transform=ax_w.transAxes, fontsize=7, va="top",
                      color=ARM_COLORS[arm], fontweight="bold",
                      bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.7))

        # ── Spider plot ────────────────────────────────────────────────────
        ax_s = axes[1][col]
        n_spider = 25   # patients per arm
        for arm in ["IRd", "Rd"]:
            subj_ids = adsl[adsl["ARMCD"] == arm]["USUBJID"].sample(
                min(n_spider, len(adsl[adsl["ARMCD"] == arm])),
                random_state=42
            ).values
            for uid in subj_ids:
                traj = mp[(mp["USUBJID"] == uid) & (mp["EPOCH"] != "BASELINE")][
                    ["VISITNUM", "AVAL"]].dropna()
                if uid not in bl.index or len(traj) < 2:
                    continue
                b = bl[uid]
                if b <= 0:
                    continue
                traj = traj.groupby("VISITNUM")["AVAL"].first().reset_index()
                pct  = (traj["AVAL"] - b) / b * 100
                ax_s.plot(traj["VISITNUM"].values, pct.values,
                          color=ARM_COLORS[arm], alpha=0.25, lw=0.8)

        ax_s.axhline(-50,  color="black",  lw=0.8, ls="--", alpha=0.5)
        ax_s.axhline(-90,  color="#333333",lw=0.8, ls=":",  alpha=0.5)
        ax_s.axhline( 25,  color="#CC0000",lw=0.8, ls="--", alpha=0.5)
        ax_s.axhline(  0,  color="black",  lw=0.5, alpha=0.4)
        ax_s.set_xlabel("Cycle", fontsize=9)
        ax_s.set_ylabel("M-protein % change from baseline", fontsize=9)
        ax_s.set_title(f"{sk} — M-protein Trajectories\n(n={n_spider}/arm, sample)",
                       fontsize=10, fontweight="bold")
        ax_s.tick_params(labelsize=8)
        ax_s.spines[["top", "right"]].set_visible(False)
        handles = [Line2D([0],[0], color=c, lw=2, label=a)
                   for a, c in ARM_COLORS.items()]
        ax_s.legend(handles=handles, fontsize=8, framealpha=0.7)

    path = os.path.join(OUT_FIG, "mprotein.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {path}")
    return response_results

# ── Figure 3: PLT within-cycle nadir profile ──────────────────────────────────

def plot_plt_profile(studies):
    fig, axes = plt.subplots(1, len(studies), figsize=(7 * len(studies), 5),
                             gridspec_kw={"wspace": 0.35})
    if len(studies) == 1:
        axes = [axes]

    plt_results = {}

    for ax, sk in zip(axes, studies):
        _, _, adlb, _ = load_study(sk)
        adsl, *_ = load_study(sk)

        plt_lb = adlb[(adlb["PARAMCD"] == "PLT") & adlb["AVAL"].notna()].copy()
        # adlb already contains ARMCD; only merge if missing
        if "ARMCD" not in plt_lb.columns:
            plt_lb = plt_lb.merge(adsl[["USUBJID", "ARMCD"]], on="USUBJID")

        # Spec §5B: Days 1, 8, 15 correspond to WEEKNUM 1, 2, 3
        wk_col = "WEEKNUM" if "WEEKNUM" in plt_lb.columns else None

        plt_results[sk] = {}
        tgt = TARGETS[sk]

        for arm in ["IRd", "Rd"]:
            sub = plt_lb[plt_lb["ARMCD"] == arm]
            days = [1, 8, 15]
            means = []; sems = []
            for wk in [1, 2, 3]:
                if wk_col:
                    v = sub[sub[wk_col] == wk]["AVAL"]
                else:
                    v = sub[sub["LBDY"] % 28 == (wk - 1) * 7 + 1]["AVAL"] \
                        if "LBDY" in sub.columns else pd.Series(dtype=float)
                if len(v) > 5:
                    means.append(v.mean())
                    sems.append(v.sem() * 1.96)
                else:
                    means.append(np.nan); sems.append(np.nan)

            color = ARM_COLORS[arm]
            means_arr = np.array(means); sems_arr = np.array(sems)
            valid = ~np.isnan(means_arr)
            ax.plot(np.array(days)[valid], means_arr[valid],
                    color=color, marker="o", ms=7, lw=2,
                    label=f"{arm} (n={len(sub['USUBJID'].unique())})")
            ax.fill_between(np.array(days)[valid],
                            (means_arr - sems_arr)[valid],
                            (means_arr + sems_arr)[valid],
                            color=color, alpha=0.15)

            plt_results[sk][arm] = means

        # Published reference lines
        ax.axhline(50,  color="red",    lw=0.8, ls="--", alpha=0.5, label="Grade 3 (<50)")
        ax.axhline(25,  color="darkred",lw=0.8, ls=":",  alpha=0.5, label="Grade 4 (<25)")
        ax.axhline(130, color="grey",   lw=0.8, ls="--", alpha=0.4, label="Spec nadir target (~130)")

        ax.set_xticks([1, 8, 15])
        ax.set_xticklabels(["Day 1\n(pre-dose)", "Day 8\n(post-dose)", "Day 15\n(nadir)"])
        ax.set_ylabel("Platelets (×10⁹/L)", fontsize=9)
        ax.set_title(f"{sk} — PLT Within-Cycle Profile\n(mean ± 95% CI across all cycles)",
                     fontsize=10, fontweight="bold")
        ax.legend(fontsize=8, framealpha=0.7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)

        # Grade 3 thrombocytopenia rate annotation
        for arm in ["IRd", "Rd"]:
            sub_nadir = plt_lb[(plt_lb["ARMCD"] == arm)]
            if wk_col:
                nadir_vals = sub_nadir[sub_nadir[wk_col] == 3]["AVAL"]
            else:
                nadir_vals = sub_nadir["AVAL"]
            if len(nadir_vals):
                g3 = (nadir_vals < 50).mean() * 100
                pub_g3 = tgt[f"g3_plt_{arm.lower()[:3].replace('ird','ird').replace('rd_','rd')}"]
                pub_key = "g3_plt_ird" if arm == "IRd" else "g3_plt_rd"
                pub_g3  = tgt[pub_key]
                ax.text(0.98, 0.95 if arm == "IRd" else 0.85,
                        f"{arm} G3 PLT: {g3:.0f}% (pub {pub_g3:.0f}%)",
                        transform=ax.transAxes, fontsize=8, ha="right",
                        color=ARM_COLORS[arm], fontweight="bold")

    path = os.path.join(OUT_FIG, "plt_within_cycle.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {path}")
    return plt_results


# ── Figure 4: PK VPC (Ixazomib + Lenalidomide) ────────────────────────────────

def plot_pk_vpc(studies):
    for drug, lloq, y_label in [
        ("IXAZOMIB",     0.5,  "Ixazomib concentration (ng/mL)"),
        ("LENALIDOMIDE", 2.0,  "Lenalidomide concentration (ng/mL)"),
    ]:
        fig, axes = plt.subplots(1, len(studies), figsize=(7 * len(studies), 5),
                                 gridspec_kw={"wspace": 0.35})
        if len(studies) == 1:
            axes = [axes]

        for ax, sk in zip(axes, studies):
            pc = pd.read_csv(os.path.join(BASE, sk, "sdtm_pc.csv"))
            sub = pc[(pc["PCTESTCD"] == drug) & (pc["BLQ"] == "N")].copy()
            sub["CONC"] = pd.to_numeric(sub["PCSTRESN"], errors="coerce")
            sub = sub.dropna(subset=["CONC", "PCTPTNUM"])

            if sub.empty:
                ax.text(0.5, 0.5, f"No {drug} PK data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=10)
                continue

            # Group by dose-day timepoint (PCTPTNUM) and compute percentiles
            grp = sub.groupby("PCTPTNUM")["CONC"]
            t_pts = sorted(grp.groups.keys())
            p05 = [grp.get_group(t).quantile(0.05) for t in t_pts]
            p50 = [grp.get_group(t).quantile(0.50) for t in t_pts]
            p95 = [grp.get_group(t).quantile(0.95) for t in t_pts]

            # Scatter observed
            ax.scatter(sub["PCTPTNUM"], sub["CONC"], color="#AAAAAA", s=4, alpha=0.3,
                       zorder=1, label="Observed")
            # Prediction intervals
            ax.plot(t_pts, p50, color="#E05C3A", lw=2.0, zorder=4, label="Median (sim)")
            ax.fill_between(t_pts, p05, p95, color="#E05C3A", alpha=0.15,
                            label="5th–95th pct (sim)")
            ax.axhline(lloq, color="black", lw=0.8, ls="--", alpha=0.5,
                       label=f"LLOQ = {lloq} ng/mL")

            # Cmax annotation
            cmax_med = sub.groupby("USUBJID")["CONC"].max().median()
            pub_cmax = 41 if drug == "IXAZOMIB" else 425
            ax.text(0.97, 0.95, f"Median Cmax: {cmax_med:.1f} ng/mL\nPublished: ~{pub_cmax} ng/mL",
                    transform=ax.transAxes, fontsize=8, ha="right", va="top",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

            ax.set_xlabel("Time post-dose (h)", fontsize=9)
            ax.set_ylabel(y_label, fontsize=9)
            ax.set_title(f"{sk} — {drug.title()} PK VPC", fontsize=10, fontweight="bold")
            ax.legend(fontsize=7, framealpha=0.7)
            ax.spines[["top", "right"]].set_visible(False)
            ax.tick_params(labelsize=8)

            cycles = sorted(sub["VISITNUM"].unique()) if "VISITNUM" in sub else []
            if cycles:
                ax.text(0.02, 0.98, f"Cycles: {cycles}", transform=ax.transAxes,
                        fontsize=7, va="top", color="#555555")

        path = os.path.join(OUT_FIG, f"pk_vpc_{drug.lower()}.png")
        fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"  Saved: {path}")


# ── Figure 5: Covariate distributions ─────────────────────────────────────────

def plot_covariates(studies):
    fig, axes = plt.subplots(3, 4, figsize=(18, 12),
                             gridspec_kw={"hspace": 0.55, "wspace": 0.40})
    axes = axes.flatten()

    all_data = {}
    for sk in studies:
        adsl, *_ = load_study(sk)
        all_data[sk] = adsl

    plots = [
        ("AGE",          "Age (years)",                None),
        ("WEIGHT",       "Weight (kg)",                None),
        ("BSA",          "BSA (m²)",                   None),
        ("BASE_CREACL",  "CrCL (mL/min)",              None),
        ("BASE_SPEP_MPROT","Baseline M-protein (g/L)", None),
        ("BASE_HGB",     "Baseline HGB (g/L)",         None),
        ("BASE_PLT",     "Baseline Platelets (×10⁹/L)",None),
        ("BASE_B2MG",    "Baseline β2M (mg/L)",        None),
    ]

    for i, (col, label, _) in enumerate(plots):
        ax = axes[i]
        for sk in studies:
            adsl = all_data[sk]
            if col not in adsl.columns:
                continue
            vals = pd.to_numeric(adsl[col], errors="coerce").dropna()
            ax.hist(vals, bins=30, alpha=0.55, color=STUDY_COLORS[sk],
                    label=sk, density=True)
            ax.axvline(vals.median(), color=STUDY_COLORS[sk], lw=1.5, ls="--",
                       alpha=0.8)
        ax.set_xlabel(label, fontsize=8)
        ax.set_ylabel("Density", fontsize=8)
        ax.legend(fontsize=7, framealpha=0.7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=7)

    # Categorical: ISS stage
    for i, (col, label, cats) in enumerate([
        ("ISSSTAGE",  "ISS Stage",           [1, 2, 3]),
        ("RISS",      "R-ISS Stage",         ["I", "II", "III"]),
        ("CYTOGR",    "Cytogenetic Risk",     ["HIGH RISK", "STANDARD RISK"]),
        ("ARMCD",     "Treatment Arm",        ["IRd", "Rd"]),
    ]):
        ax = axes[len(plots) + i]
        x = np.arange(len(cats))
        w = 0.35
        for j, sk in enumerate(studies):
            adsl = all_data[sk]
            if col not in adsl.columns:
                continue
            pcts = [(adsl[col] == c).mean() * 100 for c in cats]
            ax.bar(x + j * w, pcts, w, label=sk, color=STUDY_COLORS[sk], alpha=0.75)

            # Published targets overlay for ISS
            if col == "ISSSTAGE":
                tgt_pcts = [TARGETS[sk][f"iss_{['i','ii','iii'][k]}"]
                            for k in range(len(cats))]
                ax.plot(x + j * w, tgt_pcts, "k+", ms=8, markeredgewidth=1.5)

        ax.set_xticks(x + w / 2)
        ax.set_xticklabels([str(c) for c in cats], fontsize=8)
        ax.set_ylabel("%", fontsize=8)
        ax.set_title(label, fontsize=9, fontweight="bold")
        ax.legend(fontsize=7, framealpha=0.7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=7)

    fig.suptitle("Covariate Distributions — MM1 vs MM2\n(+ marks = published targets for ISS)",
                 fontsize=11, y=0.99)
    path = os.path.join(OUT_FIG, "covariate_distributions.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {path}")


# ── Tables ─────────────────────────────────────────────────────────────────────

def write_validation_table(studies, km_results, resp_results):
    rows = []

    def row(study, metric, simulated, target, unit=""):
        tol = abs(target) * 0.15 if target != 0 else 1.0
        status = "PASS" if abs(simulated - target) <= tol else "FAIL"
        rows.append({
            "Study": study, "Metric": metric,
            "Simulated": round(simulated, 2),
            "Target": round(target, 2),
            "Unit": unit,
            "Diff%": round((simulated - target) / target * 100, 1) if target != 0 else 0,
            "Status": status,
        })

    for sk in studies:
        adsl, adtte, adlb, pc = load_study(sk)
        tgt = TARGETS[sk]

        # PFS / OS medians
        for param in ["PFS", "OS"]:
            df = adtte[adtte["PARAMCD"] == param]
            for arm in ["IRd", "Rd"]:
                sub = df[df["ARMCD"] == arm]
                med = km_median(sub["AVAL"].values, (sub["CNSR"]==0).values)
                t_key = f"{param.lower()}_{arm.lower()}"
                row(sk, f"{param} median {arm}", med, tgt[t_key], "months")

        # Arm balance
        for arm in ["IRd", "Rd"]:
            n = (adsl["ARMCD"] == arm).sum()
            row(sk, f"N {arm}", n, tgt[f"n_{arm.lower()}"], "subjects")

        # Demographics
        row(sk, "Median age", adsl["AGE"].median(), tgt["median_age"], "years")
        pct_f = (adsl["SEX"] == "F").mean() * 100
        row(sk, "% female", pct_f, tgt["pct_female"], "%")

        # ISS
        for stage, key in [(1,"iss_i"),(2,"iss_ii"),(3,"iss_iii")]:
            pct = (adsl["ISSSTAGE"] == stage).mean() * 100
            row(sk, f"ISS Stage {stage}", pct, tgt[key], "%")

        # Cytogenetics
        row(sk, "High-risk cyto",
            (adsl["CYTOGR"] == "HIGH RISK").mean() * 100, tgt["cyto_hr"], "%")
        row(sk, "DEL17P",
            (adsl["DEL17P"] == "Y").mean() * 100, tgt["del17p"], "%")
        row(sk, "T(4;14)",
            (adsl["T414"] == "Y").mean() * 100, tgt["t414"], "%")

        # Renal function
        crcl = adsl["BASE_CREACL"] if "BASE_CREACL" in adsl.columns else None
        if crcl is not None:
            row(sk, "% CrCL ≤60", (crcl <= 60).mean() * 100, tgt["crcl_le60"], "%")

        # Response rates from resp_results
        if sk in resp_results:
            rr = resp_results[sk]
            for arm in ["IRd", "Rd"]:
                a = arm.lower()[:3]
                row(sk, f"ORR  {arm}", rr["orr"][arm],  tgt[f"orr_{a}"],  "%")
                row(sk, f"VGPR {arm}", rr["vgpr"][arm], tgt[f"vgpr_{a}"], "%")
                row(sk, f"CR   {arm}", rr["cr"][arm],   tgt[f"cr_{a}"],   "%")

        # PK Cmax (Ixazomib) — read NCA-derived Cmax from sdtm_pp.csv (Cycle 1),
        # not raw max from sdtm_pc.csv.  sdtm_pc contains multi-dose accumulated
        # concentrations; the published 41 ng/mL is a single-dose Cycle 1 NCA value.
        pp_path = os.path.join(BASE, sk, "sdtm_pp.csv")
        if os.path.exists(pp_path):
            pp = pd.read_csv(pp_path)
            ixaz_cmax = pp[
                (pp["PPTESTCD"] == "CMAX") &
                (pp["PPCAT"] == "IXAZOMIB") &
                (pp["VISITNUM"] == 1)
            ]["PPSTRESN"]
            if len(ixaz_cmax):
                cmax_median = pd.to_numeric(ixaz_cmax, errors="coerce").dropna().median()
                if not np.isnan(cmax_median):
                    row(sk, "Ixazomib Cmax (median)", cmax_median, 41.0, "ng/mL")

        # PLT Grade 3 — per-patient worst PLT across all cycles (matches clinical reporting)
        _, _, adlb, _ = load_study(sk)
        plt_lb = adlb[(adlb["PARAMCD"] == "PLT") & adlb["AVAL"].notna()].copy()
        if "ARMCD" not in plt_lb.columns:
            plt_lb = plt_lb.merge(adsl[["USUBJID","ARMCD"]], on="USUBJID")
        for arm in ["IRd", "Rd"]:
            sub = plt_lb[plt_lb["ARMCD"] == arm]
            if len(sub):
                worst_plt = sub.groupby("USUBJID")["AVAL"].min()
                g3 = (worst_plt < 50).mean() * 100
                pub_key = "g3_plt_ird" if arm == "IRd" else "g3_plt_rd"
                row(sk, f"Grade 3 PLT {arm}", g3, tgt[pub_key], "%")

    df = pd.DataFrame(rows)
    path = os.path.join(OUT_TAB, "validation_summary.csv")
    df.to_csv(path, index=False)
    print(f"  Saved: {path}")
    return df


def write_response_table(studies, resp_results):
    rows = []
    for sk in studies:
        tgt = TARGETS[sk]
        if sk not in resp_results:
            continue
        rr = resp_results[sk]
        for arm in ["IRd", "Rd"]:
            a = arm.lower()[:3]
            rows.append({
                "Study": sk, "Arm": arm,
                "ORR_sim":  round(rr["orr"][arm],  1),
                "ORR_pub":  tgt[f"orr_{a}"],
                "VGPR_sim": round(rr["vgpr"][arm], 1),
                "VGPR_pub": tgt[f"vgpr_{a}"],
                "CR_sim":   round(rr["cr"][arm],   1),
                "CR_pub":   tgt[f"cr_{a}"],
            })
    df = pd.DataFrame(rows)
    path = os.path.join(OUT_TAB, "response_rates.csv")
    df.to_csv(path, index=False)
    print(f"  Saved: {path}")
    return df


def write_table1(studies):
    rows = []
    for sk in studies:
        adsl, *_ = load_study(sk)
        tgt = TARGETS[sk]
        for arm in ["IRd", "Rd", "Overall"]:
            sub = adsl if arm == "Overall" else adsl[adsl["ARMCD"] == arm]
            n   = len(sub)
            def pct(mask): return f"{mask.sum()} ({mask.mean()*100:.1f}%)"
            def med_iqr(col):
                v = pd.to_numeric(sub[col], errors="coerce").dropna()
                return f"{v.median():.1f} ({v.quantile(0.25):.1f}–{v.quantile(0.75):.1f})"

            r = {"Study": sk, "Arm": arm, "N": n}
            r["Age median (IQR)"] = med_iqr("AGE")
            r["Female, n (%)"]    = pct(sub["SEX"] == "F")
            r["Weight kg (IQR)"]  = med_iqr("WEIGHT")
            r["CrCL mL/min (IQR)"]= med_iqr("BASE_CREACL") if "BASE_CREACL" in sub.columns else ""
            r["% CrCL ≤60"]       = f"{(pd.to_numeric(sub.get('BASE_CREACL',pd.Series()), errors='coerce')<=60).mean()*100:.1f}%"
            r["ISS I, n (%)"]     = pct(sub["ISSSTAGE"] == 1)
            r["ISS II, n (%)"]    = pct(sub["ISSSTAGE"] == 2)
            r["ISS III, n (%)"]   = pct(sub["ISSSTAGE"] == 3)
            r["R-ISS I, n (%)"]   = pct(sub["RISS"] == "I")   if "RISS" in sub else ""
            r["R-ISS II, n (%)"]  = pct(sub["RISS"] == "II")  if "RISS" in sub else ""
            r["R-ISS III, n (%)"] = pct(sub["RISS"] == "III") if "RISS" in sub else ""
            r["High-risk cyto, n (%)"] = pct(sub["CYTOGR"] == "HIGH RISK") if "CYTOGR" in sub else ""
            r["DEL17P, n (%)"]    = pct(sub["DEL17P"] == "Y") if "DEL17P" in sub else ""
            r["T(4;14), n (%)"]   = pct(sub["T414"] == "Y")   if "T414"   in sub else ""
            r["GAIN1Q21, n (%)"]  = pct(sub["GAIN1Q21"] == "Y") if "GAIN1Q21" in sub else ""
            r["IgG, n (%)"]       = pct(sub["IGTYPE"] == "IgG") if "IGTYPE" in sub else ""
            r["Baseline M-protein (IQR)"] = med_iqr("BASE_SPEP_MPROT") if "BASE_SPEP_MPROT" in sub.columns else ""
            r["Baseline HGB (IQR)"]       = med_iqr("BASE_HGB") if "BASE_HGB" in sub.columns else ""
            r["Baseline PLT (IQR)"]       = med_iqr("BASE_PLT") if "BASE_PLT" in sub.columns else ""
            if sk == "MM1":
                r["Prior lines 1, n (%)"] = pct(sub["NPRIORLINE"] == 1) if "NPRIORLINE" in sub else ""
                r["Prior PI, n (%)"]      = pct(sub["PRIORPHI"] == "Y") if "PRIORPHI"   in sub else ""
                r["Prior IMiD, n (%)"]    = pct(sub["PRIORIMID"] == "Y") if "PRIORIMID"  in sub else ""
            rows.append(r)

    df = pd.DataFrame(rows)
    path = os.path.join(OUT_TAB, "covariate_table1.csv")
    df.to_csv(path, index=False)
    print(f"  Saved: {path}")


# ── Summary print ──────────────────────────────────────────────────────────────

def print_summary(val_df):
    print(f"\n{'='*65}")
    print("  VALIDATION SUMMARY")
    print(f"{'='*65}")
    n_pass = (val_df["Status"] == "PASS").sum()
    n_fail = (val_df["Status"] == "FAIL").sum()
    print(f"  PASS: {n_pass}   FAIL: {n_fail}   Total: {len(val_df)}")
    print()

    fails = val_df[val_df["Status"] == "FAIL"]
    if len(fails):
        print("  ── Failing checks ──")
        for _, r in fails.iterrows():
            print(f"  {r['Study']:4s}  {r['Metric']:35s}  "
                  f"sim={r['Simulated']:.1f}  pub={r['Target']:.1f}  "
                  f"diff={r['Diff%']:+.0f}%")
    else:
        print("  All checks PASSED ✓")

    print()
    passes = val_df[val_df["Status"] == "PASS"]
    print("  ── Passing checks ──")
    for _, r in passes.iterrows():
        print(f"  {r['Study']:4s}  {r['Metric']:35s}  "
              f"sim={r['Simulated']:.1f}  pub={r['Target']:.1f}  "
              f"diff={r['Diff%']:+.0f}%  ✓")


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    studies = sys.argv[1:] or ["MM2", "MM1"]

    print(f"\nValidating: {studies}")
    print(f"Figures → {OUT_FIG}")
    print(f"Tables  → {OUT_TAB}")

    print("\n[1/6] KM Survival Curves...")
    km_results = plot_km(studies)

    print("\n[2/6] M-protein Waterfall + Spider...")
    resp_results = plot_mprotein(studies)

    print("\n[3/6] PLT Within-Cycle Profile...")
    plot_plt_profile(studies)

    print("\n[4/6] PK VPC...")
    plot_pk_vpc(studies)

    print("\n[5/6] Covariate Distributions...")
    plot_covariates(studies)

    print("\n[6/6] Writing validation tables...")
    val_df = write_validation_table(studies, km_results, resp_results)
    write_response_table(studies, resp_results)
    write_table1(studies)

    print_summary(val_df)
    print(f"\n✓  Validation complete.")
