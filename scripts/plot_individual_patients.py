"""
Individual Patient PK + PD Visualization
=========================================
Creates one figure per patient per trial (MM1, MM2) showing:
  - Ixazomib PK concentration (Cycles 1 & 3) — only for PK-sampled subjects
  - Lenalidomide PK concentration (Cycles 1 & 3) — only for PK-sampled subjects
  - Serum M-protein (SPEP_MPROT, g/L) and urine M-protein (UPEP_MPROT, mg/day)
  - Haemoglobin (HGB, g/L)
  - Neutrophils (NEUT, 10^9/L) and Platelets (PLT, 10^9/L) on dual y-axes

Dose times are shown as vertical lines coloured by drug and shaded by dose amount.
Time axis: study day (day 1 = first dose).

Output: outputs/figures/MM2/<USUBJID>.png
         outputs/figures/MM1/<USUBJID>.png
"""

import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # non-interactive — safe for batch generation
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = os.path.join(os.path.dirname(__file__), "..")   # Takeda-data/

# reference date: all EX/LB dates are relative to 2015-01-01 (study day 1)
REF_DATE = pd.Timestamp("2015-01-01")

# ── PK time helpers ──────────────────────────────────────────────────────────
# Cycle N starts at study day (N-1)*28 + 1
def cycle_start_day(visitnum: int) -> float:
    return (visitnum - 1) * 28 + 1

def pk_to_study_day(visitnum_series, pctptnum_series) -> pd.Series:
    """Convert PK visitnum + hours-after-dose to fractional study day."""
    return (visitnum_series - 1) * 28 + 1 + pctptnum_series / 24.0

def exstdtc_to_study_day(exstdtc_series) -> pd.Series:
    return (pd.to_datetime(exstdtc_series) - REF_DATE).dt.days + 1

# ── Colour palette ───────────────────────────────────────────────────────────
DRUG_COLORS = {
    "IXAZOMIB":     "#E05C3A",
    "LENALIDOMIDE": "#3A7EBB",
    "DEXAMETHASONE":"#5BAD6F",
    "PLACEBO":      "#AAAAAA",
}

NORMAL_RANGES = {        # (lo, hi) — drawn as a shaded band on PD panels
    "SPEP_MPROT": (0,   5),
    "UPEP_MPROT": (0,  80),
    "HGB":        (120, 160),
    "NEUT":       (1.8, 7.7),
    "PLT":        (150, 400),
}

# ── Data loading ─────────────────────────────────────────────────────────────

def load_study(study_key: str):
    sdir = os.path.join(BASE, study_key)

    pc   = pd.read_csv(os.path.join(sdir, "sdtm_pc.csv"))
    lb   = pd.read_csv(os.path.join(sdir, "sdtm_lb.csv"))
    ex   = pd.read_csv(os.path.join(sdir, "sdtm_ex.csv"))
    adsl = pd.read_csv(os.path.join(sdir, "adam_adsl.csv"))

    # ── normalise PK ─────────────────────────────────────────────────────────
    # DOSE_DAY is the absolute study day of each dose event (added in v2 PK gen).
    # STUDY_DAY = DOSE_DAY + hours-post-dose / 24 → precise fractional study day.
    # Fall back to cycle-based formula if DOSE_DAY column is absent (old files).
    pc["CONC"] = pd.to_numeric(pc["PCSTRESN"], errors="coerce")
    if "DOSE_DAY" in pc.columns:
        pc["STUDY_DAY"] = pd.to_numeric(pc["DOSE_DAY"], errors="coerce") \
                          + pc["PCTPTNUM"] / 24.0
    else:
        pc["STUDY_DAY"] = pk_to_study_day(pc["VISITNUM"], pc["PCTPTNUM"])

    # ── normalise LB ─────────────────────────────────────────────────────────
    lb["AVAL"] = pd.to_numeric(lb["LBSTRESN"], errors="coerce")

    # ── normalise EX ─────────────────────────────────────────────────────────
    ex["DOSE_DAY"] = exstdtc_to_study_day(ex["EXSTDTC"])
    ex["EXDOSE"]   = pd.to_numeric(ex["EXDOSE"], errors="coerce")

    # ── PK subject set ────────────────────────────────────────────────────────
    pk_subj = set(pc["USUBJID"].unique())

    return pc, lb, ex, adsl, pk_subj

# ── Per-subject data extraction ───────────────────────────────────────────────

def get_lb_series(lb_subj, testcd):
    d = lb_subj[lb_subj["LBTESTCD"] == testcd][["LBDY", "AVAL"]].dropna()
    return d["LBDY"].values, d["AVAL"].values

def get_pk_series(pc_subj, drug):
    d = pc_subj[(pc_subj["PCTESTCD"] == drug) & (pc_subj["BLQ"] == "N")][
        ["STUDY_DAY", "CONC", "VISITNUM"]].dropna()
    return d.sort_values("STUDY_DAY")

def get_ex_series(ex_subj, drug):
    d = ex_subj[ex_subj["EXTRT"] == drug][["DOSE_DAY", "EXDOSE"]].dropna()
    return d.sort_values("DOSE_DAY")

# ── Dose line helper ──────────────────────────────────────────────────────────

def draw_dose_lines(ax, dose_df, color, nominal_dose, alpha_max=0.35, lw=0.8):
    """
    Draw a vertical line for every dose event. Lines become more transparent
    when the dose is the nominal dose; reduced doses are shown at full alpha
    to highlight modifications.
    """
    if dose_df.empty:
        return
    for _, row in dose_df.iterrows():
        day  = row["DOSE_DAY"]
        dose = row["EXDOSE"]
        is_full = abs(dose - nominal_dose) < 0.05
        a = alpha_max * 0.5 if is_full else alpha_max
        lw_use = lw if is_full else lw * 1.6
        ax.axvline(day, color=color, alpha=a, linewidth=lw_use, zorder=1)

# ── Normal-range shading ──────────────────────────────────────────────────────

def shade_normal(ax, testcd, xlim, color="#90EE90", alpha=0.08):
    lo, hi = NORMAL_RANGES.get(testcd, (None, None))
    if lo is None:
        return
    ax.axhspan(lo, hi, color=color, alpha=alpha, zorder=0)

# ── Main figure builder ───────────────────────────────────────────────────────

NOMINAL_DOSES = {"IXAZOMIB": 4.0, "LENALIDOMIDE": 25.0,
                 "DEXAMETHASONE": 40.0, "PLACEBO": 0.0}

def make_figure(usubjid, pc_subj, lb_subj, ex_subj, adsl_row, has_pk: bool):
    """
    Build the per-patient figure. Returns the figure object.
    """
    arm    = adsl_row["ARMCD"] if "ARMCD" in adsl_row else "?"
    igtype = adsl_row["IGTYPE"] if "IGTYPE" in adsl_row else "?"
    age    = adsl_row["AGE"]   if "AGE"   in adsl_row else "?"
    sex    = adsl_row["SEX"]   if "SEX"   in adsl_row else "?"
    bm_mp  = adsl_row.get("BASE_SPEP_MPROT", np.nan)
    bm_hgb = adsl_row.get("BASE_HGB", np.nan)
    riss   = adsl_row.get("RISS", "?")

    # ── Layout ────────────────────────────────────────────────────────────────
    if has_pk:
        n_rows   = 5
        pk_rows  = [0, 1]
        pd_rows  = [2, 3, 4]
        h_ratios = [1.4, 1.4, 1.6, 1.2, 1.2]
    else:
        n_rows   = 3
        pk_rows  = []
        pd_rows  = [0, 1, 2]
        h_ratios = [1.6, 1.2, 1.2]

    fig, axes = plt.subplots(
        n_rows, 1, figsize=(14, 3.2 * n_rows),
        gridspec_kw={"height_ratios": h_ratios},
        sharex=False
    )
    if n_rows == 1:
        axes = [axes]

    fig.subplots_adjust(hspace=0.45, left=0.07, right=0.93, top=0.91, bottom=0.06)

    # ── Title ─────────────────────────────────────────────────────────────────
    bm_mp_s  = f"{bm_mp:.1f}" if not np.isnan(float(bm_mp))  else "N/A"
    bm_hgb_s = f"{bm_hgb:.0f}" if not np.isnan(float(bm_hgb)) else "N/A"
    fig.suptitle(
        f"{usubjid}   |   Arm: {arm}   |   Ig type: {igtype}   |   "
        f"Age: {age}  {sex}   |   R-ISS: {riss}\n"
        f"Baseline M-protein: {bm_mp_s} g/L   |   Baseline HGB: {bm_hgb_s} g/L   |   "
        f"Nominal doses — Ixazomib: 4 mg (weekly), Lenalidomide: 25 mg (daily d1-21), Dex: 40 mg (weekly)",
        fontsize=9, y=0.97, linespacing=1.4
    )

    # ── Dose event series ─────────────────────────────────────────────────────
    # Lenalidomide is given in BOTH arms (Days 1-21 per cycle) and is used as
    # the primary dosing marker on PD panels so Rd-arm patients show visible
    # dose rhythm. Ixazomib (or PLACEBO for Rd arm) is overlaid in addition on
    # IRd-arm panels. Placebo lines are suppressed on PD panels entirely.
    ex_ixaz  = get_ex_series(ex_subj, "IXAZOMIB" if arm == "IRd" else "PLACEBO")
    ex_lena  = get_ex_series(ex_subj, "LENALIDOMIDE")

    if has_pk:
        for i, drug in enumerate(["IXAZOMIB", "LENALIDOMIDE"]):
            ax = axes[pk_rows[i]]
            if drug == "IXAZOMIB" and arm != "IRd":
                ax.text(0.5, 0.5, "Rd arm — no Ixazomib (Placebo)",
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=10, color="grey")
                ax.set_ylabel("Ixazomib\n(ng/mL)", fontsize=8)
                ax.set_yticks([])
            else:
                pk_d = get_pk_series(pc_subj, drug)
                if pk_d.empty:
                    ax.text(0.5, 0.5, f"No {drug.title()} PK data",
                            ha="center", va="center", transform=ax.transAxes,
                            fontsize=10, color="grey")
                else:
                    color = DRUG_COLORS[drug]
                    # Plot per-cycle
                    for cyc, grp in pk_d.groupby("VISITNUM"):
                        label = f"Cycle {cyc}"
                        marker = "o" if cyc == 1 else "s"
                        ax.plot(grp["STUDY_DAY"], grp["CONC"],
                                color=color, marker=marker, markersize=4,
                                linewidth=1.2, label=label, zorder=4)
                    ax.set_yscale("log")
                    ax.set_ylabel(f"{drug.title()}\n(ng/mL, log)", fontsize=8)
                    ax.legend(fontsize=7, loc="upper right", framealpha=0.6)

            # dose lines on PK panel — use actual drug's EX for Ixazomib/Lena;
            # for the "IXAZOMIB" slot in Rd arm the placebo marker is suppressed
            # and we draw Lenalidomide lines instead so there's visible timing.
            if drug == "IXAZOMIB" and arm == "IRd":
                draw_dose_lines(ax, ex_ixaz, DRUG_COLORS["IXAZOMIB"],
                                NOMINAL_DOSES["IXAZOMIB"])
            else:
                # Lenalidomide panel (or Ixazomib slot for Rd arm): use Lena
                draw_dose_lines(ax, ex_lena, DRUG_COLORS["LENALIDOMIDE"],
                                NOMINAL_DOSES["LENALIDOMIDE"], alpha_max=0.25)

            ax.set_xlabel("")
            ax.tick_params(axis="both", labelsize=8)
            ax.spines[["top", "right"]].set_visible(False)
            # Cycle boundary annotations at top
            _annotate_cycles(ax, ex_subj)

    # ── PD: M-protein ─────────────────────────────────────────────────────────
    ax_mp = axes[pd_rows[0]]
    t_smp, v_smp = get_lb_series(lb_subj, "SPEP_MPROT")
    t_ump, v_ump = get_lb_series(lb_subj, "UPEP_MPROT")
    shade_normal(ax_mp, "SPEP_MPROT", None)

    if len(t_smp):
        ax_mp.plot(t_smp, v_smp, color="#B22222", marker="o", markersize=3.5,
                   linewidth=1.3, label="Serum M-protein (g/L)", zorder=4)

    ax_mp2 = ax_mp.twinx()
    if len(t_ump):
        ax_mp2.plot(t_ump, v_ump, color="#CC7722", marker="^", markersize=3.5,
                    linewidth=1.1, linestyle="--", label="Urine M-protein (mg/day)", zorder=4)
        ax_mp2.set_ylabel("Urine M-protein\n(mg/day)", fontsize=8, color="#CC7722")
        ax_mp2.tick_params(axis="y", labelsize=7, colors="#CC7722")
        ax_mp2.spines["right"].set_edgecolor("#CC7722")

    ax_mp.set_ylabel("Serum M-protein\n(g/L)", fontsize=8, color="#B22222")
    ax_mp.tick_params(axis="y", labelsize=8, colors="#B22222")
    ax_mp.set_xlabel("")
    ax_mp.spines[["top"]].set_visible(False)

    # combined legend
    h1, l1 = ax_mp.get_legend_handles_labels()
    h2, l2 = ax_mp2.get_legend_handles_labels()
    ax_mp.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper right", framealpha=0.6)
    # PD dose lines: Lenalidomide lines visible in BOTH arms (primary dose marker).
    # Ixazomib lines overlaid for IRd arm only (faint, so they don't dominate).
    _draw_pd_dose_lines(ax_mp, ex_ixaz, ex_lena, arm)
    _annotate_cycles(ax_mp, ex_subj)

    # ── PD: HGB ───────────────────────────────────────────────────────────────
    ax_hgb = axes[pd_rows[1]]
    t_hgb, v_hgb = get_lb_series(lb_subj, "HGB")
    shade_normal(ax_hgb, "HGB", None)
    if len(t_hgb):
        ax_hgb.plot(t_hgb, v_hgb, color="#4B0082", marker="o", markersize=3.5,
                    linewidth=1.3, label="HGB (g/L)", zorder=4)
    ax_hgb.axhline(80, color="#4B0082", linewidth=0.6, linestyle=":", alpha=0.5,
                   label="Grade 3 threshold (80 g/L)")
    ax_hgb.set_ylabel("HGB (g/L)", fontsize=8)
    ax_hgb.legend(fontsize=7, loc="upper right", framealpha=0.6)
    ax_hgb.spines[["top", "right"]].set_visible(False)
    ax_hgb.tick_params(axis="both", labelsize=8)
    _draw_pd_dose_lines(ax_hgb, ex_ixaz, ex_lena, arm)
    _annotate_cycles(ax_hgb, ex_subj)

    # ── PD: NEUT + PLT ────────────────────────────────────────────────────────
    ax_neut = axes[pd_rows[2]]
    t_neut, v_neut = get_lb_series(lb_subj, "NEUT")
    t_plt,  v_plt  = get_lb_series(lb_subj, "PLT")
    shade_normal(ax_neut, "NEUT", None, color="#ADD8E6")

    if len(t_neut):
        ax_neut.plot(t_neut, v_neut, color="#006994", marker="o", markersize=3.5,
                     linewidth=1.3, label="NEUT (10⁹/L)", zorder=4)
    ax_neut.axhline(0.5, color="#006994", linewidth=0.6, linestyle=":", alpha=0.5,
                    label="G4 neutropenia (<0.5)")

    ax_plt = ax_neut.twinx()
    if len(t_plt):
        ax_plt.plot(t_plt, v_plt, color="#8B0000", marker="s", markersize=3.5,
                    linewidth=1.1, linestyle="--", label="PLT (10⁹/L)", zorder=4)
    ax_plt.axhline(25, color="#8B0000", linewidth=0.6, linestyle=":", alpha=0.4)
    ax_plt.set_ylabel("Platelets (10⁹/L)", fontsize=8, color="#8B0000")
    ax_plt.tick_params(axis="y", labelsize=7, colors="#8B0000")
    ax_plt.spines["right"].set_edgecolor("#8B0000")

    ax_neut.set_ylabel("Neutrophils (10⁹/L)", fontsize=8, color="#006994")
    ax_neut.tick_params(axis="y", labelsize=8, colors="#006994")
    ax_neut.set_xlabel("Study Day", fontsize=9)
    ax_neut.spines["top"].set_visible(False)

    h1, l1 = ax_neut.get_legend_handles_labels()
    h2, l2 = ax_plt.get_legend_handles_labels()
    ax_neut.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper right", framealpha=0.6)
    _draw_pd_dose_lines(ax_neut, ex_ixaz, ex_lena, arm)
    _annotate_cycles(ax_neut, ex_subj)

    # ── Dose legend footer ────────────────────────────────────────────────────
    legend_elements = [
        Line2D([0], [0], color=DRUG_COLORS["LENALIDOMIDE"], linewidth=1.0, alpha=0.5,
               label="Lenalidomide dose event (both arms)"),
        mpatches.Patch(facecolor="#90EE90", alpha=0.3, label="Normal range"),
    ]
    if arm == "IRd":
        legend_elements.append(
            Line2D([0], [0], color=DRUG_COLORS["IXAZOMIB"], linewidth=1.0, alpha=0.5,
                   label="Ixazomib dose event (IRd arm only)")
        )
        legend_elements.append(
            Line2D([0], [0], color=DRUG_COLORS["IXAZOMIB"], linewidth=1.6, alpha=0.35,
                   label="Ixazomib dose reduction")
        )
    fig.legend(handles=legend_elements, loc="lower center", ncol=4,
               fontsize=7.5, bbox_to_anchor=(0.5, 0.01), framealpha=0.5)

    return fig


def _draw_pd_dose_lines(ax, ex_ixaz, ex_lena, arm):
    """
    Draw dose-timing markers on PD panels.

    Strategy (fixes invisible Rd-arm dose lines):
    - Lenalidomide lines are drawn for BOTH arms — they are the primary rhythm
      marker since all patients receive Lenalidomide Days 1-21 per cycle.
    - Ixazomib lines are additionally drawn for IRd arm (lighter, so they
      complement rather than dominate the Lenalidomide lines).
    - Placebo lines are NOT drawn — they carry no pharmacological signal and
      the old grey-at-0.175-alpha was essentially invisible.
    """
    # Lenalidomide: visible blue lines for both arms (Days 1, 8, 15 within cycle
    # plus daily through day 21 — downsampled to every-dose in EX records)
    draw_dose_lines(ax, ex_lena, DRUG_COLORS["LENALIDOMIDE"],
                    NOMINAL_DOSES["LENALIDOMIDE"], alpha_max=0.22)
    # Ixazomib: IRd arm only — orange lines at lower alpha so they show weekly
    # dosing days (1, 8, 15) without overplotting the Lena lines
    if arm == "IRd" and not ex_ixaz.empty:
        draw_dose_lines(ax, ex_ixaz, DRUG_COLORS["IXAZOMIB"],
                        NOMINAL_DOSES["IXAZOMIB"], alpha_max=0.30)


def _annotate_cycles(ax, ex_subj):
    """
    Mark cycle starts (every 28 days) with a light dashed grey line.
    Label every other cycle at the top of the axis.
    """
    # Determine max study day from EX
    if ex_subj.empty:
        return
    max_day = ex_subj["DOSE_DAY"].max()
    n_cycles = int(np.ceil(max_day / 28))

    for cyc in range(1, n_cycles + 2):
        day = (cyc - 1) * 28 + 1
        if day > max_day + 30:
            break
        ax.axvline(day, color="#CCCCCC", linewidth=0.6, linestyle="--", zorder=0, alpha=0.7)
        if cyc % 3 == 1:   # label every 3rd cycle to avoid clutter
            ylim = ax.get_ylim()
            y_top = ylim[1]
            ax.text(day, y_top, f"C{cyc}", fontsize=5.5, color="#999999",
                    ha="left", va="top", clip_on=True)


# ── Main generation loop ──────────────────────────────────────────────────────

def generate_study(study_key: str):
    print(f"\n{'='*65}")
    print(f"  Generating individual patient figures: {study_key}")
    print(f"{'='*65}")

    out_dir = os.path.join(BASE, "outputs", "figures", study_key)
    os.makedirs(out_dir, exist_ok=True)

    pc, lb, ex, adsl, pk_subj = load_study(study_key)

    subjects = adsl["USUBJID"].tolist()
    n = len(subjects)

    for i, usubjid in enumerate(subjects):
        if (i + 1) % 100 == 0 or i == 0:
            print(f"  [{i+1}/{n}] {usubjid}")

        pc_s  = pc[pc["USUBJID"]   == usubjid]
        lb_s  = lb[lb["USUBJID"]   == usubjid]
        ex_s  = ex[ex["USUBJID"]   == usubjid]
        row   = adsl[adsl["USUBJID"] == usubjid].iloc[0]

        has_pk = usubjid in pk_subj

        try:
            fig = make_figure(usubjid, pc_s, lb_s, ex_s, row, has_pk)
            out_path = os.path.join(out_dir, f"{usubjid}.png")
            fig.savefig(out_path, dpi=130, bbox_inches="tight",
                        facecolor="white", edgecolor="none")
        except Exception as e:
            print(f"    WARNING: {usubjid} failed — {e}")
        finally:
            plt.close("all")

    print(f"  Done. {n} figures saved to {out_dir}")


if __name__ == "__main__":
    studies = sys.argv[1:] or ["MM2", "MM1"]
    for sk in studies:
        generate_study(sk)
    print("\n✓ All figures generated.")
