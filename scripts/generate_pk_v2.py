"""
TOURMALINE PK Data Generator v2
Ixazomib: 3-compartment model per Gupta et al. Clin Pharmacokinet 2017
"""
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import os

RNG = np.random.default_rng(1234)
OUT_BASE = os.path.join(os.path.dirname(__file__), "..")  # Takeda-data/

PK_DRUGS = {
    "IXAZOMIB": {
        "model": "3cmt_oral",
        "dose_mg": 4.0, "route": "ORAL", "conc_unit": "ng/mL", "time_unit": "h",
        # ─── 3-cmt params (Gupta 2017 TOURMALINE-MM1 labeling paper) ───
        # CL=1.86 L/h, Vss=543 L, t½=9.5 d=228 h, F=58%, Ka=0.5 h⁻¹
        # V2=14.3 L (central, Gupta 2015 Table 3)
        # V3+V4 = 543-14.3 = 528.7 L split as V3=200, V4=328.7
        # Q3=8 L/h (fast distribution), Q4=0.5 L/h (slow/deep)
        "Ka": 0.50, "ALAG1": 0.15,
        "CL": 1.86, "V2": 14.3,
        "Q3": 8.00, "V3": 200.0,
        "Q4": 0.50, "V4": 328.7,
        "F": 0.58,
        "theta_BSA_V4": 0.70,
        "iiv_CL": 0.36, "iiv_V2": 0.30, "iiv_Q3": 0.45,
        "iiv_V3": 0.50, "iiv_Q4": 0.40, "iiv_V4": 0.45,
        "iiv_Ka": 0.55, "iiv_F": 0.30,
        "sigma_prop": 0.18, "sigma_add": 0.3,
        "nominal_times": [0.0,0.5,1.0,2.0,4.0,8.0,24.0,48.0,72.0,168.0,336.0],
        "sparse_times":  [0.0,1.0,4.0,24.0,168.0,336.0],  # 336h = 14d; needed to capture terminal t½=228h
        "cycles_assessed": [1,3], "lloq": 0.5,
    },
    "LENALIDOMIDE": {
        "model": "1cmt_oral",
        "dose_mg": 25.0, "route": "ORAL", "conc_unit": "ng/mL", "time_unit": "h",
        "Ka": 1.10, "CL_F": 11.36, "Vd_F": 80.0, "F": 1.0,  # CL/F=dose/AUCinf=25000/2200; Chen 2012
        "cl_renal_slope": 0.007,
        "iiv_Ka": 0.50, "iiv_CL_F": 0.32, "iiv_Vd_F": 0.38,
        "sigma_prop": 0.18, "sigma_add": 8.0,
        "nominal_times": [0.0,0.5,1.0,1.5,2.0,3.0,4.0,6.0,8.0,12.0,24.0],
        "sparse_times":  [0.0,1.0,2.0,6.0,24.0],
        "cycles_assessed": [1,3], "lloq": 2.0,
    },
    "DEXAMETHASONE": {
        "model": "1cmt_oral",
        "dose_mg": 40.0, "route": "ORAL", "conc_unit": "ng/mL", "time_unit": "h",
        "Ka": 1.50, "CL_F": 22.2, "Vd_F": 208.0, "F": 1.0,  # CL/F=dose/AUCinf=40000/1800; Vd/F=CL/F*t½/ln2=22.2*6.5/0.693
        "iiv_Ka": 0.45, "iiv_CL_F": 0.28, "iiv_Vd_F": 0.32,
        "sigma_prop": 0.20, "sigma_add": 3.0,
        "nominal_times": [0.0,0.5,1.0,2.0,3.0,4.0,6.0,8.0,12.0,24.0],
        "sparse_times":  [0.0,1.0,3.0,8.0,24.0],
        "cycles_assessed": [1,3], "lloq": 0.2,
    },
}

PP_META = {
    "CMAX":          ("CMAX",  "Maximum Observed Concentration",            "ng/mL"),
    "TMAX":          ("TMAX",  "Time of Maximum Observed Concentration",    "h"),
    "AUC_LAST":      ("AUCLST","AUC from Time 0 to Last Measurable Conc",  "ng*h/mL"),
    "AUC_INF":       ("AUCINF","AUC from Time 0 Extrapolated to Infinity", "ng*h/mL"),
    "T_HALF":        ("THALF", "Apparent Terminal Half-Life",               "h"),
    "CL_F":          ("CLF",   "Apparent Oral Clearance",                   "L/h"),
    "VZ_F":          ("VZF",   "Apparent Volume of Distribution",           "L"),
    "AUCPCT_EXTRAP": ("AUCPEX","Percent AUC Extrapolated",                  "%"),
}

# ── PK Models ────────────────────────────────────────────────────────────────

def pk_3cmt_oral(t_eval, dose, Ka, ALAG, CL, V2, Q3, V3, Q4, V4, F=0.58):
    ke  = CL / V2
    k23 = Q3 / V2; k32 = Q3 / V3
    k24 = Q4 / V2; k42 = Q4 / V4
    def odes(t, y):
        depot, central, p1, p2 = y
        ab = Ka * depot if t >= ALAG else 0.0
        dd = -Ka * depot if t >= ALAG else 0.0
        dc = F * ab - (ke + k23 + k24) * central + k32 * p1 + k42 * p2
        dp1 = k23 * central - k32 * p1
        dp2 = k24 * central - k42 * p2
        return [dd, dc, dp1, dp2]
    sol = solve_ivp(odes, [0.0, t_eval[-1]+1e-6], [dose,0,0,0],
                    t_eval=t_eval, method='RK45', rtol=1e-6, atol=1e-9)
    return np.maximum(sol.y[1] / V2 * 1000.0, 0.0)

def pk_1cmt_oral(t_eval, dose, Ka, CL_F, Vd_F):
    ke = CL_F / Vd_F
    if abs(Ka - ke) < 1e-6: ke *= 1.001
    C = (dose * Ka) / (Vd_F * (Ka - ke)) * (np.exp(-ke*t_eval) - np.exp(-Ka*t_eval))
    return np.maximum(C * 1000.0, 0.0)

def sim_profile(drug, row, times):
    p = PK_DRUGS[drug]
    if p["model"] == "3cmt_oral":
        return pk_3cmt_oral(times, p["dose_mg"],
                            row["Ka_i"], row["ALAG"],
                            row["CL_i"], row["V2_i"],
                            row["Q3_i"], row["V3_i"],
                            row["Q4_i"], row["V4_i"], row["F_i"])
    return pk_1cmt_oral(times, p["dose_mg"], row["Ka_i"], row["CL_F_i"], row["Vd_F_i"])

def add_res(c_true, drug):
    p = PK_DRUGS[drug]
    c_obs = c_true * (1 + RNG.normal(0, p["sigma_prop"], len(c_true))) \
            + RNG.normal(0, p["sigma_add"], len(c_true))
    c_obs = np.maximum(c_obs, 0.0)
    blq = c_obs < p["lloq"]
    c_obs[blq] = p["lloq"] / 2
    return c_obs, blq

def eta(cv, n):
    return np.exp(RNG.normal(0, cv, n))

# ── Multi-dose superposition ──────────────────────────────────────────────────
# Reference date: all EX/LB dates are relative to 2015-01-01 (study day 1)
_REF_DATE = pd.Timestamp("2015-01-01")

# Post-dose sampling timepoints (hours) — 4 measurements per dosing interval
# Per Spec §8A: pre-dose trough / Cmax (~1h) / mid-interval / late
#
# Ixazomib (Days 1,8,15 weekly):  [0=trough, 1=Cmax, 4=mid, 8=late]
#   Day 1 Cycle 1: no pre-dose (first dose) → 0h is "time 0" (pre-dosing reference)
# Lenalidomide (daily D1-21):     [0=trough, 2=Tmax, 4=mid, 8=late]
#   (Tmax ~2h for lenalidomide; Chen 2012)
# Dexamethasone (Days 1,8,15,22): [0=trough, 1.5=Cmax, 4=mid, 8=late]
#   (Tmax ~1-1.5h for dexamethasone)
SAMPLING_POST_DOSE_H = {
    "IXAZOMIB":      [0.0, 1.0,  4.0,  8.0],   # Spec §8A — Cmax at 1h
    "LENALIDOMIDE":  [0.0, 2.0,  4.0,  8.0],   # Tmax ~2h
    "DEXAMETHASONE": [0.0, 1.5,  4.0,  8.0],   # Tmax ~1-1.5h
}
# For Lenalidomide (21 daily doses/cycle) sample only 4 representative days
# (0-indexed offset from cycle-start day): day 1, 8, 15, 21
_LENA_SAMPLE_OFFSETS = {0, 7, 14, 20}


def _unit_dose_grid(drug: str, row, t_max: float) -> tuple:
    """
    Precompute concentration profile for ONE nominal dose on a fine time grid.
    Returns (t_grid, c_grid).  Actual dose contributions scale linearly:
        C_actual(t) = C_unit(t) * (actual_dose / nominal_dose)
    t_max covers ~22 × Ixazomib t½ = 5 000 h to capture >99 % elimination.
    """
    t_grid = np.linspace(0.0, t_max, 10_000)
    c_grid = sim_profile(drug, row, t_grid)
    return t_grid, c_grid


def _superpose(t_grid, c_grid, dose_times_h, dose_amounts_mg,
               nominal_mg, eval_times_h) -> np.ndarray:
    """
    Linear superposition of scaled single-dose profiles.
        C(t) = Σ_i  (d_i / nominal)  ×  C_unit(t − t_i)
    dose_times_h   : absolute hours of every dose event (from study day 1)
    dose_amounts_mg: administered mg for each dose event
    nominal_mg     : dose used when computing c_grid
    eval_times_h   : absolute hours at which to evaluate total concentration
    """
    from scipy.interpolate import interp1d
    c_unit = interp1d(t_grid, c_grid, kind="linear",
                      bounds_error=False, fill_value=0.0)
    c_total = np.zeros(len(eval_times_h))
    for t_d, d_mg in zip(dose_times_h, dose_amounts_mg):
        scale = float(d_mg) / float(nominal_mg) if nominal_mg > 0 else 0.0
        c_total += scale * c_unit(eval_times_h - t_d)
    return np.maximum(c_total, 0.0)

def sample_ixaz(adsl):
    """
    Sample individual Ixazomib PK parameters.
    Covariates applied:
      - BSA on V4: V4 = TV_V4 * (BSA/1.73)^0.70  [Gupta 2017]
      - CYP3A4 DDI on CL: CL *= CL_DDI_MULT from ADSL

    Track A2 — OMEGA Cholesky for correlated PK etas (CL, V2, V4):
      Off-diagonals from Cross_Correlations_Synthetic_Data_Guide.md §2B (assumed):
        Cov(CL,V2)=0.042 (r=0.40), Cov(CL,V4)=0.046 (r=0.29), Cov(V2,V4)=0.027 (r=0.20)

    Track A1/A3 — If IXAZ_CL_I exists in ADSL (set by generate_v2.py make_adsl()),
      use that pre-computed CL_i directly instead of re-sampling, ensuring both
      generators use the EXACT same per-patient clearance value.
    """
    p = PK_DRUGS["IXAZOMIB"]; n = len(adsl)

    # BSA covariate on deep peripheral volume V4
    bsa_ref = 1.73  # reference BSA (m²)
    bsa_i   = adsl["BSA"].values if "BSA" in adsl.columns else np.full(n, bsa_ref)
    bsa_i   = np.where(np.isnan(bsa_i.astype(float)), bsa_ref, bsa_i.astype(float))
    v4_bsa  = p["V4"] * (bsa_i / bsa_ref) ** p["theta_BSA_V4"]

    # CYP3A4 DDI multiplier on CL
    ddi_mult = adsl["CL_DDI_MULT"].values if "CL_DDI_MULT" in adsl.columns else np.ones(n)
    ddi_mult = np.where(np.isnan(ddi_mult.astype(float)), 1.0, ddi_mult.astype(float))

    # ── Track A2: OMEGA Cholesky for correlated CL, V2, V4 etas ─────────────
    # OMEGA matrix (variances + covariances):
    #   ω_CL=0.35, ω_V2=0.30, ω_V4=0.45 [Gupta 2017 IIV estimates]
    #   Off-diagonal covariances assumed (correlations not published for Ixazomib)
    _OMEGA = np.array([
        [0.1225, 0.042,  0.046],   # CL: ω²=0.35², Cov(CL,V2), Cov(CL,V4)
        [0.042,  0.09,   0.027],   # V2: ω²=0.30²
        [0.046,  0.027,  0.2025],  # V4: ω²=0.45²
    ])
    _L_omega = np.linalg.cholesky(_OMEGA)

    if "IXAZ_CL_I" in adsl.columns:
        # ── Track A1/A3: read pre-computed CL_i from ADSL ─────────────────────
        # CL_i already includes DDI correction (computed in make_adsl()).
        # For V2 and V4 etas: draw from conditional distribution given eta_CL.
        # Simplified: recover eta_CL from stored CL_i, then use OMEGA Cholesky
        # to get correlated eta_V2, eta_V4 conditional on eta_CL.
        cl_base   = adsl["IXAZ_CL_I"].values.astype(float)
        cl_base   = np.where(np.isnan(cl_base) | (cl_base <= 0), p["CL"], cl_base)
        eta_cl    = np.log(cl_base / (p["CL"] * ddi_mult))  # recover eta_CL

        # Conditional draw for eta_V2, eta_V4 given eta_CL:
        # Cholesky: eta_raw = L^{-1} @ eta → eta_raw[0] = eta_CL / L[0,0]
        # Remaining rows: eta_raw[1,2] ~ N(0,1) independent → correlated etas
        eta_raw_0 = eta_cl / _L_omega[0, 0]           # standardized CL component
        eta_raw_12 = RNG.standard_normal((n, 2))       # independent V2, V4 components
        eta_v2 = _L_omega[1, 0] * eta_raw_0 + _L_omega[1, 1] * eta_raw_12[:, 0]
        eta_v4 = _L_omega[2, 0] * eta_raw_0 + _L_omega[2, 1] * eta_raw_12[:, 0] \
               + _L_omega[2, 2] * eta_raw_12[:, 1]
    else:
        # Fallback: draw all 3 etas jointly from OMEGA Cholesky
        eta_raw = RNG.standard_normal((n, 3))
        etas    = ((_L_omega @ eta_raw.T).T)   # shape (n, 3): [η_CL, η_V2, η_V4]
        cl_base = p["CL"] * np.exp(etas[:, 0]) * ddi_mult
        eta_v2  = etas[:, 1]
        eta_v4  = etas[:, 2]

    return pd.DataFrame({
        "USUBJID": adsl["USUBJID"].values,
        "CL_i":  cl_base,
        "V2_i":  p["V2"] * np.exp(eta_v2),
        "Q3_i":  p["Q3"] * eta(p["iiv_Q3"], n),
        "V3_i":  p["V3"] * eta(p["iiv_V3"], n),
        "Q4_i":  p["Q4"] * eta(p["iiv_Q4"], n),
        "V4_i":  v4_bsa  * np.exp(eta_v4),
        "Ka_i":  p["Ka"] * eta(p["iiv_Ka"], n),
        "F_i":   np.clip(p["F"] * eta(p["iiv_F"], n), 0.1, 1.0),
        "ALAG":  p["ALAG1"],
    })

def sample_1cmt(drug, adsl):
    """
    Sample individual 1-compartment PK parameters.
    Lenalidomide: CrCL (Cockcroft-Gault) covariate on CL/F.
    CrCL covariate: CL = TV_CL * (CrCL/80)^0.60  (Chen et al. 2012)
    """
    p = PK_DRUGS[drug]; n = len(adsl)
    cl = p["CL_F"] * eta(p["iiv_CL_F"], n)
    if drug == "LENALIDOMIDE" and "BASE_CREACL" in adsl.columns:
        crcl = adsl["BASE_CREACL"].values.astype(float)
        crcl = np.where(np.isnan(crcl), 60.0, crcl).clip(10.0, 150.0)
        # Power-function covariate: (CrCL/80)^0.60  [Chen 2012 J Clin Pharmacol]
        cl *= (crcl / 80.0) ** 0.60
    return pd.DataFrame({
        "USUBJID": adsl["USUBJID"].values,
        "Ka_i":   p["Ka"]   * eta(p["iiv_Ka"],    n),
        "CL_F_i": cl,
        "Vd_F_i": p["Vd_F"] * eta(p["iiv_Vd_F"], n),
    })

# ── NCA ──────────────────────────────────────────────────────────────────────

def nca(times, concs, dose, lloq):
    mask = concs >= lloq
    t, c = times[mask], concs[mask]
    nan_d = {k: np.nan for k in ["CMAX","TMAX","AUC_LAST","AUC_INF","T_HALF","CL_F","VZ_F","AUCPCT_EXTRAP"]}
    if len(t) < 3: return nan_d
    cmax = c.max(); tmax = t[c.argmax()]
    auc_last = sum(
        (t[i]-t[i-1])*(c[i-1]-c[i])/np.log(c[i-1]/c[i])
        if c[i-1]>0 and c[i]>0 and c[i]<c[i-1]
        else (t[i]-t[i-1])*(c[i-1]+c[i])/2
        for i in range(1, len(t))
    )
    tidx = np.where(t > tmax)[0]
    lam_z = np.nan
    if len(tidx) >= 3:
        tt = t[tidx[-min(5,len(tidx)):]]
        cc = np.log(np.maximum(c[tidx[-min(5,len(tidx)):]], 1e-12))
        if len(tt) >= 2:
            s, _ = np.polyfit(tt, cc, 1)
            if s < 0: lam_z = -s
    t12 = 0.693/lam_z if not np.isnan(lam_z) else np.nan
    auc_ex = c[-1]/lam_z if not np.isnan(lam_z) else np.nan
    auc_inf = auc_last + auc_ex if not np.isnan(auc_ex) else np.nan
    pct = 100*auc_ex/auc_inf if not np.isnan(auc_inf) else np.nan
    clf = dose*1000/auc_inf if not np.isnan(auc_inf) else np.nan
    vzf = dose*1000/(auc_inf*lam_z) if (not np.isnan(auc_inf) and not np.isnan(lam_z)) else np.nan
    return {
        "CMAX": round(float(cmax),3), "TMAX": round(float(tmax),2),
        "AUC_LAST": round(float(auc_last),2),
        "AUC_INF":  round(float(auc_inf),2) if not np.isnan(auc_inf) else np.nan,
        "T_HALF":   round(float(t12),2)      if not np.isnan(t12)    else np.nan,
        "CL_F":     round(float(clf),4)      if not np.isnan(clf)    else np.nan,
        "VZ_F":     round(float(vzf),2)      if not np.isnan(vzf)    else np.nan,
        "AUCPCT_EXTRAP": round(float(pct),1) if not np.isnan(pct)    else np.nan,
    }

# ── Builders ─────────────────────────────────────────────────────────────────

def get_pk_subj(adsl):
    ird = adsl[adsl["ARMCD"]=="IRd"].sample(min(75,(adsl["ARMCD"]=="IRd").sum()), random_state=42)
    rd  = adsl[adsl["ARMCD"]=="Rd" ].sample(min(75,(adsl["ARMCD"]=="Rd" ).sum()), random_state=42)
    return pd.concat([ird, rd]).reset_index(drop=True)


def make_pc(adsl, ex, study_key):
    """
    Generate sdtm_pc aligned with the TOURMALINE Schedule of Assessments.

    Per SoA: PK sampling is a sparse-sampling substudy in a SUBSET of patients.
    "Timing: Selected cycles (typically Cycle 1, Day 1). Samples: Pre-dose and
    post-dose time points. Purpose: Population PK modelling."

    Implemented PK substudy design (3 tiers):
    ─────────────────────────────────────────
    Tier 1 — DENSE (Cycles 1 & 3): all weekly dose events × ≥4 timepoints
             Ixazomib:     [0, 4, 24, 72 h] after doses on days 1, 8, 15
                           → 12 obs / cycle (satisfies ≥4 per 7-day interval)
             Lenalidomide: [0, 2, 8, 24 h] on representative days 1, 8, 15, 21
                           → 16 obs / cycle (satisfies ≥4 per 24-h interval)
             Dexamethasone:[0, 4, 24, 72 h] after doses on days 1, 8, 15, 22
                           → 16 obs / cycle

    Tier 2 — TROUGH (Cycles 5, 7, 9, … every 2nd cycle after cycle 3):
             Pre-dose (0 h) trough only at Day 1 of each sampled cycle.
             Purpose: characterise inter-cycle accumulation in population PK.

    Tier 3 — NOT SAMPLED: all other cycles (most patients in real trial)

    Accumulation is captured via linear superposition of scaled single-dose
    profiles (valid for all linear compartment models).
    """
    sid      = adsl["STUDYID"].iloc[0]
    pk_subj  = get_pk_subj(adsl)

    # ── Pre-process EX: convert calendar dates → absolute study hours ─────────
    ex = ex.copy()
    ex["DOSE_DAY"] = (pd.to_datetime(ex["EXSTDTC"]) - _REF_DATE).dt.days + 1  # 1-indexed
    ex["DOSE_H"]   = (ex["DOSE_DAY"] - 1) * 24.0   # absolute hours from study start
    ex["EXDOSE"]   = pd.to_numeric(ex["EXDOSE"], errors="coerce").fillna(0.0)
    ex["DAY_IN_CYCLE_0IDX"] = (ex["DOSE_DAY"] - 1) % 28   # 0-indexed day within cycle

    rows = []
    for drug in ["IXAZOMIB", "LENALIDOMIDE", "DEXAMETHASONE"]:
        p         = PK_DRUGS[drug]
        nom_dose  = p["dose_mg"]
        post_h    = SAMPLING_POST_DOSE_H[drug]

        # Ixazomib only in IRd arm; others in both arms
        subj_df = pk_subj[pk_subj["ARMCD"] == "IRd"].copy() if drug == "IXAZOMIB" \
                  else pk_subj.copy()
        ip      = sample_ixaz(subj_df) if drug == "IXAZOMIB" else sample_1cmt(drug, subj_df)
        subj_df = subj_df.merge(ip, on="USUBJID")

        for _, row in subj_df.iterrows():
            uid = row["USUBJID"]
            arm = row["ARMCD"]
            # In Rd arm, Ixazomib dose events are labelled PLACEBO in EX
            extrt = drug if (arm == "IRd" or drug != "IXAZOMIB") else "PLACEBO"

            ex_s = ex[(ex["USUBJID"] == uid) & (ex["EXTRT"] == extrt)].copy()
            if ex_s.empty:
                continue

            # All dose times and amounts (used for superposition)
            all_dose_h  = ex_s["DOSE_H"].values
            all_dose_mg = ex_s["EXDOSE"].values

            # Subset of dose events to *sample* at
            if drug == "LENALIDOMIDE":
                # Sample at 4 representative days per cycle; use all 21 for superposition
                sample_ex = ex_s[ex_s["DAY_IN_CYCLE_0IDX"].isin(_LENA_SAMPLE_OFFSETS)]
            else:
                sample_ex = ex_s   # sample at every dose event (3–4 per cycle)

            if sample_ex.empty:
                continue

            # Precompute single-dose unit profile on a fine time grid
            # t_max covers last dose + 5 × Ixazomib t½  (conservative for shorter drugs)
            t_max  = float(all_dose_h.max()) + 5.0 * 228.0
            t_grid, c_grid = _unit_dose_grid(drug, row, t_max=t_max)

            # ── PK substudy sampling tier for each dose event ────────────────
            # Dense cycles: Cycles 1 & 3 — full multi-timepoint profile
            # Trough cycles: every 2nd odd cycle >3 (5, 7, 9, …) — pre-dose only
            # All other cycles: not sampled
            DENSE_CYCLES  = {1, 3}
            TROUGH_CYCLES = {c for c in range(5, 27, 2)}  # 5, 7, 9, …, 25

            for _, dose_row in sample_ex.iterrows():
                t_dose   = float(dose_row["DOSE_H"])
                d_amount = float(dose_row["EXDOSE"])
                cyc_num  = int(dose_row["VISITNUM"])
                dose_num = int(dose_row["DAY_IN_CYCLE_0IDX"]) // 7 + 1

                # Determine timepoints to sample based on SoA tier
                if cyc_num in DENSE_CYCLES:
                    post_h_use = post_h                  # full profile
                elif cyc_num in TROUGH_CYCLES and dose_num == 1:
                    post_h_use = [0.0]                   # pre-dose trough only
                else:
                    continue                              # not sampled this cycle

                eval_abs = np.array([t_dose + pt for pt in post_h_use])

                # Causal: only include doses up to this event in superposition
                mask    = all_dose_h <= t_dose
                c_total = _superpose(t_grid, c_grid,
                                     all_dose_h[mask], all_dose_mg[mask],
                                     nom_dose, eval_abs)
                c_obs, blq = add_res(c_total, drug)

                for i, pt in enumerate(post_h_use):
                    rows.append({
                        "STUDYID":   sid,
                        "DOMAIN":    "PC",
                        "USUBJID":   uid,
                        "PCTESTCD":  drug,
                        "PCTEST":    f"{drug.title()} Concentration",
                        "PCCAT":     "PLASMA",
                        "PCORRES":   f"{c_obs[i]:.4f}" if not blq[i] else f"<{p['lloq']}",
                        "PCORRESU":  p["conc_unit"],
                        "PCSTRESC":  f"{c_obs[i]:.4f}",
                        "PCSTRESN":  round(float(c_obs[i]), 4),
                        "PCSTRESU":  p["conc_unit"],
                        "PCLLOQ":    p["lloq"],
                        "BLQ":       "Y" if blq[i] else "N",
                        "PCBLFL":    "Y" if pt == 0.0 else "",
                        "VISITNUM":  cyc_num,
                        "DOSE_NUM":  dose_num,
                        "VISIT":     f"CYCLE {cyc_num} DOSE {dose_num}",
                        "EPOCH":     "TREATMENT",
                        "PCTPT":     f"{pt:.1f}H",
                        "PCTPTNUM":  round(pt, 2),
                        "DOSE_DAY":  int(dose_row["DOSE_DAY"]),
                        "DOSE_AMT":  round(d_amount, 2),
                        "PK_TIER":   "DENSE"  if cyc_num in DENSE_CYCLES else "TROUGH",
                    })

    pc = pd.DataFrame(rows)
    if len(pc):
        pc["PCSEQ"] = pc.groupby("USUBJID").cumcount() + 1
    return pc

def make_pp(adsl, study_key):
    sid = adsl["STUDYID"].iloc[0]
    pk_subj = get_pk_subj(adsl)
    rows = []; seq = 1
    for drug in ["IXAZOMIB","LENALIDOMIDE","DEXAMETHASONE"]:
        p = PK_DRUGS[drug]
        subj = pk_subj[pk_subj["ARMCD"]=="IRd"].copy() if drug=="IXAZOMIB" else pk_subj.copy()
        ip = sample_ixaz(subj) if drug=="IXAZOMIB" else sample_1cmt(drug, subj)
        subj = subj.merge(ip, on="USUBJID")
        # Use actual scheduled sampling times for NCA so that Cmax reflects
        # what the real trial NCA would compute from the clinical sample schedule
        # (e.g. Ixazomib: 0, 1, 4, 8h — not the denser 0.5h grid used for VPC).
        # Terminal timepoints are appended for AUC/t½ estimation.
        sched = SAMPLING_POST_DOSE_H[drug]
        terminal = [t for t in p["nominal_times"] if t > max(sched)]
        dt = np.array(sorted(set(sched) | set(terminal)))
        for cycle in p["cycles_assessed"]:
            for _, row in subj.iterrows():
                c_obs, _ = add_res(sim_profile(drug, row, dt), drug)
                res = nca(dt, c_obs, p["dose_mg"], p["lloq"])
                for key, val in res.items():
                    if np.isnan(val): continue
                    ptcd, ptest, pu = PP_META[key]
                    rows.append({
                        "STUDYID": sid, "DOMAIN": "PP", "USUBJID": row["USUBJID"],
                        "PPSEQ": seq, "PPTESTCD": ptcd, "PPTEST": ptest, "PPCAT": drug,
                        "PPORRES": f"{val:.4f}", "PPORRESU": pu,
                        "PPSTRESC": f"{val:.4f}", "PPSTRESN": round(float(val),4),
                        "PPSTRESU": pu, "PPSPEC": "PLASMA",
                        "VISITNUM": cycle, "VISIT": f"CYCLE {cycle} DAY 1",
                        "EPOCH": "TREATMENT", "DOSE_MG": p["dose_mg"],
                    })
                    seq += 1
    return pd.DataFrame(rows)

def make_adpc(pc, pp, adsl, study_key):
    adpc = pc.copy()
    adpc["STUDYID"] = adsl["STUDYID"].iloc[0]
    cov = [c for c in ["USUBJID","AGE","SEX","ARMCD","TRT01PN","IGTYPE",
                        "ISSSTAGE","BASE_CREACL","BASE_CREAT","RISKGR","CYTOGR"]
           if c in adsl.columns]
    adpc = adpc.merge(adsl[cov], on="USUBJID", how="left")
    for ptcd, rn in [("CMAX","CMAX_PP"),("AUCINF","AUCINF_PP")]:
        tmp = pp[pp["PPTESTCD"]==ptcd][["USUBJID","PPCAT","PPSTRESN","VISITNUM"]]\
              .rename(columns={"PPSTRESN":rn,"PPCAT":"PCTESTCD"})
        adpc = adpc.merge(tmp, on=["USUBJID","PCTESTCD","VISITNUM"], how="left")
    adpc["AVAL"]    = adpc["PCSTRESN"]
    adpc["AVALU"]   = adpc["PCSTRESU"]
    adpc["ANL01FL"] = np.where(adpc["BLQ"]=="N","Y","")
    adpc["ABLFL"]   = adpc["PCBLFL"]
    adpc["PARAM"]   = adpc["PCTEST"]
    adpc["PARAMCD"] = adpc["PCTESTCD"]
    dm = {d: PK_DRUGS[d]["dose_mg"] for d in PK_DRUGS}
    adpc["DOSE"]   = adpc["PCTESTCD"].map(dm)
    adpc["DNCONC"] = adpc["AVAL"] / adpc["DOSE"]
    adpc["LNAVAL"] = np.log(np.where(adpc["AVAL"]>0, adpc["AVAL"], np.nan))
    if "BASE_CREACL" in adpc.columns:
        adpc["RENGRP"] = pd.cut(adpc["BASE_CREACL"].fillna(90),
            bins=[0,30,60,90,999],labels=["SEVERE","MODERATE","MILD","NORMAL"])
    adpc["PKWINDFL"] = np.select(
        [adpc["PCTPTNUM"]==0, adpc["PCTPTNUM"]<=2], ["PRE-DOSE","ABSORPTION"], default="POST-DIST")
    adpc["PK_MODEL"] = adpc["PCTESTCD"].map({
        "IXAZOMIB":      "3-COMPARTMENT (Gupta et al. Clin Pharmacokinet 2017;56:1355)",
        "LENALIDOMIDE":  "1-COMPARTMENT (Chen et al. J Clin Pharmacol 2012;52:1776)",
        "DEXAMETHASONE": "1-COMPARTMENT (Frey et al. J Clin Pharmacol 1990;30:690)",
    })
    return adpc

# ── Main ─────────────────────────────────────────────────────────────────────

def generate_pk():
    ref = {
        "IXAZOMIB":     {"CMAX":"~41 ng/mL",     "THALF":"~228h (9.5d)",  "AUCINF":"~1247 ng·h/mL"},
        "LENALIDOMIDE": {"CMAX":"350-500 ng/mL",  "THALF":"~3-5h",         "AUCINF":"~2200 ng·h/mL"},
        "DEXAMETHASONE":{"CMAX":"120-250 ng/mL",  "THALF":"~5-8h",         "AUCINF":"~1800 ng·h/mL"},
    }
    for study_key in ["MM2","MM1"]:
        sdir  = os.path.join(OUT_BASE, study_key)
        adsl  = pd.read_csv(os.path.join(sdir, "adam_adsl.csv"))
        ex    = pd.read_csv(os.path.join(sdir, "sdtm_ex.csv"))
        print(f"\n{'='*65}\n  PK: {study_key} — Ixazomib: 3-COMPARTMENT [Gupta 2017]\n{'='*65}")

        print("  [1/3] PC ..."); pc = make_pc(adsl, ex, study_key)
        pc.to_csv(os.path.join(sdir, "sdtm_pc.csv"), index=False); print(f"        {len(pc):,} records")

        print("  [2/3] PP ..."); pp = make_pp(adsl, study_key)
        pp.to_csv(os.path.join(sdir, "sdtm_pp.csv"), index=False); print(f"        {len(pp):,} records")

        print("  [3/3] ADPC ..."); adpc = make_adpc(pc, pp, adsl, study_key)
        adpc.to_csv(os.path.join(sdir, "adam_adpc.csv"), index=False); print(f"        {len(adpc):,} records")

        print(f"\n  {'Drug':<18} {'Param':<10} {'Published':<18} {'Simulated'}")
        print(f"  {'-'*65}")
        for drug in ["IXAZOMIB","LENALIDOMIDE","DEXAMETHASONE"]:
            ppd = pp[pp["PPCAT"]==drug]
            if not len(ppd): continue
            for pcd, stat, unit in [("CMAX","Cmax","ng/mL"),("THALF","t½","h"),("AUCINF","AUCinf","ng·h/mL")]:
                v = ppd[ppd["PPTESTCD"]==pcd]["PPSTRESN"].median()
                s = f"{v:.1f} {unit}" if not np.isnan(v) else "NaN"
                print(f"  {drug:<18} {stat:<10} {ref[drug][pcd]:<18} {s}")
    print("\n✓ Done.")

if __name__ == "__main__":
    generate_pk()
