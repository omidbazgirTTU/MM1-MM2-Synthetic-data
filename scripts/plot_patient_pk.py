import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from scipy.integrate import solve_ivp
import warnings; warnings.filterwarnings('ignore')

# ─── Load data ────────────────────────────────────────────────────────────────
MM2  = '/mnt/user-data/outputs/tourmaline_synthetic/MM2'
adsl = pd.read_csv(f'{MM2}/adam_adsl.csv')
pc   = pd.read_csv(f'{MM2}/sdtm_pc.csv')
ds   = pd.read_csv(f'{MM2}/sdtm_ds.csv')
ae   = pd.read_csv(f'{MM2}/sdtm_ae.csv')
lb   = pd.read_csv(f'{MM2}/sdtm_lb.csv')

SUBJ = 'TOURMALINE-MM2-0094'
si   = adsl[adsl['USUBJID']==SUBJ].iloc[0]
spc  = pc[pc['USUBJID']==SUBJ]
sae  = ae[ae['USUBJID']==SUBJ].copy()
sds  = ds[ds['USUBJID']==SUBJ]
slb  = lb[lb['USUBJID']==SUBJ]

pfs  = sds[sds['PARAMCD']=='PFS'].iloc[0]
oss  = sds[sds['PARAMCD']=='OS'].iloc[0]
pfs_day = float(pfs['DSSTDY'])   # day of PFS event
os_day  = float(oss['DSSTDY'])   # day of OS event

# ─── PK models ────────────────────────────────────────────────────────────────
def pk_3cmt(t, dose, Ka, ALAG, CL, V2, Q3, V3, Q4, V4, F=0.58):
    ke=CL/V2; k23=Q3/V2; k32=Q3/V3; k24=Q4/V2; k42=Q4/V4
    def odes(t_, y):
        d,c,p1,p2=y
        ab=Ka*d if t_>=ALAG else 0.
        return [-Ka*d if t_>=ALAG else 0, F*ab-(ke+k23+k24)*c+k32*p1+k42*p2,
                k23*c-k32*p1, k24*c-k42*p2]
    sol=solve_ivp(odes,[0,max(t[-1]+0.01,0.02)],[dose,0,0,0],t_eval=t,
                  method='RK45',rtol=1e-7,atol=1e-10)
    return np.maximum(sol.y[1]/V2*1000, 0.)

def pk_1cmt(t, dose, Ka, CL_F, Vd_F):
    ke=CL_F/Vd_F
    if abs(Ka-ke)<1e-6: ke*=1.001
    return np.maximum((dose*Ka)/(Vd_F*(Ka-ke))*(np.exp(-ke*t)-np.exp(-Ka*t))*1000, 0.)

# Individual params — estimated to approximate observed Cmax ~106 ng/mL ixazomib
# Use lower CL (=1.0 L/h) to match this high-exposure patient
IND = {
    'IXAZOMIB':     dict(Ka=0.52, ALAG=0.12, CL=1.0, V2=12.0, Q3=8.0, V3=190., Q4=0.45, V4=341., F=0.58),
    'LENALIDOMIDE': dict(Ka=0.9,  CL_F=7.0,  Vd_F=75.),
    'DEXAMETHASONE':dict(Ka=1.4,  CL_F=13.5, Vd_F=160.),
}
DOSES   = {'IXAZOMIB':4., 'LENALIDOMIDE':25., 'DEXAMETHASONE':40.}
SCHED   = {'IXAZOMIB':[0,7,14], 'LENALIDOMIDE':list(range(21)), 'DEXAMETHASONE':[0,7,14,21]}
COLORS  = {'IXAZOMIB':'#E05C5C', 'LENALIDOMIDE':'#4D9FEC', 'DEXAMETHASONE':'#F5A623'}
LEND_C, DEX_C = '#4D9FEC', '#F5A623'

n_cycles = int(np.ceil(pfs_day / 28)) + 1
t_end_day = pfs_day + 56   # show 2 cycles past PFS

def superpose(drug, t_hours):
    conc = np.zeros_like(t_hours, dtype=float)
    p    = IND[drug]; dose=DOSES[drug]
    for cyc in range(n_cycles):
        for day in SCHED[drug]:
            t_dose = (cyc*28 + day)*24
            tr     = t_hours - t_dose
            mask   = tr > 0
            if not mask.any(): continue
            if drug=='IXAZOMIB':
                c = pk_3cmt(tr[mask], dose, **p)
            else:
                c = pk_1cmt(tr[mask], dose, **p)
            conc[mask] += c
    return conc

# Fine grid for smooth curves (0.25h steps)
t_h  = np.arange(0, t_end_day*24, 0.25)
t_d  = t_h / 24

print("Simulating profiles...")
curves = {d: superpose(d, t_h) for d in IND}
for d,c in curves.items():
    print(f"  {d}: Cmax={c.max():.1f} ng/mL")

# ─── Observed samples → absolute time (days) ─────────────────────────────────
def get_obs(drug):
    rows = spc[spc['PCTESTCD']==drug]
    ts, cs, bq = [], [], []
    for _,r in rows.iterrows():
        t_abs = (r['VISITNUM']-1)*28 + r['PCTPTNUM']/24
        ts.append(t_abs); cs.append(r['PCSTRESN']); bq.append(r['BLQ']=='Y')
    return np.array(ts), np.array(cs), np.array(bq)

# ─── Dose tick times (days) ───────────────────────────────────────────────────
def dose_days(drug):
    return [cyc*28+day for cyc in range(n_cycles) for day in SCHED[drug]
            if cyc*28+day <= t_end_day]

# ─── Lab data ─────────────────────────────────────────────────────────────────
mp  = slb[slb['LBTESTCD']=='SPEP_MPROT'][['LBDY','LBSTRESN']].dropna().sort_values('LBDY')
hgb = slb[slb['LBTESTCD'].isin(['HGB','HGBHGB'])][['LBDY','LBSTRESN']].dropna().sort_values('LBDY')
plt_df = slb[slb['LBTESTCD']=='PLT'][['LBDY','LBSTRESN']].dropna().sort_values('LBDY')

# ═══════════════════════════════════════════════════════════════════════════════
#  FIGURE
# ═══════════════════════════════════════════════════════════════════════════════
BG  = '#0E1117'
PAN = '#13161E'
GRD = '#1E2233'
TXT = '#C8CCE0'
DIM = '#5A5F7A'

fig = plt.figure(figsize=(20, 24), facecolor=BG)
gs  = gridspec.GridSpec(5, 1, figure=fig,
                         height_ratios=[0.22, 2.2, 2.0, 1.4, 1.4],
                         hspace=0.0)

ax0 = fig.add_subplot(gs[0])   # header
ax1 = fig.add_subplot(gs[1])   # ixazomib PK
ax2 = fig.add_subplot(gs[2])   # len + dex PK
ax3 = fig.add_subplot(gs[3])   # biomarkers
ax4 = fig.add_subplot(gs[4])   # event / dose timeline

for ax in [ax0,ax1,ax2,ax3,ax4]:
    ax.set_facecolor(PAN)
    for sp in ax.spines.values():
        sp.set_color(GRD); sp.set_linewidth(0.8)

# helper — draw shared vertical event lines
def draw_events(ax, xscale='days'):
    x_pfs = pfs_day if xscale=='days' else pfs_day/28
    x_os  = os_day  if xscale=='days' else os_day/28
    ax.axvline(x_pfs, color='#FF5555', lw=1.8, ls='--', zorder=10)
    ax.axvline(x_os,  color='#9999FF', lw=1.5, ls=':',  zorder=10)

def grid(ax):
    ax.grid(True, color=GRD, lw=0.5, ls='--', zorder=0)
    ax.tick_params(colors=DIM, which='both', labelsize=9)
    ax.xaxis.label.set_color(TXT); ax.yaxis.label.set_color(TXT)

# ── [0] Header ────────────────────────────────────────────────────────────────
ax0.axis('off')
pfs_str = f"{pfs['AVAL']:.1f} mo ({'Event' if pfs['CNSR']==0 else 'Censored'})"
os_str  = f"{oss['AVAL']:.1f} mo ({'Event' if oss['CNSR']==0 else 'Censored'})"
header  = (
    f"Patient PK / PD Profile  —  {SUBJ}\n"
    f"Age {int(si['AGE'])} {si['SEX']}  ·  {si['IGTYPE']} Myeloma  ·  "
    f"ISS Stage {int(si['ISSSTAGE'])}  ·  Arm: {si['ARMCD']}  ·  "
    f"CrCL: {int(si['BASE_CREACL'])} mL/min  ·  "
    f"PFS: {pfs_str}  ·  OS: {os_str}"
)
ax0.text(0.5, 0.52, header, transform=ax0.transAxes,
         ha='center', va='center', color=TXT, fontsize=12,
         fontweight='bold', linespacing=1.7)

# ── [1] Ixazomib PK (log scale) ───────────────────────────────────────────────
col = COLORS['IXAZOMIB']
c1  = curves['IXAZOMIB']
t_obs1, c_obs1, blq1 = get_obs('IXAZOMIB')
dd1 = dose_days('IXAZOMIB')

# Dose tick verticals (faint)
for dt in dd1:
    ax1.axvline(dt, color=col, lw=0.35, alpha=0.20, zorder=1)

ax1.semilogy(t_d, np.where(c1>0.2, c1, np.nan),
             color=col, lw=2.0, alpha=0.92, zorder=4, label='Predicted curve (3-cmt)')
ax1.fill_between(t_d, 0.2, np.where(c1>0.2, c1, np.nan), color=col, alpha=0.07)

# LLOQ line
ax1.axhline(0.5, color=DIM, lw=0.9, ls=':', alpha=0.7, zorder=3)
ax1.text(t_end_day*0.99, 0.52, 'LLOQ', ha='right', va='bottom',
         color=DIM, fontsize=7.5)

# Observed
ax1.scatter(t_obs1[~blq1], c_obs1[~blq1],
            color=col, s=75, zorder=8, edgecolors='white', lw=1.0,
            label='Observed (sparse sampling)')
if blq1.any():
    ax1.scatter(t_obs1[blq1], [0.35]*blq1.sum(),
                color='#AAAAAA', s=50, marker='v', zorder=7, alpha=0.8,
                label='BLQ (plotted at LLOQ/2)')

# Dose amount annotations (first 4 cycles)
for i, dt in enumerate(dd1[:12]):
    if i % 3 == 0:  # only Day 1 of each cycle
        cyc = i // 3 + 1
        ax1.annotate(f'C{cyc}\n4mg', xy=(dt, 0.25), xytext=(dt+0.4, 0.21),
                     fontsize=6.5, color=col, alpha=0.7,
                     arrowprops=dict(arrowstyle='-', color=col, alpha=0.3, lw=0.7))

draw_events(ax1)
ax1.set_xlim(0, t_end_day)
ax1.set_ylim(0.18, 400)
ax1.set_ylabel('Conc. (ng/mL)\n[log scale]', color=TXT, fontsize=10)
ax1.set_title('Ixazomib  ·  3-Compartment Population PK  [Gupta et al., Clin Pharmacokinet 2017;56:1355]',
              color='#8899DD', fontsize=9.5, loc='left', pad=6)
grid(ax1)
ax1.tick_params(labelbottom=False)

# Cycle number labels
for cyc in range(min(n_cycles, 12)):
    mid = cyc*28 + 14
    if mid < pfs_day:
        ax1.text(mid, 0.23, f'C{cyc+1}', ha='center', fontsize=7.5,
                 color='#3A3D55', fontweight='bold')

leg1 = ax1.legend(loc='upper right', fontsize=8.5, framealpha=0.45,
                  facecolor='#1A1E30', edgecolor='#33365A', labelcolor=TXT)

# ── [2] Lenalidomide + Dexamethasone ─────────────────────────────────────────
ax2r = ax2.twinx()
ax2r.set_facecolor(PAN)
ax2r.spines['right'].set_color(GRD)

for drug, ax_, col_, side in [
    ('LENALIDOMIDE',  ax2,  LEND_C, 'left'),
    ('DEXAMETHASONE', ax2r, DEX_C,  'right')]:

    c_ = curves[drug]
    t_o, c_o, bq_o = get_obs(drug)
    dd_ = dose_days(drug)

    for dt in dd_:
        ax_.axvline(dt, color=col_, lw=0.25, alpha=0.12, zorder=1)

    ax_.plot(t_d, c_, color=col_, lw=1.8, alpha=0.88, zorder=4)
    ax_.fill_between(t_d, 0, c_, color=col_, alpha=0.07)
    ax_.scatter(t_o, c_o, color=col_, s=65, zorder=8,
                edgecolors='white', lw=0.9)

ax2.set_xlim(0, t_end_day)
ax2.set_ylim(bottom=0)
ax2r.set_ylim(bottom=0)
draw_events(ax2)
grid(ax2); ax2r.tick_params(colors=DIM, labelsize=9)

ax2.set_ylabel('Lenalidomide (ng/mL)', color=LEND_C, fontsize=10)
ax2r.set_ylabel('Dexamethasone (ng/mL)', color=DEX_C, fontsize=10)
ax2.yaxis.label.set_color(LEND_C)
ax2r.yaxis.label.set_color(DEX_C)
ax2.set_title('Lenalidomide  ·  1-Compartment  &  Dexamethasone  ·  1-Compartment  [1-cmt model]',
              color='#8899DD', fontsize=9.5, loc='left', pad=6)
ax2.tick_params(labelbottom=False)

# Dose labels
ax2.text(0.01, 0.95, 'Lenalidomide 25 mg  Days 1–21',
         transform=ax2.transAxes, color=LEND_C, fontsize=8.5, va='top', alpha=0.9)
ax2.text(0.01, 0.87, 'Dexamethasone 40 mg  Days 1, 8, 15, 22',
         transform=ax2.transAxes, color=DEX_C, fontsize=8.5, va='top', alpha=0.9)

l_patch = mpatches.Patch(color=LEND_C, label='Lenalidomide — predicted')
d_patch = mpatches.Patch(color=DEX_C,  label='Dexamethasone — predicted')
obs_dot = Line2D([0],[0], marker='o', color='w', markerfacecolor='white',
                 markersize=7, lw=0, label='Observed points')
ax2.legend(handles=[l_patch, d_patch, obs_dot], loc='upper right',
           fontsize=8.5, framealpha=0.45, facecolor='#1A1E30',
           edgecolor='#33365A', labelcolor=TXT)

# ── [3] Disease Biomarkers ────────────────────────────────────────────────────
ax3r = ax3.twinx(); ax3r.set_facecolor(PAN)
ax3r.spines['right'].set_color(GRD)

if len(mp):
    t_mp = mp['LBDY'].values
    v_mp = mp['LBSTRESN'].values
    ax3.plot(t_mp, v_mp, color='#50C878', lw=2.0, marker='o',
             ms=6, markeredgecolor='white', mew=0.8, zorder=5, label='M-Protein (g/L)')
    ax3.fill_between(t_mp, 0, v_mp, color='#50C878', alpha=0.12)
    ax3.axhspan(0, 5,  color='#50C878', alpha=0.05)  # response zone

if len(hgb):
    t_hb = hgb['LBDY'].values
    v_hb = hgb['LBSTRESN'].values
    ax3r.plot(t_hb, v_hb, color='#FFD166', lw=2.0, marker='s',
              ms=5, markeredgecolor='white', mew=0.8, zorder=5, label='Hemoglobin (g/L)')
    ax3r.axhspan(120, 160, color='#FFD166', alpha=0.05)
    ax3r.text(t_end_day*0.995, 125, 'Normal\nrange', ha='right', va='bottom',
              color='#FFD166', fontsize=7, alpha=0.6)

draw_events(ax3)
ax3.set_xlim(0, t_end_day)
ax3.set_ylim(bottom=0)
ax3r.set_ylim(bottom=0)
grid(ax3); ax3r.tick_params(colors=DIM, labelsize=9)

ax3.set_ylabel('M-Protein (g/L)', color='#50C878', fontsize=10)
ax3r.set_ylabel('Hemoglobin (g/L)', color='#FFD166', fontsize=10)
ax3.yaxis.label.set_color('#50C878')
ax3r.yaxis.label.set_color('#FFD166')
ax3.set_title('Disease Biomarkers', color='#8899DD', fontsize=9.5, loc='left', pad=6)
ax3.tick_params(labelbottom=False)

mp_line  = Line2D([0],[0],color='#50C878',lw=2,marker='o',ms=5,label='M-Protein (g/L)')
hb_line  = Line2D([0],[0],color='#FFD166',lw=2,marker='s',ms=5,label='Hemoglobin (g/L)')
pfs_line = Line2D([0],[0],color='#FF5555',lw=2,ls='--',
                  label=f'PFS Event  (Day {int(pfs_day)}, {pfs["AVAL"]:.1f} mo)')
os_line  = Line2D([0],[0],color='#9999FF',lw=1.5,ls=':',
                  label=f'OS Event   (Day {int(os_day)}, {oss["AVAL"]:.1f} mo)')
ax3.legend(handles=[mp_line,hb_line,pfs_line,os_line],
           loc='upper right', fontsize=8.5, framealpha=0.45,
           facecolor='#1A1E30', edgecolor='#33365A', labelcolor=TXT)

# ── [4] Dose & Event Timeline ─────────────────────────────────────────────────
ax4.set_xlim(0, t_end_day)
ax4.set_ylim(-0.5, 6.0)
ax4.set_yticks([])

# ── Treatment bar
ax4.broken_barh([(0, pfs_day)], (4.55, 0.7),
                facecolors='#1E3A5F', edgecolors='#3A6EA8', lw=1.3)
ax4.text(pfs_day/2, 4.9, f'IRd  ({int(pfs_day)} days)', ha='center', va='center',
         color='#7AAED4', fontsize=9.5, fontweight='bold')

# ── Dose rows (triangles = dose administration)
DOC_Y  = {'IXAZOMIB':3.65, 'DEXAMETHASONE':3.05, 'LENALIDOMIDE':2.45}
DOC_LB = {'IXAZOMIB':'Ixazomib 4 mg', 'DEXAMETHASONE':'Dexamethasone 40 mg',
           'LENALIDOMIDE':'Lenalidomide 25 mg'}

for drug, y_pos in DOC_Y.items():
    col_ = COLORS[drug]
    dd_  = [d for d in dose_days(drug) if d <= pfs_day]
    # For lenalidomide show only cycle start for clarity
    if drug == 'LENALIDOMIDE':
        dd_ = [d for d in dd_ if d % 28 == 0]
        extra = ' (cycle start shown)'
    else:
        extra = ''
    ax4.scatter(dd_, [y_pos]*len(dd_), marker='v', color=col_,
                s=50 if drug != 'LENALIDOMIDE' else 35, alpha=0.85, zorder=5)
    ax4.text(-2, y_pos, f'{DOC_LB[drug]}{extra}',
             ha='right', va='center', color=col_, fontsize=8.0, clip_on=False)
    # Draw dose amount box at first dose
    ax4.annotate(f'{int(DOSES[drug])} mg',
                 xy=(0, y_pos+0.18), fontsize=7, color=col_,
                 alpha=0.8, ha='left', va='bottom')

# ── Adverse Events
AE_COLORS = {
    'Diarrhea':      '#FFB347',
    'Heart Failure': '#FF6B6B',
    'Neutropenia':   '#CC88FF',
    'Thrombocytopenia':'#DDA0DD',
    'Rash':          '#98FB98',
    'Acute Renal Failure':'#FF4444',
    'Default':       '#AAAAAA',
}
AE_GRADES = {'1':'●','2':'●●','3':'●●●','4':'●●●●'}

y_ae_base = 1.15
for i, (_, row) in enumerate(sae.iterrows()):
    t_start = float(row['AEDY'])
    t_end_  = float(row.get('AEENDY', t_start+7))
    dur     = max(t_end_ - t_start, 3)
    name    = row['AEDECOD']
    grade   = str(int(row.get('AETOXGR', 2)))
    serious = row.get('AESER','N') == 'Y'
    col_ae  = AE_COLORS.get(name, AE_COLORS['Default'])
    y_ae    = y_ae_base + i*0.52

    ax4.broken_barh([(t_start, dur)], (y_ae-0.2, 0.4),
                    facecolors=col_ae, alpha=0.75, edgecolors='none', zorder=4)
    label_txt = f'{name}  G{grade}{"  ⚠ SERIOUS" if serious else ""}'
    ax4.text(t_start + dur/2, y_ae, label_txt, ha='center', va='center',
             fontsize=7.5, color='#111', fontweight='bold', zorder=6)

    # Action taken annotation
    action = str(row.get('AEACN',''))
    if action:
        ax4.text(t_end_ + 2, y_ae, action.title(), ha='left', va='center',
                 fontsize=6.5, color=col_ae, alpha=0.8)

ax4.text(-2, y_ae_base + (len(sae)-1)*0.52/2, 'Adverse\nEvents',
         ha='right', va='center', color=TXT, fontsize=8, clip_on=False)

# ── Event lines + labels
ax4.axvline(pfs_day, color='#FF5555', lw=2.2, ls='--', zorder=10)
ax4.axvline(os_day,  color='#9999FF', lw=1.8, ls=':',  zorder=10)

ax4.text(pfs_day+2, 5.55, f'PFS Event\nDay {int(pfs_day)}',
         color='#FF5555', fontsize=8.5, fontweight='bold', va='top')
ax4.text(os_day+2,  5.0,  f'OS Event\nDay {int(os_day)}',
         color='#9999FF', fontsize=8.5, fontweight='bold', va='top')

# ── Cycle tick marks
for cyc in range(n_cycles+2):
    x = cyc*28
    if x <= t_end_day:
        ax4.axvline(x, color='#1E2233', lw=0.8, ls='-', zorder=1, alpha=0.7)
        ax4.text(x+14, -0.3, f'C{cyc+1}', ha='center', va='top',
                 fontsize=7.5, color='#44475A')

ax4.set_xlabel('Study Day', color=TXT, fontsize=11)
ax4.set_title('Dosing Schedule  ·  Adverse Events  ·  Clinical Outcomes',
              color='#8899DD', fontsize=9.5, loc='left', pad=6)
grid(ax4)

# ── X-axis sync: all panels share study-day scale ─────────────────────────────
for ax_ in [ax1, ax2, ax3, ax4]:
    ax_.xaxis.set_major_locator(ticker.MultipleLocator(28))
    ax_.xaxis.set_major_formatter(ticker.FuncFormatter(
        lambda x, _: f'Day {int(x)}\n(Mo {x/28:.0f})' if x > 0 else 'Day 0'))
    ax_.tick_params(axis='x', colors=DIM, labelsize=8.5)
    ax_.set_xlim(0, t_end_day)

for ax_ in [ax1, ax2, ax3]:
    ax_.tick_params(labelbottom=False)

# ── Outer styling ──────────────────────────────────────────────────────────────
for ax_ in [ax0,ax1,ax2,ax3,ax4]:
    ax_.spines['top'].set_visible(False)

plt.tight_layout(rect=[0.06, 0, 1, 1])

out = '/mnt/user-data/outputs/patient_pk_profile.png'
fig.savefig(out, dpi=160, bbox_inches='tight', facecolor=BG)
print(f'Saved: {out}')
