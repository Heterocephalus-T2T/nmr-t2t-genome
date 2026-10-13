#!/usr/bin/env python3
"""
visualize_relax_atlas_v4.py

Atlas RELAX figures, rebuilt against the v4 data layout.

WHAT CHANGED FROM part10e (v3)
------------------------------
1. Scheme names. v3 expected relax_results_{nmr_only,amr_only,bats_only,
   bathyergidae}.tsv. v4 has relax_atlas_{scheme}_vs_nonmolerat.tsv. Schemes are
   auto-discovered.

2. Categories. v4 splits ARTIFACT into ARTIFACT_K0 and ARTIFACT_KHIGH. v3's
   CAT_COLORS knew only 'ARTIFACT', so BOTH fell through to grey and every
   `(Category=='ARTIFACT').sum()` returned 0. Fixed.

3. Counting rules. K~0 is the boundary of RELAX's relaxation parameter
   (omega_T = omega_R ^ K, K in [0,50]); a significant K~0 is maximal
   relaxation, not a numerical failure. Two rules are plotted side by side:
       strict    relaxed = Category == RELAXED
       promoted  relaxed = RELAXED or (ARTIFACT_K0 and p < 0.05)
   K > 5 (ARTIFACT_KHIGH) is excluded under both.

4. Naming. v3's `eusocial_specific` = NMR & AMR & !BATS. AMR (Fukomys anselli)
   is a COOPERATIVE BREEDER, not eusocial. Renamed `social_specific`.

5. Reproducibility. Single RELAX fits are not reproducible for genes with p near
   alpha (observed: p = 0.034-0.160 across 8 fits of identical input). The
   volcano panels shade that band and the counts are labelled provisional.

Usage:
    python3 visualize_relax_atlas_v4.py
    python3 visualize_relax_atlas_v4.py --base results/10_extended_sperm_genes
    python3 visualize_relax_atlas_v4.py --rule strict     # only one rule
"""

import os, sys, glob, argparse, itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.lines as mlines

# ---------------------------------------------------------------- config
ALPHA = 0.05
RISK_LO, RISK_HI = 0.005, 0.25      # band where single-fit calls are unreliable
K_FLOOR = 1e-4                       # below this, K is at the boundary

# Which columns carry significance. --fdr switches these to the BH columns
# written by part11_fdr.py. fig_gradient_intersection is deliberately NOT
# switched: it computes its own BH from raw p on the shared gene set, and
# feeding it q-values would correct twice.
CAT_COL = 'Category'
SIG_COL = 'p_value'

SCHEME_META = {
    'nmr':        ('HetGla (eusocial)',              '#E63946'),
    'amr':        ('FukAns (cooperative breeder)',   '#457B9D'),
    'bats':       ('BatSui (solitary)',              '#2A9D8F'),
    'social':     ('Social (HetGla + FukAns)',       '#6D597A'),
    'allmolerat': ('All Bathyergidae',               '#8D6E63'),
    'highcomp':   ('High sperm competition',         '#C77DFF'),
    'lowcomp':    ('Low sperm competition',          '#FFB703'),
}
GRADIENT_ORDER = ['nmr', 'amr', 'bats']     # the sociality gradient

CAT_COLORS = {
    'RELAXED':         '#3182CE',
    'INTENSIFIED':     '#E53E3E',
    'NOT_SIGNIFICANT': '#A0AEC0',
    'ARTIFACT_K0':     '#2C7A7B',   # teal: boundary = maximal relaxation
    'ARTIFACT_KHIGH':  '#ED8936',   # orange: genuinely unreliable
    'NO_DATA':         '#E2E8F0',
}


from scipy import stats as _stats


def bh_fdr(p):
    """Benjamini-Hochberg q-values (same as gradient_intersection_test_fdr.py)."""
    p = np.asarray(p, float); n = len(p)
    if n == 0: return p
    o = np.argsort(p); q = np.empty(n)
    q[o] = p[o] * n / np.arange(1, n + 1)
    q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.minimum(q, 1.0)


def mcnemar(a_flags, b_flags):
    """Paired McNemar on two boolean arrays (same gene order)."""
    a = np.asarray(a_flags, bool); b = np.asarray(b_flags, bool)
    a_only = int((a & ~b).sum()); b_only = int((~a & b).sum())
    n = a_only + b_only
    if n == 0: return a_only, b_only, 0.0, 1.0
    chi2 = (abs(a_only - b_only) - 1) ** 2 / n
    return a_only, b_only, chi2, float(_stats.chi2.sf(chi2, 1))


def _stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


def save(fig, outdir, name):
    p = os.path.join(outdir, name)
    fig.savefig(p + '.pdf', bbox_inches='tight')
    fig.savefig(p + '.png', dpi=600, bbox_inches='tight')
    fig.savefig(p + '.svg', bbox_inches='tight')
    print(f'  saved {name}.(pdf|png|svg)')
    plt.close(fig)


# ---------------------------------------------------------------- data
def load(base, gene_set='sperm', use_fdr=False):
    """v5 layout: results/10_summary/{set}_summary.tsv holds one row per
    gene x scheme (Gene, K, p_value, LRT, Scheme, Set, Category). Split it by
    Scheme and key on the short scheme name so SCHEME_META still matches.
    With use_fdr, read {set}_summary_fdr.tsv (part11_fdr.py) instead, which
    adds q_value and Category_FDR; p_value is kept intact so the intersection
    panel can still run its own correction."""
    suffix = '_summary_fdr.tsv' if use_fdr else '_summary.tsv'
    path = os.path.join(base, f'{gene_set}{suffix}')
    if not os.path.exists(path):
        if use_fdr:
            sys.exit(f'ERROR: {path} not found. Run first:\n'
                     f'    python3 part11_fdr.py --set {gene_set}')
        return {}
    df = pd.read_csv(path, sep='\t')
    df['K'] = pd.to_numeric(df['K'], errors='coerce')
    df['p_value'] = pd.to_numeric(df['p_value'], errors='coerce')
    if use_fdr:
        df['q_value'] = pd.to_numeric(df['q_value'], errors='coerce')
    out = {}
    for scheme, d in df.groupby('Scheme'):
        s = scheme.replace('_vs_nonmolerat', '').replace('_vs_others', '')
        out[s] = d.reset_index(drop=True)
        print(f'  {s:11s} {len(d):4d} genes')
    return out


def is_relaxed(row, rule):
    if row[CAT_COL] == 'RELAXED':
        return True
    if rule == 'promoted' and row[CAT_COL] == 'ARTIFACT_K0' \
       and pd.notna(row[SIG_COL]) and row[SIG_COL] < ALPHA:
        return True
    return False


def at_risk(p):
    return pd.notna(p) and RISK_LO < p < RISK_HI


# ---------------------------------------------------------------- fig 1: volcano
def fig_volcano(data, outdir, rule):
    schemes = [s for s in GRADIENT_ORDER if s in data]
    schemes += [s for s in data if s not in schemes]
    n = len(schemes)
    ncol = min(3, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 4.6 * nrow), squeeze=False)

    for i, s in enumerate(schemes):
        ax = axes[i // ncol][i % ncol]
        d = data[s].dropna(subset=['K', 'p_value']).copy()
        d = d[d['K'] > 0]

        # K==0 genes: park at the left edge, mark with a triangle
        d['log2K'] = np.log2(d['K'].clip(lower=K_FLOOR))
        d['x'] = d['log2K'].clip(-6, 6)
        d['y'] = (-np.log10(d[SIG_COL].clip(lower=1e-300))).clip(0, 30)

        boundary = d['K'] < K_FLOOR
        colors = d[CAT_COL].map(CAT_COLORS).fillna('#A0AEC0')

        # shade the unreliable p band (defined on raw p; omitted on a q axis)
        if SIG_COL == 'p_value':
            ax.axhspan(-np.log10(RISK_HI), -np.log10(RISK_LO),
                       color='#FFF3CD', alpha=0.55, zorder=0)

        ax.scatter(d.loc[~boundary, 'x'], d.loc[~boundary, 'y'],
                   c=colors[~boundary], s=16, alpha=0.7,
                   edgecolors='white', linewidths=0.3, zorder=3)
        ax.scatter(d.loc[boundary, 'x'], d.loc[boundary, 'y'],
                   c=CAT_COLORS['ARTIFACT_K0'], s=34, marker='<',
                   edgecolors='black', linewidths=0.4, zorder=4)

        ax.axhline(-np.log10(ALPHA), color='grey', ls='--', lw=0.9)
        ax.axvline(0, color='grey', ls='--', lw=0.9)

        n_rel = int(d.apply(lambda r: is_relaxed(r, rule), axis=1).sum())
        n_int = int((d[CAT_COL] == 'INTENSIFIED').sum())
        n_k0 = int((d[CAT_COL] == 'ARTIFACT_K0').sum())
        n_kh = int((d[CAT_COL] == 'ARTIFACT_KHIGH').sum())
        n_risk = int(d['p_value'].apply(at_risk).sum())

        ax.text(0.03, 0.97,
                f'relaxed {n_rel}\nintensified {n_int}\n'
                f'K$\\approx$0 {n_k0}   K>5 {n_kh}\n'
                f'unreliable band {n_risk}',
                transform=ax.transAxes, va='top', fontsize=7.5,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85))

        lab, col = SCHEME_META.get(s, (s, '#333'))
        ax.set_title(lab, fontsize=11, fontweight='bold', color=col)
        ax.set_xlabel('log$_2$(K)   \u2190 relaxed | intensified \u2192', fontsize=9)
        ax.set_ylabel('-log$_{10}$(q)' if SIG_COL == 'q_value'
                      else '-log$_{10}$(p)', fontsize=9)
        ax.set_xlim(-6.5, 6.5)
        ax.tick_params(labelsize=8)

    for j in range(n, nrow * ncol):
        axes[j // ncol][j % ncol].axis('off')

    leg = [Patch(facecolor=CAT_COLORS['RELAXED'], label='Relaxed (K<1, p<0.05)'),
           Patch(facecolor=CAT_COLORS['INTENSIFIED'], label='Intensified (K>1, p<0.05)'),
           Patch(facecolor=CAT_COLORS['NOT_SIGNIFICANT'], label='Not significant'),
           mlines.Line2D([], [], marker='<', color='w',
                         markerfacecolor=CAT_COLORS['ARTIFACT_K0'],
                         markeredgecolor='k', markersize=8,
                         label='K$\\approx$0 (boundary: maximal relaxation)'),
           Patch(facecolor=CAT_COLORS['ARTIFACT_KHIGH'], label='K>5 (unreliable, excluded)')]
    if SIG_COL == 'p_value':
        leg.append(Patch(facecolor='#FFF3CD',
                         label=f'p in [{RISK_LO}, {RISK_HI}] \u2014 single fit unreliable'))
    fig.legend(handles=leg, loc='lower center', ncol=3, fontsize=8,
               frameon=True, bbox_to_anchor=(0.5, -0.04))
    _thr = 'BH q<0.05' if SIG_COL == 'q_value' else 'raw LRT p<0.05'
    fig.suptitle(f'RELAX across sperm-gene atlas \u2014 counting rule: {rule} '
                 f'({_thr})',
                 fontsize=13, fontweight='bold', y=1.00)
    plt.tight_layout()
    save(fig, outdir, f'atlas_volcano_{rule}')


# ---------------------------------------------------------------- fig 2: the gradient
def fig_gradient(data, outdir):
    """The money figure: relaxed % per focal species, both rules, side by side."""
    schemes = [s for s in GRADIENT_ORDER if s in data]
    if len(schemes) < 2:
        return

    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    width = 0.36
    x = np.arange(len(schemes))

    for i, rule in enumerate(('strict', 'promoted')):
        pcts, ns, ks = [], [], []
        for s in schemes:
            d = data[s]
            d = d[d[CAT_COL] != 'ARTIFACT_KHIGH']
            n = len(d)
            k = int(d.apply(lambda r: is_relaxed(r, rule), axis=1).sum())
            pcts.append(100 * k / n if n else 0)
            ns.append(n); ks.append(k)

        off = -width / 2 + i * width
        cols = [SCHEME_META[s][1] for s in schemes]
        bars = ax.bar(x + off, pcts, width,
                      color=cols, edgecolor='black', linewidth=0.6,
                      alpha=1.0 if rule == 'strict' else 0.55,
                      hatch='' if rule == 'strict' else '///')
        for b, p_, k_, n_ in zip(bars, pcts, ks, ns):
            ax.annotate(f'{p_:.1f}%\n{k_}/{n_}',
                        (b.get_x() + b.get_width() / 2, p_),
                        xytext=(0, 3), textcoords='offset points',
                        ha='center', va='bottom', fontsize=7.5)

    ax.set_xticks(x)
    ax.set_xticklabels([SCHEME_META[s][0].split(' (')[0] for s in schemes],
                       fontsize=11, fontweight='bold')
    ax.set_ylabel('Sperm genes under relaxed selection (%)', fontsize=11)
    ax.set_title('Relaxation scales with sociality\n'
                 'solid = strict (RELAXED only)   hatched = promoted (+ significant K\u22480)',
                 fontsize=11)
    ax.yaxis.grid(True, ls=':', alpha=0.4)
    ax.set_axisbelow(True)
    ax.text(0.99, 0.02,
            'Counts are single-fit; genes near \u03b1 are not reliably classified.',
            transform=ax.transAxes, ha='right', fontsize=7, style='italic', color='#666')
    plt.tight_layout()
    save(fig, outdir, 'atlas_gradient')


# ------------------------------------------------ fig 2b: INTERSECTION gradient
def _intersection_genes(data):
    """Genes analysed in ALL of nmr/amr/bats, excluding boundary-K/no-data in any.
    Matches gradient_intersection_test_fdr.py."""
    schemes = [s for s in GRADIENT_ORDER if s in data]
    if len(schemes) != 3:
        return None, None
    common = set(data[schemes[0]].Gene)
    for s in schemes[1:]:
        common &= set(data[s].Gene)
    bad = set()
    idx = {s: data[s].set_index('Gene') for s in schemes}
    for s in schemes:
        for g in list(common):
            cat = idx[s].loc[g, 'Category'] if g in idx[s].index else 'NO_DATA'
            if cat in ('ARTIFACT_K0', 'ARTIFACT_KHIGH', 'NO_DATA'):
                bad.add(g)
    genes = sorted(common - bad)
    return genes, idx


def fig_gradient_intersection(data, outdir):
    """The STATISTICALLY CORRECT gradient: intersection of the 3 schemes, FDR-
    corrected, with paired McNemar brackets. Numbers match the stats script.
    Also plots raw-vs-FDR on the same shared genes."""
    genes, idx = _intersection_genes(data)
    if not genes:
        return
    schemes = [s for s in GRADIENT_ORDER if s in data]
    n = len(genes)

    def flags(rule):
        out = {}
        for s in schemes:
            d = idx[s].loc[genes]
            K = d['K'].values; p = d['p_value'].values
            sig = (bh_fdr(p) < ALPHA) if rule == 'fdr' else (p < ALPHA)
            out[s] = sig & (K < 1)
        return out

    fdr = flags('fdr'); raw = flags('raw')

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.0))

    # --- panel A: FDR gradient with McNemar brackets ---
    ax = axes[0]
    x = np.arange(len(schemes))
    pct = [100 * fdr[s].sum() / n for s in schemes]
    cnt = [int(fdr[s].sum()) for s in schemes]
    cols = [SCHEME_META[s][1] for s in schemes]
    ax.bar(x, pct, color=cols, edgecolor='black', linewidth=1.0, width=0.6, zorder=3)
    for xi, p_, c_ in zip(x, pct, cnt):
        ax.annotate(f'{p_:.1f}%\n{c_}/{n}', (xi, p_), xytext=(0, 3),
                    textcoords='offset points', ha='center', va='bottom',
                    fontsize=9, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([SCHEME_META[s][0].split(' (')[0] for s in schemes],
                       fontsize=10, fontweight='bold')
    ax.set_ylabel('Sperm genes relaxed (%)', fontsize=11)
    ax.set_title(f'Relaxed genes across mole-rats - FDR (q<0.05), n={n} shared genes',
                 fontsize=11, fontweight='bold')
    ax.set_ylim(0, max(pct) * 1.5)
    ax.yaxis.grid(True, ls=':', alpha=0.4); ax.set_axisbelow(True)
    # McNemar brackets
    pairs = [('nmr','amr',0,1), ('amr','bats',1,2), ('nmr','bats',0,2)]
    lvl = max(pct) * 1.10; step = max(pct) * 0.13
    for i,(s1,s2,xi,xj) in enumerate(pairs):
        if s1 in fdr and s2 in fdr:
            _,_,_,pv = mcnemar(fdr[s1], fdr[s2])
            y = lvl + i*step
            ax.plot([xi,xi,xj,xj],[y,y+step*0.22,y+step*0.22,y], color='black', lw=1.0)
            ax.text((xi+xj)/2, y+step*0.24, f'{_stars(pv)} p={pv:.1e}',
                    ha='center', va='bottom', fontsize=8)

    # --- panel B: raw vs FDR on the shared genes ---
    ax = axes[1]
    w = 0.38
    raw_pct = [100*raw[s].sum()/n for s in schemes]
    fdr_pct = [100*fdr[s].sum()/n for s in schemes]
    ax.bar(x - w/2, raw_pct, w, label='raw LRT (p<0.05)', color='#b8b8b8',
           edgecolor='black', zorder=3)
    ax.bar(x + w/2, fdr_pct, w, label='FDR (q<0.05)', color=cols,
           edgecolor='black', zorder=3)
    for xi,r_,f_ in zip(x, raw_pct, fdr_pct):
        ax.annotate(f'{r_:.1f}', (xi-w/2, r_), xytext=(0,2), textcoords='offset points',
                    ha='center', fontsize=8)
        ax.annotate(f'{f_:.1f}', (xi+w/2, f_), xytext=(0,2), textcoords='offset points',
                    ha='center', fontsize=8, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([SCHEME_META[s][0].split(' (')[0] for s in schemes],
                       fontsize=10, fontweight='bold')
    ax.set_ylabel('Sperm genes relaxed (%)', fontsize=11)
    ax.set_title(f'Raw LRT vs FDR on the shared {n} genes', fontsize=11, fontweight='bold')
    ax.legend(frameon=False, fontsize=9)
    ax.set_ylim(0, max(raw_pct) * 1.25)
    ax.yaxis.grid(True, ls=':', alpha=0.4); ax.set_axisbelow(True)

    plt.tight_layout()
    save(fig, outdir, 'atlas_gradient_intersection')

    # save numbers so tables/text match the figure
    rows = []
    for s in schemes:
        rows.append({'scheme': s, 'n_shared': n,
                     'raw_relaxed': int(raw[s].sum()), 'raw_pct': round(100*raw[s].sum()/n,1),
                     'fdr_relaxed': int(fdr[s].sum()), 'fdr_pct': round(100*fdr[s].sum()/n,1)})
    pd.DataFrame(rows).to_csv(os.path.join(outdir,'gradient_intersection_numbers.tsv'),
                             sep='\t', index=False)
    mc = []
    for s1,s2 in itertools.combinations(schemes,2):
        for rule,fl in (('fdr',fdr),('raw',raw)):
            a,b,chi2,pv = mcnemar(fl[s1], fl[s2])
            mc.append({'rule':rule,'pair':f'{s1}_vs_{s2}',
                       f'{s1}_only':a,f'{s2}_only':b,'chi2':round(chi2,3),
                       'p_value':pv,'sig':_stars(pv)})
    pd.DataFrame(mc).to_csv(os.path.join(outdir,'gradient_intersection_mcnemar.tsv'),
                           sep='\t', index=False)
    print(f'  intersection gradient: n={n}; '
          + ', '.join(f'{s}={100*fdr[s].sum()/n:.1f}%' for s in schemes))


# ---------------------------------------------------------------- fig 3: overlap
def fig_overlap(data, outdir, rule):
    schemes = [s for s in GRADIENT_ORDER if s in data]
    if len(schemes) != 3:
        return
    sets = {}
    for s in schemes:
        d = data[s]
        d = d[d[CAT_COL] != 'ARTIFACT_KHIGH']
        sets[s] = set(d[d.apply(lambda r: is_relaxed(r, rule), axis=1)].Gene)

    regions, labels = [], []
    a, b, c = schemes
    # species display names (HetGla / FukAns / BatSui) instead of raw keys nmr/amr/bats
    nm = {s: SCHEME_META[s][0].split(' (')[0] for s in schemes}
    combos = [
        (f'{nm[a]} only', sets[a] - sets[b] - sets[c]),
        (f'{nm[b]} only', sets[b] - sets[a] - sets[c]),
        (f'{nm[c]} only', sets[c] - sets[a] - sets[b]),
        (f'{nm[a]} & {nm[b]}', (sets[a] & sets[b]) - sets[c]),
        (f'{nm[a]} & {nm[c]}', (sets[a] & sets[c]) - sets[b]),
        (f'{nm[b]} & {nm[c]}', (sets[b] & sets[c]) - sets[a]),
        ('all three', sets[a] & sets[b] & sets[c]),
    ]
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    names = [n for n, _ in combos]
    vals = [len(v) for _, v in combos]
    cols = [SCHEME_META[a][1], SCHEME_META[b][1], SCHEME_META[c][1],
            '#8899aa', '#8899aa', '#8899aa', '#444']
    bars = ax.bar(range(len(vals)), vals, color=cols, edgecolor='black', linewidth=0.5)
    for bar, v in zip(bars, vals):
        ax.annotate(str(v), (bar.get_x() + bar.get_width() / 2, v),
                    xytext=(0, 2), textcoords='offset points',
                    ha='center', fontsize=9, fontweight='bold')
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=30, ha='right', fontsize=9)
    ax.set_ylabel('Relaxed genes', fontsize=11)
    ax.set_title(f'Relaxed-gene overlap between focal species ({rule})', fontsize=11)
    ax.yaxis.grid(True, ls=':', alpha=0.4); ax.set_axisbelow(True)
    plt.tight_layout()
    save(fig, outdir, f'atlas_overlap_{rule}')

    # write the gene lists
    with open(os.path.join(outdir, f'overlap_genes_{rule}.txt'), 'w') as f:
        for name, genes in combos:
            f.write(f'## {name} (n={len(genes)})\n')
            f.write('\n'.join(sorted(genes)) + '\n\n')


# ---------------------------------------------------------------- summary table
def summary_table(data, outdir):
    schemes = [s for s in GRADIENT_ORDER if s in data]
    all_genes = sorted(set().union(*(set(d.Gene) for d in data.values())))
    rows = []
    for g in all_genes:
        r = {'Gene': g}
        for s, d in data.items():
            m = d[d.Gene == g]
            if len(m) == 1:
                r[f'K_{s}'] = m.iloc[0]['K']
                r[f'p_{s}'] = m.iloc[0]['p_value']
                r[f'cat_{s}'] = m.iloc[0][CAT_COL]
                if SIG_COL == 'q_value':
                    r[f'q_{s}'] = m.iloc[0]['q_value']
                for rule in ('strict', 'promoted'):
                    r[f'relaxed_{rule}_{s}'] = is_relaxed(m.iloc[0], rule)
                r[f'at_risk_{s}'] = at_risk(m.iloc[0]['p_value'])
            else:
                r[f'cat_{s}'] = 'NO_DATA'
        # NOTE: NMR is eusocial, AMR is a COOPERATIVE BREEDER. NMR&AMR = "social",
        # not "eusocial". v3 mislabelled this.
        if all(s in data for s in GRADIENT_ORDER):
            for rule in ('strict', 'promoted'):
                nmr = r.get(f'relaxed_{rule}_nmr', False)
                amr = r.get(f'relaxed_{rule}_amr', False)
                bats = r.get(f'relaxed_{rule}_bats', False)
                r[f'social_specific_{rule}'] = bool(nmr and amr and not bats)
                r[f'hetgla_only_{rule}'] = bool(nmr and not amr and not bats)
        rows.append(r)
    df = pd.DataFrame(rows)
    out = os.path.join(outdir, 'atlas_relax_summary.tsv')
    df.to_csv(out, sep='\t', index=False)
    print(f'  saved atlas_relax_summary.tsv  ({len(df)} genes)')

    for rule in ('strict', 'promoted'):
        col = f'social_specific_{rule}'
        if col in df:
            print(f'    social-specific ({rule}): {int(df[col].sum())}   '
                  f'HetGla-only: {int(df[f"hetgla_only_{rule}"].sum())}')
    return df


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--base', default='results/10_summary')
    ap.add_argument('--set', dest='gene_set', default='sperm',
                    help='which summary to plot: sperm | ddr | targets17')
    ap.add_argument('--rule', choices=['strict', 'promoted', 'both'], default='both')
    ap.add_argument('--fdr', action='store_true',
                    help='use Benjamini-Hochberg q-values from '
                         '{set}_summary_fdr.tsv (written by part11_fdr.py) '
                         'instead of raw LRT p-values; output goes to a '
                         'separate atlas_{set}_fdr directory')
    args = ap.parse_args()

    global CAT_COL, SIG_COL
    if args.fdr:
        CAT_COL, SIG_COL = 'Category_FDR', 'q_value'

    _sfx = '_fdr' if args.fdr else ''
    outdir = os.path.join('results/11_figures', f'atlas_{args.gene_set}{_sfx}')
    os.makedirs(outdir, exist_ok=True)

    print('loading atlas results...')
    data = load(args.base, args.gene_set, args.fdr)
    if not data:
        sys.exit(f'no {args.gene_set}_summary.tsv under {args.base}')

    rules = ['strict', 'promoted'] if args.rule == 'both' else [args.rule]

    print('\nfigures:')
    for rule in rules:
        fig_volcano(data, outdir, rule)
        fig_overlap(data, outdir, rule)
    fig_gradient(data, outdir)                 # all-genes barplot (full set per scheme)
    fig_gradient_intersection(data, outdir)   # intersection + FDR + McNemar (the stats)
    summary_table(data, outdir)

    print(f'\n-> {outdir}/')
    print('\nREMINDER: counts are single-fit. Genes with p in '
          f'[{RISK_LO}, {RISK_HI}] are not reliably classified from one RELAX run.')


if __name__ == '__main__':
    main()
