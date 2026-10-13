#!/usr/bin/env python3
"""
visualize_17gene_relax_v7.py

Focal 17-gene RELAX figures. Keeps the v5 layout (horizontal dot plot, linear
0-5 K axis, K>5 capped with arrow, K~0 pinned near zero), with four corrections.

1. SCHEME NAMES. v5 expected results/09_relax_results/{nmr_only,amr_only,...}.
   v4 writes {nmr,amr,bats,social,highcomp,lowcomp}_vs_nonmolerat. Auto-detected.

2. SOCIAL SYSTEMS. v5 line 61 labelled AMR "eusocial" and the NMR+AMR scheme
   "Eusocial". Fukomys anselli is a COOPERATIVE BREEDER, not eusocial. That
   error would have printed on a figure. Corrected here:
       NMR  Heterocephalus glaber  eusocial
       AMR  Fukomys anselli        cooperative breeder
       BATS Bathyergus suillus     solitary

3. SIGNIFICANCE - TWO TIERS (v7). Stars mark the raw per-gene LRT p-value
   (* p<0.05, ** p<0.01, *** p<0.001), as in every previous version. A dagger
   after the stars marks survival of Benjamini-Hochberg correction at q < 0.05,
   computed WITHIN each scheme across the testable genes only, matching the
   statistics reported in the results text and in the atlas. Both tiers appear
   on the figure and in the exported TSV (q_value column).

4. BOUNDARY K. K ~ 0 is the lower bound of RELAX's parameter range
   (omega_T = omega_R ^ K, K in [0,50]). When the MLE sits on a boundary the
   likelihood-ratio statistic does not follow its nominal chi-square, so the
   reported p-value is invalid. Such genes are EXCLUDED from the hypothesis set
   (and therefore from the BH denominator) and drawn as open markers labelled
   "boundary - excluded", not as "artifacts" with stars.

5. LEGACY GENE NAMES. Files on disk carry current HGNC symbols; figures and
   tables print the symbols the candidate list was designed with
   (GARIN3 -> FAM71B, GARIN1A -> FAM71F2).

Usage:
    python3 visualize_17gene_relax_v7.py [output_dir]
"""

import os, sys, json, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

RELAX_DIR = 'results/09_relax/targets17'
DEFAULT_OUT = 'results/11_figures/targets17'

K_ZERO = 1e-4
K_HIGH = 5.0
ALPHA = 0.05

# current HGNC symbol on disk -> legacy symbol used in the manuscript
DISPLAY_NAME = {'GARIN3': 'FAM71B', 'GARIN1A': 'FAM71F2'}

SCHEME_NAMES = {
    'nmr_vs_nonmolerat':      'HetGla - naked mole-rat (eusocial)',
    'amr_vs_nonmolerat':      'FukAns - Ansell\u2019s mole-rat (cooperative breeder)',
    'bats_vs_nonmolerat':     'BatSui - Cape dune mole-rat (solitary)',
    'social_vs_nonmolerat':   'Eusocial + coop. breeder (HetGla + FukAns) vs non-mole-rats',
    'highcomp_vs_nonmolerat': 'High sperm competition vs other non-mole-rats',
    'lowcomp_vs_nonmolerat':  'Low sperm competition vs other non-mole-rats',
    'allmolerat_vs_nonmolerat': 'All Bathyergidae vs non-mole-rats',
}
PRIORITY = ['nmr_vs_nonmolerat', 'amr_vs_nonmolerat', 'bats_vs_nonmolerat',
            'social_vs_nonmolerat', 'allmolerat_vs_nonmolerat',
            'highcomp_vs_nonmolerat', 'lowcomp_vs_nonmolerat']

# colours aligned with the atlas figures
COLORS = {
    'relaxed':         '#3182CE',
    'intensified':     '#E53E3E',
    'not_significant': '#A0AEC0',
    'boundary':        '#2C7A7B',   # K~0, excluded (teal, as in atlas)
    'k_high':          '#ED8936',   # K>5, excluded
    'failed':          '#D3D3D3',
}


def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = p[o] * n / np.arange(1, n + 1)
    q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.minimum(q, 1.0)


def parse_json(path):
    try:
        d = json.load(open(path))
        tr = d.get('test results', {})
        K = tr.get('relaxation or intensification parameter')
        p = tr.get('p-value')
        return {'K': K, 'p_value': p, 'LRT': tr.get('LRT')}
    except Exception:
        return {'K': None, 'p_value': None, 'LRT': None}


def load_scheme(scheme):
    d = os.path.join(RELAX_DIR, scheme)
    rows = []
    for jf in sorted(glob.glob(os.path.join(d, 'json', '*_RELAX.json'))):
        r = parse_json(jf)
        _g = os.path.basename(jf).replace('_RELAX.json', '')
        r['gene'] = DISPLAY_NAME.get(_g, _g)
        rows.append(r)
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)

    # classify structure first, THEN compute FDR over testable genes only
    def struct(r):
        if pd.isna(r['K']) or pd.isna(r['p_value']):
            return 'failed'
        if r['K'] < K_ZERO:
            return 'boundary'
        if r['K'] >= K_HIGH:
            return 'k_high'
        return 'testable'
    df['status'] = df.apply(struct, axis=1)

    testable = df.status == 'testable'
    df['n_tested'] = int(testable.sum())

    # BH q-values within this scheme, over testable genes only. Boundary/K-high
    # categories are excluded because an MLE on a parameter boundary invalidates
    # the chi-square LRT itself (Chernoff; Self & Liang 1987) -- an exclusion
    # about LRT validity, not multiple-testing. Categories and colours below
    # stay raw-p based; the dagger annotation carries the FDR tier.
    df['q_value'] = np.nan
    df.loc[testable, 'q_value'] = bh_fdr(df.loc[testable, 'p_value'].values)
    def interp(r):
        if r['status'] != 'testable':
            return r['status']
        if pd.notna(r['p_value']) and r['p_value'] < ALPHA:
            return 'relaxed' if r['K'] < 1 else 'intensified'
        return 'not_significant'
    df['interpretation'] = df.apply(interp, axis=1)

    df['sort_key'] = df['K'].fillna(999)
    return df.sort_values('sort_key').drop(columns='sort_key')


def stars(p):
    if pd.isna(p):
        return ''
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return ''


def k_value_plot(df, scheme, outdir):
    if df.empty:
        return
    XMAX, CAP = 5.25, 5.0

    order = pd.concat([
        df[df.interpretation == 'boundary'].sort_values('K'),
        df[df.interpretation.isin(['relaxed', 'not_significant', 'intensified'])].sort_values('K'),
        df[df.interpretation == 'k_high'].sort_values('K'),
        df[df.interpretation == 'failed'],
    ], ignore_index=True)

    n = len(order)
    fig, ax = plt.subplots(figsize=(10, max(4, n * 0.42)))
    labels = []

    for y, (_, r) in enumerate(order.iterrows()):
        it = r['interpretation']

        if it == 'failed':
            labels.append(f"{r['gene']} (no result)")
            ax.scatter(1.0, y, c=COLORS['failed'], s=90, marker='x', linewidths=2, zorder=3)
            continue

        if it == 'boundary':
            labels.append(f"{r['gene']} (K\u22480, excluded)")
            ax.scatter(0.06, y, s=130, marker='<', facecolors='none',
                       edgecolors=COLORS['boundary'], linewidths=1.8, zorder=3)
            ax.annotate('boundary \u2014 LRT invalid', (0.06, y), xytext=(0.28, y),
                        fontsize=6.5, color=COLORS['boundary'], va='center',
                        style='italic')
            continue

        if it == 'k_high':
            labels.append(f"{r['gene']} (K={r['K']:.1f}, excluded)")
            ax.scatter(CAP, y, c=COLORS['k_high'], s=100, marker='D',
                       edgecolors='black', linewidths=0.5, zorder=3)
            ax.annotate(f"K={r['K']:.1f}", (CAP, y), xytext=(CAP + 0.06, y),
                        fontsize=6.5, color=COLORS['k_high'], va='center',
                        arrowprops=dict(arrowstyle='->', color=COLORS['k_high'], lw=1.0))
            continue

        labels.append(r['gene'])
        ax.scatter(min(r['K'], CAP), y, c=COLORS[it], s=100, marker='o',
                   edgecolors='black', linewidths=0.5, zorder=3)
        s = stars(r['p_value'])
        if pd.notna(r.get('q_value')) and r['q_value'] < ALPHA:
            s += '\u2020'                       # dagger: survives BH at q<0.05
        if s:
            ax.annotate(s, (min(r['K'], CAP), y), xytext=(6, 0),
                        textcoords='offset points', fontsize=10,
                        fontweight='bold', va='center', color=COLORS[it])

    ax.axvline(1.0, color='black', ls='--', lw=1.4, alpha=0.7, zorder=1)
    ax.axvspan(0, 1, color='#3182CE', alpha=0.05, zorder=0)

    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlim(-0.1, XMAX)
    ax.set_ylim(-0.8, n - 0.2)
    ax.invert_yaxis()
    # Axis label in two pieces: the direction cue is anchored to K = 1 (data
    # coordinate) so its divider sits on the neutral line, and the parameter
    # name sits to its right. Offsets are in points, so they hold whatever the
    # figure height ends up being for this gene count.
    ax.set_xlabel('')
    ax.annotate('\u2190 relaxed | intensified \u2192',
                xy=(1.0, 0), xycoords=('data', 'axes fraction'),
                xytext=(8.5, -30), textcoords='offset points',
                ha='center', va='top', fontsize=11, annotation_clip=False)
    ax.annotate('K (relaxation parameter)',
                xy=(0.70, 0), xycoords='axes fraction',
                xytext=(0, -30), textcoords='offset points',
                ha='center', va='top', fontsize=11, annotation_clip=False)
    ntest = int(df.n_tested.iloc[0]) if 'n_tested' in df else 0
    ax.set_title(f"{SCHEME_NAMES.get(scheme, scheme)}\n"
                 f"per-gene LRT, {ntest} testable genes  "
                 f"(* p<0.05  ** p<0.01  *** p<0.001  \u2020 BH q<0.05)",
                 fontsize=11, fontweight='bold')
    ax.xaxis.grid(True, ls=':', alpha=0.4)
    ax.set_axisbelow(True)

    leg = [
        Line2D([], [], marker='o', ls='', markerfacecolor=COLORS['relaxed'],
               markeredgecolor='k', markersize=8, label='Relaxed (p<0.05, K<1)'),
        Line2D([], [], marker='o', ls='', markerfacecolor=COLORS['intensified'],
               markeredgecolor='k', markersize=8, label='Intensified (p<0.05, K>1)'),
        Line2D([], [], marker='o', ls='', markerfacecolor=COLORS['not_significant'],
               markeredgecolor='k', markersize=8, label='Not significant'),
        Line2D([], [], marker='<', ls='', markerfacecolor='none',
               markeredgecolor=COLORS['boundary'], markersize=9,
               label='K\u22480: MLE on boundary, excluded'),
        Line2D([], [], marker='D', ls='', markerfacecolor=COLORS['k_high'],
               markeredgecolor='k', markersize=8, label='K>5: unreliable, excluded'),
    ]
    ax.legend(handles=leg, loc='upper right', fontsize=7.5, framealpha=0.92)
    plt.tight_layout()
    base = os.path.join(outdir, f'{scheme}_K_values')
    for ext in ('pdf', 'png', 'svg'):
        fig.savefig(f'{base}.{ext}', dpi=600 if ext == 'png' else None, bbox_inches='tight')
    plt.close(fig)
    print(f'    {os.path.basename(base)}.(pdf|png|svg)')


def summary_counts(all_res, outdir):
    schemes = [s for s in PRIORITY if s in all_res and not all_res[s].empty]
    if not schemes:
        return
    cats = ['relaxed', 'intensified', 'not_significant', 'boundary', 'k_high']
    fig, ax = plt.subplots(figsize=(max(9, len(schemes) * 1.9), 5.4))
    x = np.arange(len(schemes))
    w = 0.16
    for i, c in enumerate(cats):
        vals = [int((all_res[s].interpretation == c).sum()) for s in schemes]
        b = ax.bar(x + (i - 2) * w, vals, w, label=c.replace('_', ' ').title(),
                   color=COLORS[c], edgecolor='black', linewidth=0.5)
        for bar, v in zip(b, vals):
            if v:
                ax.annotate(str(v), (bar.get_x() + bar.get_width() / 2, v),
                            xytext=(0, 2), textcoords='offset points',
                            ha='center', fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels([SCHEME_NAMES.get(s, s).split(' \u2014 ')[0] for s in schemes],
                       rotation=25, ha='right', fontsize=9)
    ax.set_ylabel('Genes', fontsize=11)
    ax.set_title('Focal 17-gene RELAX \u2014 per-gene LRT significance (p<0.05)', fontsize=12)
    ax.legend(fontsize=8)
    ax.yaxis.grid(True, ls=':', alpha=0.4); ax.set_axisbelow(True)
    plt.tight_layout()
    base = os.path.join(outdir, 'summary_counts')
    for ext in ('pdf', 'png', 'svg'):
        fig.savefig(f'{base}.{ext}', dpi=600 if ext == 'png' else None, bbox_inches='tight')
    plt.close(fig)
    print(f'    summary_counts.(pdf|png|svg)')


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT
    if not os.path.isdir(RELAX_DIR):
        sys.exit(f'ERROR: {RELAX_DIR} not found')
    os.makedirs(outdir, exist_ok=True)

    schemes = [d for d in sorted(os.listdir(RELAX_DIR))
               if os.path.isdir(os.path.join(RELAX_DIR, d))
               and glob.glob(os.path.join(RELAX_DIR, d, 'json', '*_RELAX.json'))]
    schemes = [s for s in PRIORITY if s in schemes] + [s for s in schemes if s not in PRIORITY]

    print('=' * 68)
    print('Focal 17-gene RELAX figures  (per-gene LRT; boundary-K excluded)')
    print('=' * 68)

    all_res = {}
    for s in schemes:
        df = load_scheme(s)
        all_res[s] = df
        if df.empty:
            continue
        c = df.interpretation.value_counts().to_dict()
        print(f"\n{s}   n={len(df)}  tested={int(df.n_tested.iloc[0])}")
        print(f"  relaxed={c.get('relaxed',0)}  intensified={c.get('intensified',0)}  "
              f"ns={c.get('not_significant',0)}  boundary={c.get('boundary',0)}  "
              f"K>5={c.get('k_high',0)}  failed={c.get('failed',0)}")
        nq = int(((df.q_value < ALPHA) & (df.status == 'testable')).sum())
        print(f"  FDR-significant (q<{ALPHA}): {nq}")
        rel = df[df.interpretation == 'relaxed']
        for _, r in rel.iterrows():
            dag = ' \u2020' if pd.notna(r['q_value']) and r['q_value'] < ALPHA else ''
            print(f"    RELAXED  {r['gene']:<10s} K={r['K']:.3f}  "
                  f"p={r['p_value']:.2e}  q={r['q_value']:.3f}{dag}")

    print('\nfigures:')
    for s in schemes:
        if not all_res[s].empty:
            k_value_plot(all_res[s], s, outdir)
    summary_counts(all_res, outdir)

    comb = []
    for s, df in all_res.items():
        if df.empty:
            continue
        d = df.copy()
        d['scheme'] = s
        d['scheme_name'] = SCHEME_NAMES.get(s, s)
        comb.append(d)
    if comb:
        out = os.path.join(outdir, 'all_17gene_relax_results.tsv')
        pd.concat(comb, ignore_index=True).to_csv(out, sep='\t', index=False)
        print(f'    all_17gene_relax_results.tsv')

    print(f'\n-> {outdir}/')
    print('\nNOTE: boundary-K genes are excluded from the BH denominator because an')
    print('MLE on the parameter boundary invalidates the chi-square LRT. They are')
    print('reported separately, not counted as relaxed.')


if __name__ == '__main__':
    main()
