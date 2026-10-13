#!/usr/bin/env python3
"""
fig_gradient_intersection_5scheme.py

The gradient figure with all five schemes, reading the merged results table
(relax_all_schemes_merged.tsv). Five bars:

  single-species (the sociality gradient, McNemar brackets among them):
    HetGla (eusocial)  |  FukAns (cooperative breeder)  |  BatSui (solitary)
  --- divider ---
  combined lineages (more test branches -> more power, shown separately):
    Social (HetGla + FukAns)  |  All mole-rats (Bathyergidae)

Design decisions:
  - Gene set = the 3-WAY intersection of the single-species schemes (genes tested
    with a valid, non-boundary result in nmr AND amr AND bats). This keeps the
    single-species numbers stable (they don't change when the combined schemes are
    added) and puts every bar on the SAME gene set so they're comparable.
  - The social/allmolerat bars are computed on those same intersection genes.
  - FDR (BH q<0.05, within scheme over the intersection genes) is the significance
    rule -- this is the cross-species comparison, so FDR applies.
  - McNemar brackets are drawn ONLY among the three single species (that is the
    paired sociality-gradient test). The combined-lineage bars are shown without
    brackets because they have more test branches (more power) and are not part of
    the same-power gradient comparison.

Usage:
  python3 fig_gradient_intersection_5scheme.py \
      --merged results/10_extended_sperm_genes/relax_all_schemes_merged.tsv \
      --outdir results/10_extended_sperm_genes/figures_gradient
"""

import os, argparse, itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

# scheme keys in the merged TSV
SINGLE = ['nmr_vs_nonmolerat', 'amr_vs_nonmolerat', 'bats_vs_nonmolerat']
COMBINED = ['social_vs_nonmolerat', 'allmolerat_vs_nonmolerat']
ALL_SCHEMES = SINGLE + COMBINED

LABELS = {
    'nmr_vs_nonmolerat':        'HetGla\n(eusocial)',
    'amr_vs_nonmolerat':        'FukAns\n(coop. breeder)',
    'bats_vs_nonmolerat':       'BatSui\n(solitary)',
    'social_vs_nonmolerat':     'Eusocial+coop. breeder\n(HetGla+FukAns)',
    'allmolerat_vs_nonmolerat': 'All mole-rats\n(Bathyergidae)',
}
COLORS = {
    'nmr_vs_nonmolerat':        '#E63946',
    'amr_vs_nonmolerat':        '#457B9D',
    'bats_vs_nonmolerat':       '#2A9D8F',
    'social_vs_nonmolerat':     '#6D597A',
    'allmolerat_vs_nonmolerat': '#8D6E63',
}
ALPHA = 0.05


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p)
    if n == 0: return p
    o = np.argsort(p); q = np.empty(n)
    q[o] = p[o] * n / np.arange(1, n + 1)
    q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.minimum(q, 1.0)


def mcnemar(a, b):
    a = np.asarray(a, bool); b = np.asarray(b, bool)
    a_only = int((a & ~b).sum()); b_only = int((~a & b).sum())
    n = a_only + b_only
    if n == 0: return a_only, b_only, 0.0, 1.0
    chi2 = (abs(a_only - b_only) - 1) ** 2 / n
    return a_only, b_only, chi2, float(stats.chi2.sf(chi2, 1))


def stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--merged',
                    default='results/10_extended_sperm_genes/relax_all_schemes_merged.tsv')
    ap.add_argument('--outdir',
                    default='results/10_extended_sperm_genes/figures_gradient')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if not os.path.exists(args.merged):
        raise SystemExit(f"merged table not found: {args.merged}")
    df = pd.read_csv(args.merged, sep='\t')
    df['K'] = pd.to_numeric(df['K'], errors='coerce')
    df['p_value'] = pd.to_numeric(df['p_value'], errors='coerce')

    # per-scheme frames indexed by gene
    by_scheme = {}
    for s in ALL_SCHEMES:
        d = df[df['Scheme'] == s].set_index('Gene')
        by_scheme[s] = d
        print(f"  {s}: {len(d)} genes")

    # ---- 3-way intersection of the single-species schemes ----
    common = set(by_scheme[SINGLE[0]].index)
    for s in SINGLE[1:]:
        common &= set(by_scheme[s].index)
    # drop genes that are boundary-K / no-data in ANY single scheme
    bad = set()
    for s in SINGLE:
        d = by_scheme[s]
        for g in list(common):
            cat = d.loc[g, 'Category'] if g in d.index else 'NO_DATA'
            if cat in ('ARTIFACT_K0', 'ARTIFACT_KHIGH', 'NO_DATA'):
                bad.add(g)
    genes = sorted(common - bad)
    n = len(genes)
    print(f"\n  3-way intersection (single species): {len(common)} -> "
          f"{n} after removing boundary-K")

    # ---- relaxed flags per scheme on the intersection genes (FDR) ----
    def flags_for(scheme):
        d = by_scheme[scheme]
        # some intersection genes may be missing / boundary in a COMBINED scheme;
        # treat those as "not relaxed" (can't call it) but keep the gene in n.
        p = np.array([d.loc[g, 'p_value'] if g in d.index else np.nan for g in genes])
        K = np.array([d.loc[g, 'K'] if g in d.index else np.nan for g in genes])
        cat = np.array([d.loc[g, 'Category'] if g in d.index else 'NO_DATA' for g in genes])
        # FDR only over genes with a valid (non-boundary) p in THIS scheme
        valid = np.array([c not in ('ARTIFACT_K0', 'ARTIFACT_KHIGH', 'NO_DATA')
                          and np.isfinite(pp) for c, pp in zip(cat, p)])
        q = np.ones(len(genes))
        if valid.sum() > 0:
            q[valid] = bh_fdr(p[valid])
        flag = valid & (q < ALPHA) & (K < 1)
        return flag, int(valid.sum())

    flags = {}
    valid_counts = {}
    raw_flags = {}
    for s in ALL_SCHEMES:
        flags[s], valid_counts[s] = flags_for(s)

    # raw-p flags (for the raw-vs-FDR panel) on the same intersection genes
    def raw_flags_for(scheme):
        d = by_scheme[scheme]
        p = np.array([d.loc[g, 'p_value'] if g in d.index else np.nan for g in genes])
        K = np.array([d.loc[g, 'K'] if g in d.index else np.nan for g in genes])
        cat = np.array([d.loc[g, 'Category'] if g in d.index else 'NO_DATA' for g in genes])
        valid = np.array([c not in ('ARTIFACT_K0', 'ARTIFACT_KHIGH', 'NO_DATA')
                          and np.isfinite(pp) for c, pp in zip(cat, p)])
        return valid & (p < ALPHA) & (K < 1)
    for s in ALL_SCHEMES:
        raw_flags[s] = raw_flags_for(s)

    # ---- figure: two panels (A = FDR gradient + McNemar, B = raw vs FDR) ----
    fig, axes = plt.subplots(1, 2, figsize=(18, 6.6),
                             gridspec_kw={'width_ratios': [1.15, 1.0]})
    ax = axes[0]
    order = SINGLE + COMBINED
    xpos = [0, 1, 2, 3.5, 4.5]

    pct = [100 * flags[s].sum() / n for s in order]
    cnt = [int(flags[s].sum()) for s in order]
    cols = [COLORS[s] for s in order]

    bars = ax.bar(xpos, pct, color=cols, edgecolor='black', linewidth=1.1,
                  width=0.72, zorder=3)
    for i, s in enumerate(order):
        if s in COMBINED:
            bars[i].set_hatch('///'); bars[i].set_alpha(0.85)
    for xi, p_, c_ in zip(xpos, pct, cnt):
        ax.annotate(f'{p_:.1f}%\n{c_}/{n}', (xi, p_), xytext=(0, 3),
                    textcoords='offset points', ha='center', va='bottom',
                    fontsize=10, fontweight='bold')
    ax.set_xticks(xpos)
    ax.set_xticklabels([LABELS[s] for s in order], fontsize=9)
    ax.set_ylabel('Sperm genes relaxed (%)', fontsize=12)
    ax.set_title(f'A) Relaxed genes across speceies/lineages (FDR q<0.05, n={n} shared genes)',
                 fontsize=12, fontweight='bold', pad=12)
    ax.set_ylim(0, max(pct) * 1.55)
    ax.yaxis.grid(True, ls=':', alpha=0.4); ax.set_axisbelow(True)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.axvline(3.0 - 0.5, color='grey', ls='--', lw=1.0, alpha=0.7)
    ax.text(1.0, max(pct) * 1.5, 'single species',
            ha='center', fontsize=8.5, style='italic', color='#555')
    ax.text(4.0, max(pct) * 1.5, 'combined lineages\n(multple test branches)',
            ha='center', fontsize=8.5, style='italic', color='#555')

    pairs = [(0, 1, SINGLE[0], SINGLE[1]),
             (1, 2, SINGLE[1], SINGLE[2]),
             (0, 2, SINGLE[0], SINGLE[2])]
    single_pct = pct[:3]
    lvl = max(single_pct) * 1.10; step = max(single_pct) * 0.13
    for i, (xi, xj, s1, s2) in enumerate(pairs):
        _, _, _, pv = mcnemar(flags[s1], flags[s2])
        y = lvl + i * step
        ax.plot([xpos[xi], xpos[xi], xpos[xj], xpos[xj]],
                [y, y + step * 0.22, y + step * 0.22, y], color='black', lw=1.0)
        ax.text((xpos[xi] + xpos[xj]) / 2, y + step * 0.24,
                f'{stars(pv)} p={pv:.1e}', ha='center', va='bottom', fontsize=8)

    # ---- panel B: raw LRT vs FDR, all five schemes ----
    axB = axes[1]
    w = 0.38
    xb = np.array([0, 1, 2, 3.5, 4.5])
    raw_pct = [100 * raw_flags[s].sum() / n for s in order]
    fdr_pct = [100 * flags[s].sum() / n for s in order]
    axB.bar(xb - w/2, raw_pct, w, label='raw LRT (p<0.05)', color='#b8b8b8',
            edgecolor='black', zorder=3)
    b2 = axB.bar(xb + w/2, fdr_pct, w, label='FDR (q<0.05)', color=cols,
                 edgecolor='black', zorder=3)
    for i, s in enumerate(order):
        if s in COMBINED:
            b2[i].set_hatch('///'); b2[i].set_alpha(0.85)
    for xi, r_, f_ in zip(xb, raw_pct, fdr_pct):
        axB.annotate(f'{r_:.1f}', (xi - w/2, r_), xytext=(0, 2),
                     textcoords='offset points', ha='center', fontsize=8)
        axB.annotate(f'{f_:.1f}', (xi + w/2, f_), xytext=(0, 2),
                     textcoords='offset points', ha='center', fontsize=8, fontweight='bold')
    axB.set_xticks(xb)
    axB.set_xticklabels([LABELS[s] for s in order], fontsize=9)
    axB.set_ylabel('Sperm genes relaxed (%)', fontsize=12)
    axB.set_title(f'B) Raw LRT vs FDR (same {n} shared genes)',
                  fontsize=12, fontweight='bold', pad=12)
    axB.legend(frameon=False, fontsize=9, loc='upper left')
    axB.set_ylim(0, max(raw_pct) * 1.28)
    axB.yaxis.grid(True, ls=':', alpha=0.4); axB.set_axisbelow(True)
    for sp in ('top', 'right'):
        axB.spines[sp].set_visible(False)
    axB.axvline(3.0 - 0.5, color='grey', ls='--', lw=1.0, alpha=0.7)

    plt.tight_layout()
    out = os.path.join(args.outdir, 'atlas_gradient_intersection_5scheme')
    plt.savefig(out + '.png', dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig(out + '.pdf', bbox_inches='tight', facecolor='white')
    print(f"\n  wrote {out}.png / .pdf")

    # ---- save numbers ----
    rows = []
    for s in order:
        rows.append({'scheme': s, 'label': LABELS[s].replace('\n', ' '),
                     'n_shared': n,
                     'raw_relaxed': int(raw_flags[s].sum()),
                     'raw_pct': round(100 * raw_flags[s].sum() / n, 1),
                     'fdr_relaxed': int(flags[s].sum()),
                     'fdr_pct': round(100 * flags[s].sum() / n, 1),
                     'valid_tested': valid_counts[s],
                     'group': 'single' if s in SINGLE else 'combined'})
    pd.DataFrame(rows).to_csv(
        os.path.join(args.outdir, 'gradient_5scheme_numbers.tsv'), sep='\t', index=False)

    mc = []
    for s1, s2 in itertools.combinations(SINGLE, 2):
        a, b, chi2, pv = mcnemar(flags[s1], flags[s2])
        mc.append({'pair': f'{s1}_vs_{s2}', f'{s1}_only': a, f'{s2}_only': b,
                   'chi2': round(chi2, 3), 'p_value': pv, 'sig': stars(pv)})
    pd.DataFrame(mc).to_csv(
        os.path.join(args.outdir, 'gradient_5scheme_mcnemar.tsv'), sep='\t', index=False)

    # ---- console summary ----
    print("\n  gradient (FDR, intersection n={}):".format(n))
    for s in order:
        tag = '' if s in SINGLE else '  [combined: more power]'
        print(f"    {LABELS[s].replace(chr(10),' '):26s}: "
              f"{int(flags[s].sum())}/{n} = {100*flags[s].sum()/n:.1f}%{tag}")
    print("\n  McNemar (single-species gradient):")
    for s1, s2 in itertools.combinations(SINGLE, 2):
        _, _, _, pv = mcnemar(flags[s1], flags[s2])
        n1 = s1.split('_')[0]; n2 = s2.split('_')[0]
        print(f"    {n1} vs {n2}: p={pv:.2e} {stars(pv)}")


if __name__ == '__main__':
    main()
