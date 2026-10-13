#!/usr/bin/env python3
"""
part11_fdr.py - Benjamini-Hochberg correction for the RELAX summaries.

Reads results/10_summary/{set}_summary.tsv (written by part10_summarise.py) and
adds a q_value column plus an FDR-based category, then reports per-scheme counts
under both the raw LRT p-value and the FDR threshold.

MULTIPLE-TESTING FAMILY
  BH is applied WITHIN each scheme, across the genes that were actually
  testable in that scheme. Boundary-value estimates (ARTIFACT_K0, K <= 0.01;
  ARTIFACT_KHIGH, K >= 5) and NO_DATA rows are excluded from the denominator,
  because an MLE on the parameter boundary invalidates the chi-square LRT
  (Chernoff; Self & Liang 1987) - that exclusion is about test validity, not
  about multiplicity. Excluded genes are reported separately, never counted as
  relaxed.

  --shared additionally repeats the correction on the intersection of genes
  testable in ALL of the single-species schemes, which is the like-for-like
  comparison for a cross-species gradient, and runs paired McNemar tests
  between schemes on that shared set.

Usage:
    python3 part11_fdr.py --set sperm
    python3 part11_fdr.py --set sperm --shared
    python3 part11_fdr.py --set targets17 --alpha 0.05
    python3 part11_fdr.py --all --shared
"""

import argparse
import itertools
import os

import numpy as np
import pandas as pd

SUMMARY_DIR = 'results/10_summary'
EXCLUDED = {'ARTIFACT_K0', 'ARTIFACT_KHIGH', 'NO_DATA'}
SINGLE_SPECIES = ['nmr_vs_nonmolerat', 'amr_vs_nonmolerat', 'bats_vs_nonmolerat']


def bh(p):
    """Benjamini-Hochberg q-values."""
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = p[order] * n / np.arange(1, n + 1)
    q[order] = np.minimum.accumulate(q[order][::-1])[::-1]
    return np.minimum(q, 1.0)


def mcnemar(a, b):
    """Paired McNemar with continuity correction on two boolean arrays."""
    from scipy import stats
    a = np.asarray(a, bool)
    b = np.asarray(b, bool)
    a_only = int((a & ~b).sum())
    b_only = int((~a & b).sum())
    n = a_only + b_only
    if n == 0:
        return a_only, b_only, 0.0, 1.0
    chi2 = (abs(a_only - b_only) - 1) ** 2 / n
    return a_only, b_only, chi2, float(stats.chi2.sf(chi2, 1))


def annotate(df, alpha):
    """Add q_value and Category_FDR, correcting within each scheme."""
    df = df.copy()
    df['q_value'] = np.nan
    for scheme, d in df.groupby('Scheme'):
        ok = d.index[~d['Category'].isin(EXCLUDED) & d['p_value'].notna()]
        if len(ok):
            df.loc[ok, 'q_value'] = bh(df.loc[ok, 'p_value'])

    def cat(r):
        if r['Category'] in EXCLUDED:
            return r['Category']
        if pd.isna(r['q_value']) or r['q_value'] >= alpha:
            return 'NOT_SIGNIFICANT'
        return 'RELAXED' if r['K'] < 1 else 'INTENSIFIED'

    df['Category_FDR'] = df.apply(cat, axis=1)
    return df


def report(df, gene_set, alpha):
    print('=' * 78)
    print(f'{gene_set}: BH within scheme, alpha = {alpha}')
    print('=' * 78)
    print(f"{'scheme':<28s} {'testable':>8s}  {'RELAXED':>16s}  {'INTENS':>14s}  "
          f"{'excluded':>8s}")
    print(f"{'':<28s} {'':>8s}  {'raw -> FDR':>16s}  {'raw -> FDR':>14s}")
    for scheme, d in df.groupby('Scheme'):
        testable = int((~d['Category'].isin(EXCLUDED)).sum())
        rel_p = int((d['Category'] == 'RELAXED').sum())
        rel_q = int((d['Category_FDR'] == 'RELAXED').sum())
        int_p = int((d['Category'] == 'INTENSIFIED').sum())
        int_q = int((d['Category_FDR'] == 'INTENSIFIED').sum())
        exc = int(d['Category'].isin(EXCLUDED).sum())
        pc_p = 100 * rel_p / testable if testable else 0
        pc_q = 100 * rel_q / testable if testable else 0
        print(f'{scheme:<28s} {testable:>8d}  {rel_p:>4d} -> {rel_q:<4d}'
              f'({pc_p:.1f}->{pc_q:.1f}%)  {int_p:>4d} -> {int_q:<4d}    {exc:>8d}')


def shared_analysis(df, alpha, gene_set):
    schemes = [s for s in SINGLE_SPECIES if s in set(df['Scheme'])]
    if len(schemes) < 2:
        print('\n(shared-gene analysis needs the single-species schemes)')
        return None
    testable = {}
    for s in schemes:
        d = df[(df['Scheme'] == s) & (~df['Category'].isin(EXCLUDED))]
        testable[s] = set(d['Gene'])
    shared = sorted(set.intersection(*testable.values()))
    n = len(shared)
    print('\n' + '=' * 78)
    print(f'{gene_set}: shared-gene analysis - {n} genes testable in all '
          f'{len(schemes)} single-species schemes')
    print('=' * 78)

    flags = {}
    rows = []
    for s in schemes:
        d = df[(df['Scheme'] == s) & (df['Gene'].isin(shared))].set_index('Gene').loc[shared]
        q = bh(d['p_value'].values)          # re-corrected on the shared set
        rel_q = (q < alpha) & (d['K'].values < 1)
        rel_p = (d['p_value'].values < alpha) & (d['K'].values < 1)
        flags[s] = rel_q
        rows.append({'scheme': s, 'n_shared': n,
                     'raw_relaxed': int(rel_p.sum()),
                     'raw_pct': round(100 * rel_p.sum() / n, 1),
                     'fdr_relaxed': int(rel_q.sum()),
                     'fdr_pct': round(100 * rel_q.sum() / n, 1)})
        print(f"  {s:<28s} raw {int(rel_p.sum()):>4d} ({100*rel_p.sum()/n:>5.1f}%)"
              f"   FDR {int(rel_q.sum()):>4d} ({100*rel_q.sum()/n:>5.1f}%)")

    print('\n  paired McNemar on the shared set (FDR calls):')
    mc = []
    for s1, s2 in itertools.combinations(schemes, 2):
        a, b, chi2, p = mcnemar(flags[s1], flags[s2])
        star = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        print(f"    {s1.split('_')[0]:>4s} vs {s2.split('_')[0]:<4s}  "
              f"{a:>4d} / {b:<4d}  chi2={chi2:7.3f}  p={p:.3g}  {star}")
        mc.append({'pair': f'{s1}_vs_{s2}', f'{s1}_only': a, f'{s2}_only': b,
                   'chi2': round(chi2, 3), 'p_value': p, 'sig': star})
    return pd.DataFrame(rows), pd.DataFrame(mc), shared


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--set', dest='gene_set', default=None,
                    help='sperm | targets17 | ddr | validation')
    ap.add_argument('--all', action='store_true', help='every summary found')
    ap.add_argument('--summary_dir', default=SUMMARY_DIR)
    ap.add_argument('--alpha', type=float, default=0.05)
    ap.add_argument('--shared', action='store_true',
                    help='also correct on genes testable in all single-species schemes')
    args = ap.parse_args()

    if args.all:
        sets = [f.replace('_summary.tsv', '')
                for f in sorted(os.listdir(args.summary_dir))
                if f.endswith('_summary.tsv')]
    elif args.gene_set:
        sets = [args.gene_set]
    else:
        raise SystemExit('give --set NAME or --all')

    for gs in sets:
        path = os.path.join(args.summary_dir, f'{gs}_summary.tsv')
        if not os.path.exists(path):
            print(f'{gs}: {path} not found -> skipped')
            continue
        df = pd.read_csv(path, sep='\t')
        df = annotate(df, args.alpha)
        report(df, gs, args.alpha)

        out = os.path.join(args.summary_dir, f'{gs}_summary_fdr.tsv')
        df.to_csv(out, sep='\t', index=False)
        print(f'\n  -> {out}')

        if args.shared:
            res = shared_analysis(df, args.alpha, gs)
            if res:
                nums, mc, shared = res
                nums.to_csv(os.path.join(args.summary_dir,
                                         f'{gs}_shared_fdr_numbers.tsv'),
                            sep='\t', index=False)
                mc.to_csv(os.path.join(args.summary_dir,
                                       f'{gs}_shared_mcnemar.tsv'),
                          sep='\t', index=False)
                with open(os.path.join(args.summary_dir,
                                       f'{gs}_shared_genes.txt'), 'w') as f:
                    f.write('\n'.join(shared) + '\n')
                print(f'  -> {gs}_shared_fdr_numbers.tsv, {gs}_shared_mcnemar.tsv, '
                      f'{gs}_shared_genes.txt')
        print()


if __name__ == '__main__':
    main()
