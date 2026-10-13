#!/usr/bin/env python3
"""Part 10: classify RELAX results and write the summary tables.

Reads every results/09_relax/{set}/{scheme}/relax_results.tsv, applies the
artifact thresholds and significance classification. 
Writes per-set tables plus one combined table.

Classification (alpha = 0.05):
  K <= 0.01              ARTIFACT_K0      (boundary; excluded from tests,
  K >= 5                 ARTIFACT_KHIGH    examined qualitatively - see the
                                           methods, Section 9.2)
  p >= 0.05              NOT_SIGNIFICANT
  p < 0.05 and K < 1     RELAXED
  p < 0.05 and K > 1     INTENSIFIED

Outputs: results/10_summary/{set}_summary.tsv       (per gene x scheme)
         results/10_summary/all_sets_combined.tsv
         per-scheme category counts on the console

Usage: python3 part10_summarise.py
"""

import glob
import os

import pandas as pd

IN_BASE = 'results/09_relax'
OUTDIR = 'results/10_summary'

K_ARTIFACT_LOW = 0.01
K_ARTIFACT_HIGH = 5
ALPHA = 0.05


def classify(K, p):
    if K is None or p is None or pd.isna(K) or pd.isna(p):
        return 'NO_DATA'
    if K <= K_ARTIFACT_LOW:
        return 'ARTIFACT_K0'
    if K >= K_ARTIFACT_HIGH:
        return 'ARTIFACT_KHIGH'
    if p >= ALPHA:
        return 'NOT_SIGNIFICANT'
    return 'RELAXED' if K < 1 else 'INTENSIFIED'


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    files = sorted(glob.glob(os.path.join(IN_BASE, '*', '*', 'relax_results.tsv')))
    if not files:
        raise SystemExit(f"ERROR: no relax_results.tsv under {IN_BASE} (run part 9)")

    frames = [pd.read_csv(f, sep='\t') for f in files]
    df = pd.concat(frames, ignore_index=True)
    df['Category'] = [classify(k, p) for k, p in zip(df['K'], df['p_value'])]

    print(f"Part 10: {len(df)} gene x scheme results from {len(files)} files "
          f"(artifact bounds: K<={K_ARTIFACT_LOW}, K>={K_ARTIFACT_HIGH})")

    for gene_set, sub in df.groupby('Set'):
        out = os.path.join(OUTDIR, f'{gene_set}_summary.tsv')
        sub.sort_values(['Scheme', 'p_value'], na_position='last').to_csv(
            out, sep='\t', index=False)
        print(f"\n{gene_set}  ->  {out}")
        for scheme, ss in sub.groupby('Scheme'):
            counts = ss['Category'].value_counts()
            line = '  '.join(f"{c}={counts.get(c, 0)}"
                             for c in ['RELAXED', 'INTENSIFIED', 'NOT_SIGNIFICANT',
                                       'ARTIFACT_K0', 'ARTIFACT_KHIGH', 'NO_DATA']
                             if counts.get(c, 0))
            print(f"  {scheme:28s} n={len(ss):4d}  {line}")

    df.to_csv(os.path.join(OUTDIR, 'all_sets_combined.tsv'), sep='\t', index=False)
    print(f"\ncombined -> {OUTDIR}/all_sets_combined.tsv")


if __name__ == '__main__':
    main()
