#!/usr/bin/env python3
"""Convert the DDR gene list (xlsx) into a SYMBOL-column TSV for the
orthogroup-discovery step.

Gene->OG mapping is already solved (human GFF + mygene synonyms + protein-membership
scanning), so the DDR genes are converted to the same input format and mapped
by the same code path rather than reimplementing it. The GENE and TIER
columns are located despite header decoration (e.g. "GENE \u25bc"); symbols
are deduplicated, first occurrence wins.

Orthogroup mapping is no longer a DDR-specific script: part 4 builds the
genome-wide OG->symbol map once, and part 5d selects this set from it with
the DDR admission filters (--min_species 15 --min_single_copy_pct 40, no
copy cap); multi-copy species are resolved later by ortholog selection
(part 9). This converter only prepares the gene list.

Usage:
  python3 part5c_ddr_gene_list.py \\
      --xlsx Balmus_DDR_genes_Apr2026.xlsx \\
      --out results/11_ddr_atlas/ddr_gene_list.tsv
  # then: python3 part5d_select_orthogroups.py \\
  #           --gene_list results/11_ddr_atlas/ddr_gene_list.tsv \\
  #           --outdir results/05_sets/ddr --min_species 15 --min_single_copy_pct 40
"""

import argparse
import os
import re
from collections import Counter

import pandas as pd


def clean_header(s):
    return re.sub(r'[^A-Za-z0-9]', '', str(s)).upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--xlsx', required=True)
    ap.add_argument('--out', default='results/11_ddr_atlas/ddr_gene_list.tsv')
    ap.add_argument('--gene_col', default='GENE')
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)

    df = pd.read_excel(args.xlsx)

    gene_c = tier_c = None
    for c in df.columns:
        cc = clean_header(c)
        if cc == clean_header(args.gene_col):
            gene_c = c
        if cc == 'TIER':
            tier_c = c
    if gene_c is None:
        gene_c = df.columns[0]
        print(f"  WARNING: no '{args.gene_col}' column; using first column '{gene_c}'")

    rows = []
    for _, r in df.iterrows():
        g = str(r[gene_c]).strip()
        if g and g.lower() != 'nan':
            rows.append({
                'SYMBOL': g,        # the column the mapping step reads
                'TIER': str(r[tier_c]).strip() if tier_c else '',
            })

    out = pd.DataFrame(rows).drop_duplicates('SYMBOL')
    out.to_csv(args.out, sep='\t', index=False)

    tc = Counter(out['TIER']) if 'TIER' in out else {}
    print(f"wrote {args.out}")
    print(f"  {len(out)} unique DDR gene symbols")
    if tc:
        print(f"  tiers: {dict(tc)}")
    print("\nNext: python3 part5d_select_orthogroups.py "
          f"--gene_list {args.out} --outdir results/05_sets/ddr "
          "--min_species 15 --min_single_copy_pct 40")


if __name__ == '__main__':
    main()
