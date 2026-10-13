#!/usr/bin/env python3
"""Part 5d: select orthogroups for a gene set, from part-4 tables.

Takes a gene list (SYMBOL column) and admission filters; expands symbols to
known synonyms via MyGene.info (cached; aliases mapping to more than one
gene are discarded as ambiguous); matches them case-insensitively against
the genome-wide OG->symbol map from part 4; applies the admission filters to
the part-4 species metrics; keeps the best orthogroup per gene (most
species); and writes the per-set selection plus the per-protein detail rows
that part 6 consumes.

One script, every set - the set is the arguments:
  validation: --gene_list gene_lists/validation.tsv --outdir results/05_sets/validation \\
              --min_species 4 --min_single_copy_pct 0
  sperm:      --gene_list results/10_extended_sperm_genes/go_sperm_genes.tsv \\
              --outdir results/05_sets/sperm --min_species 15 \\
              [PENDING D2: --min_single_copy_pct 70 --max_per_species 3 (original)
               vs --min_single_copy_pct 40 (unified, no cap)]
  ddr:        --gene_list results/11_ddr_atlas/ddr_gene_list.tsv \\
              --outdir results/05_sets/ddr --min_species 15 --min_single_copy_pct 40
The 17 targets are members of the sperm set (tagged); no separate selection.

Inputs:  results/04_og_symbols/og_symbol_map.tsv, og_species_metrics.tsv
         results/02c_orthofinder/Results_*/Orthogroups/Orthogroups.tsv
Outputs (--outdir):
  gene_og_list.tsv          Gene, Orthogroup, N_Species      (-> part 6)
  filtered_og_details.tsv   Gene, Orthogroup, Species_Col, Label, Protein_ID
  all_gene_orthogroups.tsv  every matched gene-OG pair with metrics (audit)
  selected_orthogroups.tsv  the kept best-OG rows (audit)
  genes_not_in_orthofinder.txt, gene_synonyms.tsv
"""

import argparse
import glob
import json
import os
import sys
import time

import pandas as pd

try:
    import mygene
except ImportError:
    sys.exit("ERROR: pip install mygene --break-system-packages")

OG_SYMBOLS_DIR = 'results/04_og_symbols'
ORTHOFINDER_DIR = 'results/02c_orthofinder'

# OrthoFinder column -> analysis label (all 23 species)
COL_TO_LABEL = {
    'hetg': 'NMR', 'fuka': 'AMR', 'bats': 'BATS',
    'huma': 'HUMA', 'bono': 'BONO', 'gori': 'GORI', 'maca': 'MACA',
    'mous': 'MOUS', 'nrat': 'RAT', 'nade': 'NADE', 'pvol': 'PVOL', 'wwmo': 'WWMO',
    'gpig': 'GPIG', 'chin': 'CHIN', 'degu': 'DEGU',
    'rabb': 'RABB', 'squi': 'SQUI', 'elep': 'ELEP', 'csmu': 'ACAH',
    'pigs': 'PIG', 'catt': 'BOS', 'cats': 'CAT', 'dogs': 'DOG',
}


def fetch_synonyms(symbols, outdir):
    """MyGene aliases, cached per set; UPPER synonym -> gene, ambiguous dropped."""
    cache = os.path.join(outdir, "mygene_synonym_cache.json")
    if os.path.exists(cache):
        print(f"  using cached synonyms ({cache})")
        raw = json.load(open(cache))
    else:
        print(f"  querying mygene for {len(symbols)} symbols...")
        mg = mygene.MyGeneInfo()
        raw = {}
        syms = sorted(symbols)
        for i in range(0, len(syms), 200):
            batch = syms[i:i + 200]
            print(f"    batch {i // 200 + 1}/{(len(syms) + 199) // 200}")
            try:
                res = mg.querymany(batch, scopes='symbol',
                                   fields='symbol,alias', species='human',
                                   returnall=True)
                for hit in res.get('out', []):
                    q = hit.get('query', '')
                    al = hit.get('alias', [])
                    if isinstance(al, str):
                        al = [al]
                    raw[q] = list({x for x in ([q, hit.get('symbol', q)] + al) if x})
            except Exception as e:
                print(f"    warn: {e}")
                for s in batch:
                    raw.setdefault(s, [s])
            time.sleep(0.4)
        json.dump(raw, open(cache, 'w'), indent=2)

    syn_to_gene = {}
    for gene, al in raw.items():
        for s in [gene] + list(al):
            k = s.upper()
            if k in syn_to_gene and syn_to_gene[k] != gene.upper():
                syn_to_gene[k] = None   # ambiguous: same alias, two genes
            else:
                syn_to_gene[k] = gene.upper()
    return raw, {k: v for k, v in syn_to_gene.items() if v}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gene_list', required=True, help='TSV with a SYMBOL column')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--min_species', type=int, default=15)
    ap.add_argument('--min_single_copy_pct', type=float, default=40)
    ap.add_argument('--max_per_species', type=int, default=None,
                    help='per-species copy cap; omit for no cap')
    ap.add_argument('--rescue_list', default=None,
                    help='TSV with a SYMBOL column (e.g. gene_lists/targets17.tsv): '
                         'these genes skip the two admission filters (curated '
                         'targets; multi-copy is resolved by ortholog selection)')
    ap.add_argument('--results_dir', default=None,
                    help='OrthoFinder Results_* directory [latest]')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Part 5d: orthogroup selection for {args.gene_list}")
    print(f"  filters: N_Species>={args.min_species}, "
          f"single-copy>={args.min_single_copy_pct}%, "
          f"max/species{'<=' + str(args.max_per_species) if args.max_per_species else ': no cap'}")

    go = pd.read_csv(args.gene_list, sep='\t')
    symbols = {s.upper() for s in go['SYMBOL'].dropna().unique()}
    print(f"  gene symbols: {len(symbols)}")

    raw_syn, syn_to_gene = fetch_synonyms(symbols, args.outdir)
    pd.DataFrame([{'Gene': g, 'Synonyms': ';'.join(sorted(set(a) - {g}))}
                  for g, a in sorted(raw_syn.items())]
                 ).to_csv(os.path.join(args.outdir, 'gene_synonyms.tsv'),
                          sep='\t', index=False)

    sym_map = pd.read_csv(os.path.join(OG_SYMBOLS_DIR, 'og_symbol_map.tsv'), sep='\t')
    metrics = pd.read_csv(os.path.join(OG_SYMBOLS_DIR, 'og_species_metrics.tsv'), sep='\t')

    # match: OG symbol (annotation) -> list gene, via the synonym map
    sym_map['Gene'] = sym_map['Symbol'].str.upper().map(syn_to_gene)
    matched = sym_map.dropna(subset=['Gene'])
    gene_ogs = (matched[['Gene', 'Orthogroup']].drop_duplicates()
                .merge(metrics, on='Orthogroup', how='left'))
    found = set(gene_ogs['Gene'])
    print(f"  matched {len(found)} genes across "
          f"{gene_ogs['Orthogroup'].nunique()} OGs")

    not_found = sorted(symbols - found)
    with open(os.path.join(args.outdir, 'genes_not_in_orthofinder.txt'), 'w') as f:
        f.write('\n'.join(not_found) + '\n')
    print(f"  not found: {len(not_found)}")

    audit_cols = ['Gene', 'Orthogroup', 'N_Species', 'Total_Proteins',
                  'Max_Per_Species', 'Single_Copy_Pct']
    gene_ogs[audit_cols].to_csv(
        os.path.join(args.outdir, 'all_gene_orthogroups.tsv'), sep='\t', index=False)

    rescue = set()
    if args.rescue_list:
        rescue = {x.upper() for x in
                  pd.read_csv(args.rescue_list, sep='\t')['SYMBOL'].dropna()}
    keep_rescued = gene_ogs[gene_ogs['Gene'].isin(rescue)]

    f = gene_ogs[gene_ogs['N_Species'] >= args.min_species]
    print(f"  N_Species>={args.min_species}: {len(f)}")
    f = f[f['Single_Copy_Pct'] >= args.min_single_copy_pct]
    print(f"  single-copy>={args.min_single_copy_pct}%: {len(f)}")
    if args.max_per_species is not None:
        f = f[f['Max_Per_Species'] <= args.max_per_species]
        print(f"  max<={args.max_per_species}/species: {len(f)}")
    if len(keep_rescued):
        n_before = len(f)
        f = pd.concat([f, keep_rescued]).drop_duplicates(['Gene', 'Orthogroup'])
        rescued_genes = sorted(set(keep_rescued['Gene']) - set(f.head(n_before)['Gene']))
        print(f"  rescued past filters (targets): "
              f"{', '.join(rescued_genes) if rescued_genes else 'none needed'}")
    f = (f.sort_values('N_Species', ascending=False)
          .drop_duplicates('Gene', keep='first'))
    print(f"  best OG per gene: {len(f)}")
    f[audit_cols].to_csv(
        os.path.join(args.outdir, 'selected_orthogroups.tsv'), sep='\t', index=False)
    f[['Gene', 'Orthogroup', 'N_Species']].to_csv(
        os.path.join(args.outdir, 'gene_og_list.tsv'), sep='\t', index=False)

    # per-protein detail rows for part 6, from the OG table (kept OGs only)
    results_dir = args.results_dir
    if results_dir is None:
        matches = sorted(glob.glob(os.path.join(ORTHOFINDER_DIR, 'Results_*')))
        if not matches:
            sys.exit(f"ERROR: no Results_* under {ORTHOFINDER_DIR}")
        results_dir = matches[-1]
    og_df = pd.read_csv(os.path.join(results_dir, 'Orthogroups', 'Orthogroups.tsv'),
                        sep='\t', low_memory=False)
    species_cols = [c for c in og_df.columns if c != 'Orthogroup']

    og_to_genes = {}
    for _, r in f.iterrows():
        og_to_genes.setdefault(r['Orthogroup'], []).append(r['Gene'])

    details = []
    for _, row in og_df.iterrows():
        genes = og_to_genes.get(row['Orthogroup'])
        if not genes:
            continue
        for gene in genes:
            for sp in species_cols:
                v = str(row.get(sp, ''))
                if v in ('', 'nan'):
                    continue
                for pid in (x.strip() for x in v.split(',') if x.strip()):
                    details.append({'Gene': gene, 'Orthogroup': row['Orthogroup'],
                                    'Species_Col': sp,
                                    'Label': COL_TO_LABEL.get(sp, sp),
                                    'Protein_ID': pid})
    pd.DataFrame(details).to_csv(
        os.path.join(args.outdir, 'filtered_og_details.tsv'), sep='\t', index=False)

    print(f"\n  gene_og_list.tsv + filtered_og_details.tsv -> {args.outdir}/ "
          f"(ready for part 6)")


if __name__ == '__main__':
    main()
