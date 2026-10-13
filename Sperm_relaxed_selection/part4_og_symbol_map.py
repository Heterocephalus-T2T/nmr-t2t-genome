#!/usr/bin/env python3
"""Part 4: genome-wide orthogroup -> gene-symbol map, plus per-OG species
metrics. Run after OrthoFinder; every gene set selects from
these tables.

Symbols come from the human annotation only (human is the best-annotated
proteome): protein/transcript IDs are mapped to gene symbols via the GFF
(CDS -> transcript -> gene name), with the cleaned huma.fa headers as a
fallback. Every symbol found in an OG's human column is recorded - an OG
containing proteins of two genes gets two rows. 
We will match a gene list (including aliases) against this table in part 5d.

S

Inputs:  results/02c_orthofinder/Results_*/Orthogroups/Orthogroups.tsv
         annotations/Homo_sapiens.GRCh38.115.gff3
         results/02b_cleaned_proteins/huma.fa
Outputs (results/04_og_symbols/):
  og_symbol_map.tsv       Orthogroup, Symbol, Human_Protein_IDs, N_Human_Proteins
  og_species_metrics.tsv  Orthogroup, N_Species, Total_Proteins,
                          Max_Per_Species, Single_Copy_Pct, <one column per species>

Usage: python3 part4_og_symbol_map.py [--results_dir results/02c_orthofinder/Results_XXX]
       (default: the single/latest Results_* directory, printed loudly)
"""

import argparse
import glob
import os
import re
import sys

import pandas as pd

ORTHOFINDER_DIR = 'results/02c_orthofinder'
HUMAN_GFF = 'annotations/Homo_sapiens.GRCh38.115.gff3'
HUMA_FA = 'results/02b_cleaned_proteins/huma.fa'
OUTDIR = 'results/04_og_symbols'
HUMAN_COL = 'huma'


def normalize_protein(pid):
    if not isinstance(pid, str):
        return ""
    if pid.startswith('mikado') or re.match(r'^g\d+\.t\d+', pid):
        return pid
    return re.sub(r'\.\d+$', '', pid)


def parse_gff_attrs(attr):
    d = {}
    for m in re.finditer(r'(\w+)=([^;]+)', attr):
        d[m.group(1)] = m.group(2).strip('"')
    return d


def build_human_prot_to_symbol(gff):
    """Human protein & transcript IDs -> gene SYMBOL (uppercased)."""
    if not os.path.exists(gff):
        sys.exit(f"ERROR: human GFF not found ({gff})")
    gene_id_to_name, tx_to_name, prot_to_tx = {}, {}, {}
    with open(gff) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9:
                continue
            ftype, a = p[2], parse_gff_attrs(p[8])
            name = (a.get('gene_name') or a.get('Name') or a.get('gene') or '')
            if ftype == 'gene':
                gid = a.get('ID', '').replace('gene:', '')
                if gid and name:
                    gene_id_to_name[gid] = name
            elif ftype in ('mRNA', 'transcript'):
                tid = a.get('ID', '').replace('transcript:', '')
                parent = a.get('Parent', '').replace('gene:', '')
                if tid and name:
                    tx_to_name[tid] = name
                elif tid and parent in gene_id_to_name:
                    tx_to_name[tid] = gene_id_to_name[parent]
            elif ftype == 'CDS':
                prot = a.get('protein_id', '')
                parent = a.get('Parent', '').replace('transcript:', '')
                if prot and parent:
                    prot_to_tx[prot] = parent
    prot_to_symbol = {}
    for prot, tx in prot_to_tx.items():
        name = tx_to_name.get(tx)
        if name:
            for k in (prot, normalize_protein(prot)):
                prot_to_symbol[k] = name.upper()
    for tx, name in tx_to_name.items():
        for k in (tx, normalize_protein(tx)):
            prot_to_symbol.setdefault(k, name.upper())
    print(f"  human annotation: {len(prot_to_symbol)} protein/transcript -> symbol")
    return prot_to_symbol


def add_huma_header_symbols(fa, prot_to_symbol):
    """Fallback: gene_symbol=/gene= tags in the cleaned huma.fa headers."""
    if not os.path.exists(fa):
        return
    added = 0
    with open(fa) as f:
        for line in f:
            if not line.startswith('>'):
                continue
            h = line[1:].strip()
            pid = h.split()[0]
            m = re.search(r'gene_symbol[=:](\S+)', h) or re.search(r'gene[=:](\S+)', h)
            if m:
                name = m.group(1).rstrip(']')
                if len(name) < 30:
                    for k in (pid, normalize_protein(pid)):
                        if k not in prot_to_symbol:
                            prot_to_symbol[k] = name.upper()
                            added += 1
    print(f"  after huma.fa headers: +{added} mappings")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--results_dir', default=None,
                    help='OrthoFinder Results_* directory [latest]')
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)

    results_dir = args.results_dir
    if results_dir is None:
        matches = sorted(glob.glob(os.path.join(ORTHOFINDER_DIR, 'Results_*')))
        if not matches:
            sys.exit(f"ERROR: no Results_* under {ORTHOFINDER_DIR}")
        results_dir = matches[-1]
    print(f"Part 4: OG -> symbol map from {results_dir}")

    og_file = os.path.join(results_dir, 'Orthogroups', 'Orthogroups.tsv')
    if not os.path.exists(og_file):
        sys.exit(f"ERROR: {og_file} not found")

    prot_to_symbol = build_human_prot_to_symbol(HUMAN_GFF)
    add_huma_header_symbols(HUMA_FA, prot_to_symbol)

    og_df = pd.read_csv(og_file, sep='\t', low_memory=False)
    species_cols = [c for c in og_df.columns if c != 'Orthogroup']
    if HUMAN_COL not in species_cols:
        sys.exit(f"ERROR: human column '{HUMAN_COL}' not in OG table")
    print(f"  {len(og_df)} orthogroups, {len(species_cols)} species")

    symbol_rows, metric_rows = [], []
    n_no_human, n_unmapped = 0, 0

    for _, row in og_df.iterrows():
        og = row['Orthogroup']

        # ---- species metrics ----
        sp_counts = {}
        for sp in species_cols:
            v = str(row.get(sp, ''))
            if v in ('', 'nan'):
                continue
            n = len([x for x in v.split(',') if x.strip()])
            if n:
                sp_counts[sp] = n
        n_species = len(sp_counts)
        if n_species:
            single = sum(1 for c in sp_counts.values() if c == 1)
            m = {'Orthogroup': og, 'N_Species': n_species,
                 'Total_Proteins': sum(sp_counts.values()),
                 'Max_Per_Species': max(sp_counts.values()),
                 'Single_Copy_Pct': round(single / n_species * 100, 1)}
            for sp in species_cols:
                m[sp] = sp_counts.get(sp, 0)
            metric_rows.append(m)

        # ---- symbols from the human column ----
        hv = str(row.get(HUMAN_COL, ''))
        if hv in ('', 'nan'):
            n_no_human += 1
            continue
        sym_to_prots = {}
        any_mapped = False
        for pid in (x.strip() for x in hv.split(',') if x.strip()):
            sym = prot_to_symbol.get(pid) or prot_to_symbol.get(normalize_protein(pid))
            if sym:
                sym_to_prots.setdefault(sym, []).append(pid)
                any_mapped = True
        if not any_mapped:
            n_unmapped += 1
        for sym, pids in sorted(sym_to_prots.items()):
            symbol_rows.append({'Orthogroup': og, 'Symbol': sym,
                                'Human_Protein_IDs': ';'.join(pids),
                                'N_Human_Proteins': len(pids)})

    pd.DataFrame(symbol_rows).to_csv(
        os.path.join(OUTDIR, 'og_symbol_map.tsv'), sep='\t', index=False)
    pd.DataFrame(metric_rows).to_csv(
        os.path.join(OUTDIR, 'og_species_metrics.tsv'), sep='\t', index=False)

    n_ogs_with_symbol = len({r['Orthogroup'] for r in symbol_rows})
    n_symbols = len({r['Symbol'] for r in symbol_rows})
    print(f"\n  OGs with a symbol: {n_ogs_with_symbol}  "
          f"(no human member: {n_no_human}; human unmapped: {n_unmapped})")
    print(f"  distinct symbols: {n_symbols}")
    print(f"  -> {OUTDIR}/og_symbol_map.tsv, og_species_metrics.tsv")


if __name__ == '__main__':
    main()
