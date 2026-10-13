#!/usr/bin/env python3
"""Part 6: assemble per-gene CDS sets with focal replacement, for one gene set.

Reads the part-5d selection (gene_og_list.tsv + filtered_og_details.tsv) and
writes one FASTA per gene. Reference species keep ALL resolvable copies,
labelled SPECIES, SPECIES_c2, ... (ortholog selection happens in part 9).
HetGla/BatSui are handled by exclusive replacement: 
these two species are not in the FASTA; and most one curated sequence per gene is added.

Focal fallback policy when a gene has no curated-TSV symbol match
(--focal_fallback):
  keep_og   (default) the first focal orthogroup member's CDS is used,
            action 'keep_og' - the behaviour of every previous run
  synonyms  the TSV symbol match is extended with the MyGene synonym cache
            written by part 5d, converting alias-named genes into proper
            'replace'/'add'; genes still unmatched fall back to keep_og
  strict    curated-only: unmatched focal genes stay absent

Every focal decision lands in focal_replacement_log.tsv with a length-ratio
flag (SHORT_FRAGMENT <0.5x / LONG_OUTLIER >1.5x the orthogroup median, the
median taken over all retained reference copies). Genes with fewer than 4
sequences are dropped.


Inputs:  results/05_sets/{set}/gene_og_list.tsv, filtered_og_details.tsv,
         mygene_synonym_cache.json (for --focal_fallback synonyms)
         focal TSVs + focal CDS + cds/ per-species files (constants below)
Outputs: results/06_cds/{set}/cds_per_gene/{GENE}.fa
         results/06_cds/{set}/focal_replacement_log.tsv, cds_assembly_stats.tsv

Usage: python3 part6_assemble_cds.py --set sperm [--focal_fallback keep_og]
"""

import argparse
import glob
import json
import os
import re
import statistics
import sys

import pandas as pd

SETS_DIR = 'results/05_sets'
OUT_BASE = 'results/06_cds'

HETG_TSV = 'mHetGlaPri.single.copy.tsv'
BATS_TSV = 'BatSui.single.copy.tsv'
HETG_CDS = 'cds/mHetGla1.1.pri.cleaned.cds.fa'
BATS_CDS = 'cds/BatSuiPri.cds.fa'

COL_TO_LABEL = {
    'hetg': 'NMR', 'fuka': 'AMR', 'bats': 'BATS',
    'huma': 'HUMA', 'bono': 'BONO', 'gori': 'GORI', 'maca': 'MACA',
    'mous': 'MOUS', 'nrat': 'RAT', 'nade': 'NADE', 'pvol': 'PVOL', 'wwmo': 'WWMO',
    'gpig': 'GPIG', 'chin': 'CHIN', 'degu': 'DEGU',
    'rabb': 'RABB', 'squi': 'SQUI', 'elep': 'ELEP', 'csmu': 'ACAH',
    'pigs': 'PIG', 'catt': 'BOS', 'cats': 'CAT', 'dogs': 'DOG',
}
FOCAL = {'hetg': 'NMR', 'bats': 'BATS'}

CDS_DIR = 'cds'
CDS_HINTS = {
    'fuka': ['DKFZ_Fans', 'Fans'],
    'huma': ['Homo_sapiens'],
    'mous': ['Mus_musculus'],
    'nrat': ['Rattus'],
    'pvol': ['Microtus'],
    'wwmo': ['Mus_spretus'],
    'gpig': ['Cavia'],
    'chin': ['Chinchilla'],
    'degu': ['Octodon'],
    'squi': ['Sciurus'],
    'catt': ['Bos_taurus'],
    'cats': ['Felis'],
    'dogs': ['Canis'],
    'bono': ['bonobo'],
    'gori': ['gorilla'],
    'maca': ['macaque'],
    'elep': ['elephant'],
    'rabb': ['rabbit'],
    'pigs': ['pig'],
    'csmu': ['csmus'],
    'nade': ['nadeermouse'],
}

CLEANED_PROT_DIR = 'results/02b_cleaned_proteins'


def resolve_cds_file(col):
    hints = CDS_HINTS.get(col, [])
    cands = sorted(glob.glob(os.path.join(CDS_DIR, '*')))

    def cds_like(fn):
        b = os.path.basename(fn).lower()
        return ('cds' in b) and b.endswith(('.fa', '.fna', '.faa', '.fasta'))

    for fn in cands:
        b = os.path.basename(fn).lower()
        if cds_like(fn) and any(h.lower() in b for h in hints):
            return fn
    for fn in cands:
        b = os.path.basename(fn).lower()
        if b.endswith(('.fa', '.fna', '.faa', '.fasta')) and any(h.lower() in b for h in hints):
            return fn
    return None


def read_fasta(path):
    seqs, cur, buf = {}, None, []
    if not os.path.exists(path):
        return seqs
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if cur is not None:
                    seqs[cur] = ''.join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            seqs[cur] = ''.join(buf)
    return seqs


def norm_id(pid):
    if not isinstance(pid, str):
        return ""
    if pid.startswith('mikado') or re.match(r'^g\d+\.t\d+', pid):
        return pid
    return re.sub(r'\.\d+$', '', pid)


def build_norm_index(seqs):
    idx = {}
    for k in seqs:
        idx.setdefault(norm_id(k), k)
    return idx


def fasta_lookup(seqs, nidx, pid):
    if pid in seqs:
        return seqs[pid]
    n = norm_id(pid)
    if n in seqs:
        return seqs[n]
    k = nidx.get(n)
    return seqs[k] if k is not None else None


def build_ensembl_prot_to_tx(cleaned_fa):
    m = {}
    if not os.path.exists(cleaned_fa):
        return m
    with open(cleaned_fa) as f:
        for line in f:
            if not line.startswith('>'):
                continue
            h = line[1:].strip()
            pid = h.split()[0]
            tx = re.search(r'transcript:(\S+)', h)
            if tx:
                t = tx.group(1)
                m[pid] = t
                m[norm_id(pid)] = norm_id(t)
    return m


def build_ncbi_protid_to_cdskey(cds_seqs):
    m = {}
    for key in cds_seqs:
        mt = re.search(r'_cds_([A-Z]{2}_[0-9]+\.[0-9]+)_', key)
        if mt:
            pid = mt.group(1)
            m[pid] = key
            m[norm_id(pid)] = key
    return m


def make_cds_resolver(col, cds_seqs):
    nidx = build_norm_index(cds_seqs)
    ens_p2t = build_ensembl_prot_to_tx(os.path.join(CLEANED_PROT_DIR, f'{col}.fa'))
    ncbi_p2k = build_ncbi_protid_to_cdskey(cds_seqs)

    def resolve(pid):
        s = fasta_lookup(cds_seqs, nidx, pid)
        if s:
            return s
        tx = ens_p2t.get(pid) or ens_p2t.get(norm_id(pid))
        if tx:
            s = fasta_lookup(cds_seqs, nidx, tx)
            if s:
                return s
        key = ncbi_p2k.get(pid) or ncbi_p2k.get(norm_id(pid))
        if key and key in cds_seqs:
            return cds_seqs[key]
        guess = re.sub(r'(ENS\w*?)P(\d)', r'\1T\2', pid)
        if guess != pid:
            s = fasta_lookup(cds_seqs, nidx, guess)
            if s:
                return s
        return None
    return resolve


def write_fasta(path, recs, wrap=60):
    with open(path, 'w') as f:
        for sid, seq in recs:
            f.write(f">{sid}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i:i + wrap] + "\n")


def load_focal_tsv(path, unique_col, hier_col, tid_col):
    df = pd.read_csv(path, sep='\t', dtype=str, low_memory=False)
    df.columns = [c.strip().strip('"') for c in df.columns]
    sym_to_tid = {}
    for _, row in df.iterrows():
        tid = row.get(tid_col)
        if not isinstance(tid, str) or not tid:
            continue
        for col in (unique_col, hier_col):
            val = row.get(col)
            if isinstance(val, str) and val and val.upper() != 'NA':
                for sym in re.split(r'[;,]', val):
                    sym = sym.strip()
                    if sym:
                        sym_to_tid.setdefault(sym.upper(), tid)
    return sym_to_tid


def extend_with_synonyms(sym_to_tid, cache_path):
    """--focal_fallback synonyms: alias -> primary via the 5d cache, so an
    alias-named TSV entry matches the list gene (and vice versa)."""
    if not os.path.exists(cache_path):
        print(f"  WARNING: {cache_path} not found; synonyms mode adds nothing")
        return sym_to_tid
    raw = json.load(open(cache_path))
    extended = dict(sym_to_tid)
    for gene, aliases in raw.items():
        g = gene.upper()
        names = {g} | {a.upper() for a in aliases}
        tids = [sym_to_tid[n] for n in names if n in sym_to_tid]
        if tids:
            for n in names:
                extended.setdefault(n, tids[0])
    return extended


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--set', required=True, help='validation | sperm | ddr')
    ap.add_argument('--focal_fallback', default='keep_og',
                    choices=['keep_og', 'synonyms', 'strict'],
                    help='policy when a focal gene has no curated match [keep_og]')
    args = ap.parse_args()

    set_dir = os.path.join(SETS_DIR, args.set)
    outdir = os.path.join(OUT_BASE, args.set)
    os.makedirs(os.path.join(outdir, 'cds_per_gene'), exist_ok=True)

    print(f"Part 6: CDS assembly for set '{args.set}' "
          f"(focal fallback: {args.focal_fallback})")

    for p in (HETG_CDS, BATS_CDS, HETG_TSV, BATS_TSV):
        if not os.path.exists(p):
            sys.exit(f"ERROR: {p} not found - fix the path constants at the "
                     f"top of this script for this folder's layout")

    gene_og = pd.read_csv(os.path.join(set_dir, 'gene_og_list.tsv'), sep='\t')
    details = pd.read_csv(os.path.join(set_dir, 'filtered_og_details.tsv'), sep='\t')
    genes = gene_og['Gene'].tolist()
    print(f"\nGenes: {len(genes)}")

    print("Loading focal single-copy TSVs...")
    hetg_sym2tid = load_focal_tsv(HETG_TSV, 'unique_gene', 'hierarchical_gene', 'tid')
    bats_sym2tid = load_focal_tsv(BATS_TSV, 'unique_gene', 'hierarchical_gene', 'tid')
    if args.focal_fallback == 'synonyms':
        cache = os.path.join(set_dir, 'mygene_synonym_cache.json')
        hetg_sym2tid = extend_with_synonyms(hetg_sym2tid, cache)
        bats_sym2tid = extend_with_synonyms(bats_sym2tid, cache)
    print(f"  HetGla TSV: {len(hetg_sym2tid)} symbol->tid")
    print(f"  BatSui TSV: {len(bats_sym2tid)} symbol->tid")
    focal_sym2tid = {'hetg': hetg_sym2tid, 'bats': bats_sym2tid}

    print("Loading focal CDS FASTAs...")
    focal_cds = {'hetg': read_fasta(HETG_CDS), 'bats': read_fasta(BATS_CDS)}
    focal_nidx = {c: build_norm_index(s) for c, s in focal_cds.items()}
    print(f"  HetGla CDS: {len(focal_cds['hetg'])} seqs")
    print(f"  BatSui CDS: {len(focal_cds['bats'])} seqs")

    print("Resolving + loading non-focal CDS files...")
    nonfocal_resolver = {}
    for col in COL_TO_LABEL:
        if col in FOCAL:
            continue
        path = resolve_cds_file(col)
        if path and os.path.exists(path):
            seqs = read_fasta(path)
            nonfocal_resolver[col] = make_cds_resolver(col, seqs)
            print(f"  {col}: {len(seqs)} seqs ({os.path.basename(path)})")
        else:
            print(f"  {col}: NO FILE FOUND -> will be absent")
            nonfocal_resolver[col] = (lambda pid: None)

    log_rows, stats_rows = [], []

    for gene in genes:
        gd = details[details['Gene'] == gene]
        og = gene_og[gene_og['Gene'] == gene]['Orthogroup'].iloc[0]

        col_to_prots = {}
        for _, row in gd.iterrows():
            col_to_prots.setdefault(row['Species_Col'], []).append(row['Protein_ID'])

        cds_recs, og_lengths = [], []

        # reference species: ALL resolvable copies, labelled _c2, _c3, ...
        for col, prots in col_to_prots.items():
            if col in FOCAL:
                continue
            label = COL_TO_LABEL.get(col, col)
            resolve = nonfocal_resolver.get(col, (lambda pid: None))
            n_kept = 0
            for pid in prots:
                s = resolve(pid)
                if s:
                    n_kept += 1
                    uid = label if n_kept == 1 else f"{label}_c{n_kept}"
                    cds_recs.append((uid, s))
                    og_lengths.append(len(s))

        og_median = statistics.median(og_lengths) if og_lengths else 0

        # focal species: exclusive replacement with full logging
        for col, label in FOCAL.items():
            present_in_og = col in col_to_prots
            og_prot = col_to_prots.get(col, [None])[0] if present_in_og else None
            tid = focal_sym2tid[col].get(gene.upper())

            action, tsv_tid, tsv_sym, cds_len, seq = None, 'none', 'none', 0, None

            if tid:
                seq = fasta_lookup(focal_cds[col], focal_nidx[col], tid)
                if seq:
                    tsv_tid, tsv_sym, cds_len = tid, gene.upper(), len(seq)
                    action = 'replace' if present_in_og else 'add'
                else:
                    action = 'tsv_tid_no_cds'
            if action is None or seq is None:
                if present_in_og and args.focal_fallback != 'strict':
                    s = fasta_lookup(focal_cds[col], focal_nidx[col], og_prot)
                    if s:
                        seq, cds_len = s, len(s)
                        action = action or 'keep_og'
                    else:
                        action = action or 'og_member_no_cds'
                elif present_in_og:
                    action = action or 'strict_absent'
                else:
                    action = action or 'no_tsv_absent'

            ratio = (cds_len / og_median) if (og_median and cds_len) else 0.0
            if cds_len == 0:
                flag = 'NO_SEQ'
            elif ratio == 0:
                flag = 'NO_OG_MEDIAN'
            elif ratio < 0.5:
                flag = 'SHORT_FRAGMENT'
            elif ratio > 1.5:
                flag = 'LONG_OUTLIER'
            else:
                flag = 'OK'

            if seq:
                cds_recs.append((label, seq))

            log_rows.append({
                'Gene': gene, 'Orthogroup': og, 'Species': label,
                'Action': action,
                'OG_member_before': og_prot if present_in_og else 'absent',
                'TSV_transcript': tsv_tid,
                'TSV_matched_symbol': tsv_sym,
                'CDS_length': cds_len,
                'OG_median_length': round(og_median, 1),
                'Length_ratio': round(ratio, 3),
                'Flag': flag,
            })

        n_species = len({sid.split('_c')[0] for sid, _ in cds_recs})
        if n_species >= 4:
            write_fasta(os.path.join(outdir, 'cds_per_gene', f"{gene}.fa"), cds_recs)
            status = 'ready'
        else:
            status = f'too_few({n_species})'

        labels = {sid.split('_c')[0] for sid, _ in cds_recs}
        stats_rows.append({'Gene': gene, 'Orthogroup': og,
                           'N_Seqs': len(cds_recs), 'N_Species': n_species,
                           'Status': status,
                           'NMR_in': 'NMR' in labels, 'BATS_in': 'BATS' in labels})

    log_df = pd.DataFrame(log_rows)
    log_path = os.path.join(outdir, 'focal_replacement_log.tsv')
    log_df.to_csv(log_path, sep='\t', index=False)
    pd.DataFrame(stats_rows).to_csv(
        os.path.join(outdir, 'cds_assembly_stats.tsv'), sep='\t', index=False)

    ready = sum(1 for r in stats_rows if r['Status'] == 'ready')
    print(f"\nSUMMARY: genes {len(genes)} | CDS ready {ready} | "
          f"too few {len(genes) - ready}")
    print(f"\n  FOCAL ACTIONS (from {log_path}):")
    for sp in ('NMR', 'BATS'):
        sub = log_df[log_df['Species'] == sp]
        print(f"    {sp}:")
        for act, n in sub['Action'].value_counts().items():
            print(f"       {act:18s}: {n}")
    print("\n  FLAGS:")
    for flag, n in log_df['Flag'].value_counts().items():
        print(f"    {flag:16s}: {n}")
    print(f"\n  cds_per_gene/ -> {outdir}/ (ready for part 7)")


if __name__ == '__main__':
    main()
