#!/usr/bin/env python3
"""Part 8: ClipKIT codon-aware trimming for one gene set.

This script is parallel and resumable (ClipKIT is fast; this runs in minutes). 
Per gene: replace MACSE special characters ('?', '!') and 
ambiguity codes ('*', N, X) with gaps, run ClipKIT with --codon 
(maintains the reading frame), trim to a codon boundary, 
replace stop codons with '---', validate (>=4 sequences, >=50
codons, no remaining stops or special characters).

Stop replacement is a single pass: with the alignment frame-aligned,
swapping a stop triplet for '---' is length-preserving and cannot create
new stops; a residual count is still checked and fails loudly.

Inputs:  results/07_macse/{set}/macse_aligned/*_NT.aln   (from part 7)
Outputs: results/08_clipkit/{set}/clipkit_trimmed/{GENE}_NT.clipkit.fa
         results/08_clipkit/{set}/clipkit_summary.tsv

Usage:
  python3 part8_clipkit.py --set sperm --threads 44
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

def set_paths(gene_set):
    base = os.path.join('results/08_clipkit', gene_set)
    return os.path.join('results/07_macse', gene_set, 'macse_aligned'), base, os.path.join(base, 'clipkit_trimmed')

SPECIAL_CHARS = {'?', '!', '*', 'N', 'n', 'X', 'x'}
STOPS = {'TAA', 'TAG', 'TGA'}
MIN_SEQS = 4       # atlas minimum (the 17-gene part8 used 8)
MIN_CODONS = 50


def read_fasta(path):
    seqs, cur, buf = [], None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if cur is not None:
                    seqs.append((cur, ''.join(buf)))
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            seqs.append((cur, ''.join(buf)))
    return seqs


def write_fasta(path, seqs, wrap=60):
    with open(path, 'w') as f:
        for sid, seq in seqs:
            f.write(f">{sid}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i:i + wrap] + "\n")


def sanitize_id(name):
    """HyPhy-safe IDs. (Atlas IDs are already bare species labels, so this
    is a no-op in practice; kept as a guard.)"""
    s = re.sub(r'[^A-Za-z0-9_]', '_', name)
    return re.sub(r'_+', '_', s).strip('_')


def count_stops(seqs):
    n = 0
    for _, q in seqs:
        for i in range(0, len(q) - 2, 3):
            if q[i:i + 3].upper() in STOPS:
                n += 1
    return n


def process_gene(args):
    gene, aln_file, out_dir, mode, tmp_root = args
    out_file = os.path.join(out_dir, f"{gene}_NT.clipkit.fa")
    if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
        return (gene, 'exists', '')

    if not os.path.exists(aln_file):
        return (gene, 'missing', 'no aln')
    seqs = read_fasta(aln_file)
    if len(seqs) < MIN_SEQS:
        return (gene, 'too_few_in', f'{len(seqs)} seqs')

    cleaned = [(sid, ''.join('-' if c in SPECIAL_CHARS else c for c in seq))
               for sid, seq in seqs]

    tmpd = tempfile.mkdtemp(prefix=f'clipkit_{gene}_', dir=tmp_root)
    try:
        clean_in = os.path.join(tmpd, 'in.fa')
        write_fasta(clean_in, cleaned)
        clip_out = os.path.join(tmpd, 'out.fa')

        cmd = ['clipkit', clean_in, '-m', mode, '--codon', '-o', clip_out]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True,
                               stdin=subprocess.DEVNULL, timeout=600)
            if r.returncode != 0 or not os.path.exists(clip_out):
                return (gene, 'clipkit_fail', (r.stderr or '')[:150])
        except subprocess.TimeoutExpired:
            return (gene, 'clipkit_timeout', '')
        except FileNotFoundError:
            return (gene, 'clipkit_missing', 'clipkit not installed')

        seqs2 = read_fasta(clip_out)
        if not seqs2:
            return (gene, 'empty_after_clip', '')

        # codon boundary
        L = len(seqs2[0][1])
        if L % 3:
            seqs2 = [(s, q[:L - L % 3]) for s, q in seqs2]

        # single-pass stop replacement (see docstring)
        new = []
        for sid, q in seqs2:
            cods = []
            for i in range(0, len(q), 3):
                c = q[i:i + 3]
                cods.append('---' if len(c) == 3 and c.upper() in STOPS else c)
            new.append((sid, ''.join(cods)))
        seqs2 = new

        # validate
        L = len(seqs2[0][1]) if seqs2 else 0
        if L % 3 != 0:
            return (gene, 'bad_frame', f'{L}')
        if L < MIN_CODONS * 3:
            return (gene, 'too_short', f'{L}bp')
        remaining = count_stops(seqs2)
        if remaining:   # cannot happen for a frame-aligned input
            return (gene, 'stops_remain', f'{remaining}')

        final = [(sanitize_id(sid), q) for sid, q in seqs2]
        write_fasta(out_file, final)
        return (gene, 'success', f'{len(seqs)}seq {L}bp')
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--set', required=True, help='validation | sperm | ddr')
    ap.add_argument('--threads', type=int, default=44)
    ap.add_argument('--mode', default='gappy',
                    help='ClipKIT mode [gappy] (matches the 17-gene part8)')
    args = ap.parse_args()

    if not shutil.which('clipkit'):
        sys.exit("ERROR: clipkit not in PATH (activate conda env)")
    ALN_DIR, BASE, OUT_DIR = set_paths(args.set)
    os.makedirs(OUT_DIR, exist_ok=True)

    aln_files = sorted(glob.glob(os.path.join(ALN_DIR, '*_NT.aln')))
    genes = [os.path.basename(f).replace('_NT.aln', '') for f in aln_files]
    print(f"Part 10e: ClipKIT --codon trimming of {len(genes)} alignments "
          f"({args.mode} mode, {args.threads} threads)")

    tmp_root = tempfile.mkdtemp(prefix='clipkit_root_')
    jobs = [(g, os.path.join(ALN_DIR, f"{g}_NT.aln"), OUT_DIR, args.mode, tmp_root)
            for g in genes]

    rows = []
    counts = defaultdict(int)
    try:
        with ProcessPoolExecutor(max_workers=args.threads) as exe:
            futs = {exe.submit(process_gene, j): j[0] for j in jobs}
            done = 0
            for fut in as_completed(futs):
                gene, status, msg = fut.result()
                rows.append((gene, status, msg))
                counts[status] += 1
                done += 1
                if done % 100 == 0 or done == len(jobs):
                    print(f"  {done}/{len(jobs)}  "
                          f"(ok={counts['success'] + counts['exists']})")
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)

    with open(os.path.join(BASE, 'clipkit_summary.tsv'), 'w') as f:
        f.write("Gene\tStatus\tMessage\n")
        for r in sorted(rows):
            f.write('\t'.join(r) + '\n')

    n_ok = counts['success'] + counts['exists']
    print(f"\ntrimmed OK: {n_ok}/{len(genes)}")
    for st, n in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {st:18s}: {n}")
    n_out = len(glob.glob(os.path.join(OUT_DIR, '*_NT.clipkit.fa')))
    print(f"output files: {n_out} in {OUT_DIR}/")


if __name__ == '__main__':
    main()
