#!/usr/bin/env python3
"""Part 9: scheme restriction -> outlier removal -> ortholog selection ->
tree -> labelling -> RELAX, for one gene set.

Schemes come from schemes/{set}.tsv (Scheme, Gene, Test_Species,
Reference_Species): Gene '*' applies a scheme to every gene (atlas sets);
explicit gene rows give per-gene contrasts (validation). Reference_Species
is a comma list or ALL_OTHERS (every aligned species outside the Test set).

Per gene and scheme:
  1. scheme restriction: keep sequences whose species is in Test+Reference
     (copies get unique ids {SP}_{GENE}, {SP}_{GENE}_c2, ...);
  2. protein tree (FastTree -gamma) from all retained copies;
  3. outlier removal (--outlier_removal on, the default): drop non-Test tips
     whose terminal branch exceeds max(95th percentile, Q3 + 1.5*IQR) of the
     gene's terminal branches (skipped when fewer than 5 branches) - the
     part6 criterion;
  4. ortholog selection: per species with multiple surviving copies, keep
     the shortest terminal branch. Exempt from selection (--exempt):
       focal (default): only NMR/BATS keep all copies - part6's original
         scope ("degradation is the signal" belongs to the curated focal
         species; other species' extra copies are paralogs, not signal)
       test: every Test species keeps all copies (the 11f behaviour)
     Outlier removal always exempts Test AND focal species - a long test
     branch is the measurement and is never culled. INTEGRITY: after
     selection the run halts if any non-exempt species holds >1 sequence;
  5. rebuild the tree from the selected set (tips == alignment), label
     {Test}/{Reference} by exact leaf token, run RELAX (frozen invocation:
     CPU=1, --models All, defaults otherwise), resumable, 3600 s timeout.

K and p are extracted to relax_results.tsv per scheme; classification and
the artifact thresholds live in part 10 only.

Alignments: results/08_clipkit/{set}/clipkit_trimmed/. Every set, targets17
included, has its own selection and alignment chain (5d -> 6 -> 7 -> 8);
overlapping alignments can be seeded by copying from a sibling set.

Two modes, one script:
  BATCH/COLLATOR (default): runs every gene x scheme for the set, skipping
    existing JSONs, then writes each scheme's complete relax_results.tsv and
    run_status.tsv ledger.
  WORKER (HPC arrays): --scheme S --only GENE runs one job. NOTE: in worker
    mode each task rewrites the scheme's TSV/ledger with only its own row
    (last writer wins) - JSONs are the ground truth, and one batch-mode
    collation pass after the array rebuilds the tables complete. Failed
    genes leave error files but no JSON, so every collation re-attempts
    them: collate each set ONCE, then treat its tables as final.

Usage: python3 part9_relax.py --set validation --threads 20
       python3 part9_relax.py --set sperm --threads 20 [--outlier_removal off]
"""

import argparse
import glob
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

OUT_BASE = 'results/09_relax'
MODELS = 'All'        # RELAX model set: 'All' fits the general-descriptive model
                      # (per-branch omega, ~44 free parameters) and is the source
                      # of the NaN failures and unstable K; 'Minimal' fits only
                      # null + alternative and is what the whole-genome arm used.
STARTING_POINTS = 1   # RELAX optimiser restarts; >1 guards against local optima
                      # (HyPhy issue #1581). 1 = the frozen default used for all
                      # results to date; raise only for stability checks.
FOCAL_EXEMPT = {'NMR', 'BATS'}   # part6's original selection exemption
MIN_TEST = 1
MIN_REF = 4
RELAX_TIMEOUT = 3600

CODON = {
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
 'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
 'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
 'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
 'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
 'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}


def read_fasta(p):
    seqs, cur, buf = [], None, []
    with open(p) as f:
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


def write_fasta(p, seqs, wrap=60):
    with open(p, 'w') as f:
        for sid, seq in seqs:
            f.write(f">{sid}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i:i + wrap] + "\n")


def translate(seq):
    out = []
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i + 3].upper()
        out.append('-' if c == '---' else ('X' if '-' in c else CODON.get(c, 'X')))
    return ''.join(out)


def species_of(label):
    return label.split('_')[0] if '_' in label else label


def terminal_branch_lengths(newick):
    bl = {}
    for m in re.finditer(r'([A-Za-z0-9_]+):([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', newick):
        bl[m.group(1)] = float(m.group(2))
    return bl


def percentile(sorted_vals, p):
    if not sorted_vals:
        return 0.0
    k = (len(sorted_vals) - 1) * p / 100.0
    f, c = int(k), min(int(k) + 1, len(sorted_vals) - 1)
    return sorted_vals[f] + (sorted_vals[c] - sorted_vals[f]) * (k - f)


def outlier_tips(newick, kept, test_set):
    """part6 criterion: terminal branch > max(95th pct, Q3 + 1.5*IQR);
    skipped when fewer than 5 branches; Test species exempt."""
    bl = terminal_branch_lengths(newick)
    lens = sorted(bl.get(uid) for uid, _, _ in kept if bl.get(uid) is not None)
    if len(lens) < 5:
        return set()
    q1, q3 = percentile(lens, 25), percentile(lens, 75)
    threshold = max(percentile(lens, 95), q3 + 1.5 * (q3 - q1))
    return {uid for uid, _, sp in kept
            if sp not in test_set and bl.get(uid, 0.0) > threshold}


def select_orthologs(newick, kept, test_set):
    """Shortest terminal branch per non-Test species; Test keeps all copies."""
    bl = terminal_branch_lengths(newick)
    by_species = {}
    for uid, _, sp in kept:
        by_species.setdefault(sp, []).append(uid)
    keep = set()
    for sp, uids in by_species.items():
        if sp in test_set:
            keep.update(uids)
        elif len(uids) == 1:
            keep.add(uids[0])
        else:
            keep.add(sorted(uids, key=lambda u: bl.get(u, float('inf')))[0])
    return keep


def build_tree(aa_records, workdir, tag):
    aa_tmp = os.path.join(workdir, f"{tag}_AA.tmp.fa")
    write_fasta(aa_tmp, aa_records)
    try:
        r = subprocess.run(['fasttree', '-gamma', '-quiet', aa_tmp],
                           capture_output=True, text=True,
                           stdin=subprocess.DEVNULL, timeout=600)
        if r.returncode != 0 or not r.stdout.strip():
            r = subprocess.run(['FastTree', '-gamma', '-quiet', aa_tmp],
                               capture_output=True, text=True,
                               stdin=subprocess.DEVNULL, timeout=600)
        if r.returncode != 0 or not r.stdout.strip():
            return None, (r.stderr or '')[:120]
        return r.stdout.strip().replace("'", ""), None
    except Exception as e:
        return None, str(e)[:120]
    finally:
        if os.path.exists(aa_tmp):
            os.remove(aa_tmp)


def process_one(args):
    (gene, aln_file, scheme, test_set, ref_spec, dirs, remove_outliers,
     select_exempt) = args
    out_json = os.path.join(dirs['json'], f"{gene}_RELAX.json")
    err_file = os.path.join(dirs['errors'], f"{gene}_error.txt")
    if os.path.exists(out_json) and os.path.getsize(out_json) > 100:
        return (gene, 'exists', '')

    seqs = read_fasta(aln_file)
    all_species = {species_of(sid) for sid, _ in seqs}
    if ref_spec == ['ALL_OTHERS']:
        keep_set = all_species
    else:
        keep_set = test_set | set(ref_spec)

    kept, sp_seen = [], {}
    for sid, seq in seqs:
        sp = species_of(sid)
        if sp not in keep_set:
            continue
        n = sp_seen.get(sp, 0)
        sp_seen[sp] = n + 1
        uid = f"{sp}_{gene}" if n == 0 else f"{sp}_{gene}_c{n + 1}"
        kept.append((uid, seq, sp))

    if len({sp for _, _, sp in kept}) < 4:
        return (gene, 'too_few_species', f'{len({s for _, _, s in kept})}')

    tree1, err = build_tree([(u, translate(s)) for u, s, _ in kept],
                            dirs['tree'], gene)
    if tree1 is None:
        with open(err_file, 'w') as f:
            f.write(f"tree build (all copies) failed: {err}\n")
        return (gene, 'tree_fail', err)

    n_outliers = 0
    if remove_outliers:
        out_uids = outlier_tips(tree1, kept, test_set | FOCAL_EXEMPT)
        n_outliers = len(out_uids)
        kept = [k for k in kept if k[0] not in out_uids]

    selected = select_orthologs(tree1, kept, select_exempt)
    kept = [k for k in kept if k[0] in selected]

    # integrity: after selection, no non-exempt species may hold >1 sequence
    per_sp = defaultdict(int)
    for _, _, sp in kept:
        per_sp[sp] += 1
    bad = [sp for sp, n in per_sp.items() if sp not in select_exempt and n > 1]
    assert not bad, f"{gene}/{scheme}: multiple copies survive selection: {bad}"

    test_n = sum(1 for _, _, sp in kept if sp in test_set)
    ref_n = len(kept) - test_n
    if test_n < MIN_TEST:
        return (gene, 'no_test', f'{test_n}')
    if ref_n < MIN_REF:
        return (gene, 'too_few_ref', f'{ref_n}')

    pruned_nt = os.path.join(dirs['pruned'], f"{gene}_NT.fa")
    write_fasta(pruned_nt, [(u, s) for u, s, _ in kept])

    tree2, err = build_tree([(u, translate(s)) for u, s, _ in kept],
                            dirs['tree'], gene + '_2')
    if tree2 is None:
        with open(err_file, 'w') as f:
            f.write(f"tree rebuild (selected set) failed: {err}\n")
        return (gene, 'tree_fail', err)

    labeled = tree2
    for uid, _, sp in kept:
        tag = '{Test}' if sp in test_set else '{Reference}'
        labeled = re.sub(r'(?<![\w])' + re.escape(uid) + r'(?=[:,\);])',
                         uid + tag, labeled)
    labeled_file = os.path.join(dirs['labeled'], f"{gene}.tree")
    with open(labeled_file, 'w') as f:
        f.write(labeled + '\n')

    cmd = ['hyphy', 'CPU=1', 'relax',
           '--alignment', pruned_nt, '--tree', labeled_file,
           '--test', 'Test', '--reference', 'Reference',
           '--models', MODELS,
           '--output', out_json]
    if STARTING_POINTS > 1:
        cmd += ['--starting-points', str(STARTING_POINTS)]
    try:
        rr = subprocess.run(cmd, capture_output=True, text=True,
                            stdin=subprocess.DEVNULL, timeout=RELAX_TIMEOUT)
        if os.path.exists(out_json) and os.path.getsize(out_json) > 100:
            if os.path.exists(err_file):
                os.remove(err_file)      # stale receipt from a rescued earlier attempt
            return (gene, 'success', f'outliers_removed={n_outliers}')
        with open(err_file, 'w') as f:
            f.write("CMD: " + ' '.join(cmd) + f"\nReturn code: {rr.returncode}\n\n"
                    f"STDERR:\n{rr.stderr}\n\nSTDOUT:\n{rr.stdout}\n")
        return (gene, 'relax_fail', (rr.stderr or '')[:150])
    except subprocess.TimeoutExpired:
        with open(err_file, 'w') as f:
            f.write("TIMEOUT >" + str(RELAX_TIMEOUT) + "s\nCMD: " + ' '.join(cmd) + "\n")
        return (gene, 'relax_timeout', f'>{RELAX_TIMEOUT}s')
    except Exception as e:
        with open(err_file, 'w') as f:
            f.write("EXCEPTION: " + str(e) + "\n")
        return (gene, 'relax_error', str(e)[:150])


def parse_json(p):
    try:
        d = json.load(open(p))
        tr = d.get('test results', {})
        K = (tr.get('relaxation or intensification parameter')
             or tr.get('relaxation test parameter') or tr.get('K'))
        pv = tr.get('p-value', tr.get('p_value'))
        return K, pv, tr.get('LRT')
    except Exception:
        return None, None, None


def load_schemes(path):
    df = pd.read_csv(path, sep='\t')
    schemes = defaultdict(list)
    for _, r in df.iterrows():
        schemes[r['Scheme']].append(
            (str(r['Gene']),
             {s.strip() for s in str(r['Test_Species']).split(',') if s.strip()},
             [s.strip() for s in str(r['Reference_Species']).split(',') if s.strip()]))
    return schemes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--set', required=True,
                    help='validation | targets17 | sperm | ddr')
    ap.add_argument('--threads', type=int, default=20)
    ap.add_argument('--outlier_removal', default='on', choices=['on', 'off'])
    ap.add_argument('--scheme', default=None,
                    help='run only this scheme from the set config (HPC array use)')
    ap.add_argument('--only', default=None,
                    help='comma-separated gene names (HPC array use: one gene)')
    ap.add_argument('--models', default='All', choices=['All', 'Minimal'],
                    help="RELAX model set. 'All' is the frozen spec used for "
                         "every result so far; 'Minimal' drops the "
                         "general-descriptive fit (faster, far more stable)")
    ap.add_argument('--starting_points', type=int, default=1,
                    help='RELAX optimiser restarts (--starting-points). 1 keeps '
                         'the frozen spec used for every result so far; 10 '
                         'matches the whole-genome arm and is slower')
    ap.add_argument('--relax_timeout', type=int, default=3600,
                    help='seconds per RELAX run [3600]; raise for giants on HPC')
    ap.add_argument('--exempt', default='focal', choices=['focal', 'test'],
                    help="selection exemption scope: focal = only NMR/BATS keep "
                         "all copies (part6 original; default); test = every "
                         "Test species does (11f behaviour)")
    args = ap.parse_args()

    aln_set = 'sperm' if args.set == 'targets17' else args.set   # targets17 reads sperm's alignments
    aln_dir = os.path.join('results/08_clipkit', aln_set, 'clipkit_trimmed')
    schemes_file = os.path.join('schemes', f'{args.set}.tsv')
    for p, hint in ((aln_dir, 'run part 8'), (schemes_file, 'place the config')):
        if not os.path.exists(p):
            sys.exit(f"ERROR: {p} missing ({hint})")
    if not shutil.which('hyphy'):
        sys.exit("ERROR: hyphy not in PATH")
    if not shutil.which('fasttree') and not shutil.which('FastTree'):
        sys.exit("ERROR: fasttree not in PATH")

    global RELAX_TIMEOUT, STARTING_POINTS, MODELS
    RELAX_TIMEOUT = args.relax_timeout
    STARTING_POINTS = args.starting_points
    MODELS = args.models

    genes = sorted(os.path.basename(f).replace('_NT.clipkit.fa', '')
                   for f in glob.glob(os.path.join(aln_dir, '*_NT.clipkit.fa')))
    if args.only:
        wanted_only = {g.strip() for g in args.only.split(',')}
        genes = [g for g in genes if g in wanted_only]
    if args.set == 'targets17':
        wanted = set(pd.read_csv('gene_lists/targets17.tsv', sep='\t')['SYMBOL'])
        genes = [g for g in genes if g in wanted]

    schemes = load_schemes(schemes_file)
    if args.scheme:
        if args.scheme not in schemes:
            sys.exit(f"ERROR: scheme '{args.scheme}' not in {schemes_file}")
        schemes = {args.scheme: schemes[args.scheme]}
    remove_outliers = args.outlier_removal == 'on'
    print(f"Part 9 [{args.set}]: RELAX  |  outlier removal: "
          f"{args.outlier_removal}  |  selection exemption: {args.exempt}  |  "
          f"models: {args.models}  |  starting-points: {args.starting_points}  |  "
          f"{len(genes)} genes, {len(schemes)} scheme(s), "
          f"alignments from {aln_dir}")

    for scheme, rows in schemes.items():
        sdir = os.path.join(OUT_BASE, args.set, scheme)
        dirs = {k: os.path.join(sdir, k) for k in ('pruned', 'tree', 'labeled', 'json', 'errors')}
        for d in dirs.values():
            os.makedirs(d, exist_ok=True)

        star = next(((t, r) for g, t, r in rows if g == '*'), None)
        per_gene = {g: (t, r) for g, t, r in rows if g != '*'}
        jobs = []
        for g in genes:
            spec = per_gene.get(g, star)
            if spec is None:
                continue
            test_set, ref_spec = spec
            select_exempt = (FOCAL_EXEMPT if args.exempt == 'focal'
                             else test_set | FOCAL_EXEMPT)
            jobs.append((g, os.path.join(aln_dir, f"{g}_NT.clipkit.fa"),
                         scheme, test_set, ref_spec, dirs, remove_outliers,
                         select_exempt))

        done = sum(1 for j in jobs
                   if os.path.exists(os.path.join(dirs['json'], f"{j[0]}_RELAX.json"))
                   and os.path.getsize(os.path.join(dirs['json'], f"{j[0]}_RELAX.json")) > 100)
        print(f"\nScheme: {scheme}  ({len(jobs)} genes, {done} already done)")

        cnt = defaultdict(int)
        status_rows = []
        with ProcessPoolExecutor(max_workers=args.threads) as exe:
            futs = {exe.submit(process_one, j): j[0] for j in jobs}
            c = 0
            for fut in as_completed(futs):
                g, st, msg = fut.result()
                status_rows.append((g, st, msg))
                cnt[st] += 1
                c += 1
                if c % 50 == 0 or c == len(jobs):
                    ok = cnt['success'] + cnt['exists']
                    print(f"  {c}/{len(jobs)}  ok={ok} no_test={cnt['no_test']} "
                          f"fail={cnt['relax_fail'] + cnt['tree_fail']} "
                          f"timeout={cnt['relax_timeout']}")

        rows_out = []
        for j in jobs:
            g = j[0]
            jp = os.path.join(dirs['json'], f"{g}_RELAX.json")
            if os.path.exists(jp) and os.path.getsize(jp) > 100:
                K, pv, lrt = parse_json(jp)
                rows_out.append({'Gene': g, 'K': K, 'p_value': pv, 'LRT': lrt,
                                 'Scheme': scheme, 'Set': args.set})
        if rows_out:
            pd.DataFrame(rows_out).sort_values('p_value', na_position='last').to_csv(
                os.path.join(sdir, 'relax_results.tsv'), sep='\t', index=False)
        with open(os.path.join(sdir, 'run_status.tsv'), 'w') as f:
            f.write("Gene\tStatus\tMessage\n")
            for g, st, msg in sorted(status_rows):
                f.write(f"{g}\t{st}\t{msg}\n")
        print("  status: " + ', '.join(f'{k}={v}' for k, v in
                                       sorted(cnt.items(), key=lambda x: -x[1])))
        print(f"  per-gene ledger -> {sdir}/run_status.tsv  "
              f"(failures detailed in {sdir}/errors/)")
        print(f"  -> {sdir}/relax_results.tsv")


if __name__ == '__main__':
    main()
