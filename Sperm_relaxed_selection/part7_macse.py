#!/usr/bin/env python3
"""Part 7: MACSE codon-aware alignment for one gene set.

This script resumes on restart (genes with a non-empty _NT.aln
are skipped). Genes exceeding the per-gene timeout are killed and logged to
macse_timeouts.tsv so the few giants can be rerun separately (longer
--timeout_hours, --only GENE1,GENE2) without dominating wall time.
ClipKIT / trees / RELAX are separate downstream scripts (10e, 10f).

Inputs:  results/06_cds/{set}/cds_per_gene/*.fa   (from part 6)
Outputs (results/07_macse/{set}/):
  macse_aligned/{GENE}_NT.aln, {GENE}_AA.aln
  macse_aligned/{GENE}_macse.log      per-gene stderr on failure
  macse_run_log.tsv                   every gene: status, seconds, n_seq
  macse_timeouts.tsv, macse_failures.tsv

Usage:
  python3 part7_macse.py --set sperm --threads 40 --timeout_hours 3
  python3 part7_macse.py --set sperm --threads 1 --timeout_hours 24 --xmx 80g --only FSIP2,SYNE1
"""

import argparse
import glob
import re
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

def set_paths(gene_set):
    base = os.path.join('results/07_macse', gene_set)
    return os.path.join('results/06_cds', gene_set, 'cds_per_gene'), base, os.path.join(base, 'macse_aligned')

TOTAL_RAM_GB = 250
MACSE_RAM_GB = 4      # per job; caps concurrency


def count_seqs(fa):
    try:
        return sum(1 for line in open(fa) if line.startswith('>'))
    except Exception:
        return 0


def run_macse_one(args):
    """Align one gene. Status: success, exists, timeout, failed, error."""
    gene, cds_fa, out_dir, timeout_s, jar, xmx = args
    out_nt = os.path.join(out_dir, f"{gene}_NT.aln")
    out_aa = os.path.join(out_dir, f"{gene}_AA.aln")
    log_f = os.path.join(out_dir, f"{gene}_macse.log")

    n_seq = count_seqs(cds_fa)
    if os.path.exists(out_nt) and os.path.getsize(out_nt) > 0:
        return (gene, 'exists', 0.0, n_seq, '')

    cmd = ['java', f'-Xmx{xmx}', '-jar', jar,
           '-prog', 'alignSequences', '-seq', cds_fa,
           '-out_NT', out_nt, '-out_AA', out_aa]
    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           stdin=subprocess.DEVNULL, timeout=timeout_s)
        dt = time.time() - t0
        if r.returncode == 0 and os.path.exists(out_nt) and os.path.getsize(out_nt) > 0:
            return (gene, 'success', dt, n_seq, '')
        with open(log_f, 'w') as f:
            f.write("CMD: " + ' '.join(cmd) + "\n\nSTDERR:\n" + (r.stderr or '')
                    + "\n\nSTDOUT:\n" + (r.stdout or ''))
        for p in (out_nt, out_aa):   # remove partial output too, or resume skips a broken alignment
            if os.path.exists(p):
                os.remove(p)
        return (gene, 'failed', dt, n_seq, (r.stderr or 'nonzero exit')[:200])
    except subprocess.TimeoutExpired:
        dt = time.time() - t0
        for p in (out_nt, out_aa):   # clean partial output so resume reruns cleanly
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
        return (gene, 'timeout', dt, n_seq, f'>{timeout_s}s')
    except Exception as e:
        return (gene, 'error', time.time() - t0, n_seq, str(e)[:200])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--set', required=True, help='validation | sperm | ddr')
    ap.add_argument('--threads', type=int, default=40,
                    help=f'parallel MACSE jobs (~{MACSE_RAM_GB}GB each)')
    ap.add_argument('--timeout_hours', type=float, default=3.0,
                    help='per-gene timeout in hours [3]')
    ap.add_argument('--only', default=None,
                    help='comma-separated gene names to run (others skipped)')
    ap.add_argument('--xmx', default='4g',
                    help='Java heap per gene [4g]; use 80g for the giants '
                         '(the conda macse wrapper hardcodes -Xmx2g, so the '
                         'jar is driven directly)')
    ap.add_argument('--jar', default=None,
                    help='path to macse jar [auto-detect under $CONDA_PREFIX]')
    args = ap.parse_args()

    if not shutil.which('java'):
        sys.exit("ERROR: 'java' not in PATH (activate the conda env)")
    jars = sorted(glob.glob(os.path.join(os.environ.get('CONDA_PREFIX', ''),
                                         'share', 'macse-*', 'macse*.jar')))
    if args.jar:
        jar = args.jar
    elif jars:
        jar = jars[-1]
    else:
        sys.exit("ERROR: MACSE jar not found under $CONDA_PREFIX/share/macse-*/ "
                 "- pass --jar /path/to/macse_v2.07.jar")
    print(f"  MACSE jar: {jar}  heap: -Xmx{args.xmx}")
    CDS_DIR, BASE, OUT_DIR = set_paths(args.set)
    os.makedirs(OUT_DIR, exist_ok=True)
    timeout_s = int(args.timeout_hours * 3600)
    m = re.match(r'^(\d+)g$', args.xmx.lower())
    per_job_gb = int(m.group(1)) if m else MACSE_RAM_GB
    max_workers = min(args.threads, max(1, TOTAL_RAM_GB // max(per_job_gb, 1)))

    cds_files = sorted(glob.glob(os.path.join(CDS_DIR, '*.fa')))
    genes = [os.path.basename(f)[:-3] for f in cds_files]
    if args.only:
        want = {g.strip() for g in args.only.split(',')}
        genes = [g for g in genes if g in want]

    done = sum(1 for g in genes
               if os.path.exists(os.path.join(OUT_DIR, f"{g}_NT.aln"))
               and os.path.getsize(os.path.join(OUT_DIR, f"{g}_NT.aln")) > 0)

    print(f"Part 7 [{args.set}]: MACSE alignment")
    print(f"  genes: {len(genes)}   already aligned: {done}   "
          f"to align: {len(genes) - done}")
    print(f"  workers: {max_workers}   per-gene timeout: {args.timeout_hours} h")

    jobs = [(g, os.path.join(CDS_DIR, f"{g}.fa"), OUT_DIR, timeout_s, jar, args.xmx) for g in genes]

    run_rows, timeouts, failures = [], [], []
    success = exists = to = fail = 0
    t_start = time.time()

    with ProcessPoolExecutor(max_workers=max_workers) as exe:
        futs = {exe.submit(run_macse_one, j): j[0] for j in jobs}
        completed = 0
        for fut in as_completed(futs):
            gene, status, secs, n_seq, msg = fut.result()
            run_rows.append((gene, status, round(secs, 1), n_seq, msg))
            if status == 'success':
                success += 1
            elif status == 'exists':
                exists += 1
            elif status == 'timeout':
                to += 1
                timeouts.append((gene, n_seq, round(secs, 1)))
            else:
                fail += 1
                failures.append((gene, status, n_seq, msg))
            completed += 1
            if completed % 20 == 0 or completed == len(jobs):
                el = (time.time() - t_start) / 60
                print(f"  {completed}/{len(jobs)}  ok={success} skip={exists} "
                      f"timeout={to} fail={fail}  [{el:.1f} min]")

    with open(os.path.join(BASE, 'macse_run_log.tsv'), 'w') as f:
        f.write("Gene\tStatus\tSeconds\tN_seq\tMessage\n")
        for r in sorted(run_rows):
            f.write('\t'.join(str(x) for x in r) + '\n')
    if timeouts:
        with open(os.path.join(BASE, 'macse_timeouts.tsv'), 'w') as f:
            f.write("Gene\tN_seq\tSeconds\n")
            for r in sorted(timeouts):
                f.write('\t'.join(str(x) for x in r) + '\n')
    if failures:
        with open(os.path.join(BASE, 'macse_failures.tsv'), 'w') as f:
            f.write("Gene\tStatus\tN_seq\tMessage\n")
            for r in sorted(failures):
                f.write('\t'.join(str(x) for x in r) + '\n')

    el = (time.time() - t_start) / 60
    n_aln = len(glob.glob(os.path.join(OUT_DIR, '*_NT.aln')))
    print(f"\nsuccess={success} skipped={exists} timeout={to} failed={fail}  "
          f"[{el:.1f} min]   _NT.aln in dir: {n_aln}")
    if timeouts:
        print(f"{len(timeouts)} timed out (macse_timeouts.tsv) - rerun with longer "
              f"--timeout_hours and --only {','.join(g for g, _, _ in sorted(timeouts)[:8])}")
    if failures:
        print(f"{len(failures)} failed (macse_failures.tsv + per-gene _macse.log)")


if __name__ == '__main__':
    main()
