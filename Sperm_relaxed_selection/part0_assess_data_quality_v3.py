#!/usr/bin/env python3
"""Quick QC of the reference: for each species in
species_list_17missinggene.csv, locate its genome, annotation, CDS and
protein files and tabulate basic counts into data_quality_summary.tsv.

Subfolders are: genomes/, annotations/, cds/, proteins/.

Usage: python3 part0_assess_data_quality_v3.py [threads]

"""

import glob
import gzip
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

import pandas as pd
from Bio import SeqIO

GENOME_EXTS = ['.fna', '.fa', '.fasta', '.fna.gz', '.fa.gz', '.fasta.gz']
CDS_EXTS = GENOME_EXTS
ANNOT_EXTS = ['.gff', '.gff3', '.gtf', '.gff.gz', '.gff3.gz', '.gtf.gz']
PROTEIN_EXTS = ['.faa', '.fa', '.fasta', '.aa', '.faa.gz', '.fa.gz', '.fasta.gz', '.aa.gz']

COLUMNS = ['species', 'version', 'genome_sequences', 'genome_size',
           'genes', 'transcripts', 'cds_sequences', 'protein_sequences']


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def fasta_stats(path):
    """Sequence count and total length, streaming (genome files are large)."""
    n = total = 0
    with open_text(path) as fh:
        for rec in SeqIO.parse(fh, 'fasta'):
            n += 1
            total += len(rec.seq)
    return n, total


def gff_stats(path):
    counts = defaultdict(int)
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                counts[fields[2]] += 1
    return counts['gene'], counts['transcript'] + counts['mRNA']


def find_file(directory, patterns, exts):
    """Return the first file matching *pattern*ext, trying patterns in order.
    Derived files (CDS, proteins) often drop the annotation-release suffix
    from the assembly name, hence the progressively looser patterns."""
    for pattern in patterns:
        hits = []
        for ext in exts:
            hits.extend(glob.glob(f"{directory}/*{pattern}*{ext}"))
        if hits:
            return hits[0]
    return None


def survey_species(job):
    """Worker: locate and profile one species' four input files."""
    species, version = job

    patterns = [version]
    if '.' in version:
        patterns.append(version.rsplit('.', 1)[0])  # e.g. drop a trailing ".115"
    base = version.split('.')[0]
    if base not in patterns:
        patterns.append(base)

    def locate_and_count(directory, exts, parser):
        path = find_file(directory, patterns, exts)
        if path is None:
            return 'NOT FOUND', (0, 0)
        try:
            return os.path.basename(path), parser(path)
        except Exception as exc:
            print(f"WARNING {species}: failed to parse {path}: {exc}", file=sys.stderr)
            return os.path.basename(path), (0, 0)

    g_file, (g_n, g_len) = locate_and_count('genomes', GENOME_EXTS, fasta_stats)
    a_file, (genes, transcripts) = locate_and_count('annotations', ANNOT_EXTS, gff_stats)
    c_file, (c_n, _) = locate_and_count('cds', CDS_EXTS, fasta_stats)
    p_file, (p_n, _) = locate_and_count('proteins', PROTEIN_EXTS, fasta_stats)

    return {
        'species': species, 'version': version,
        'genome_file': g_file, 'genome_sequences': g_n, 'genome_size': g_len,
        'annotation_file': a_file, 'genes': genes, 'transcripts': transcripts,
        'cds_file': c_file, 'cds_sequences': c_n,
        'protein_file': p_file, 'protein_sequences': p_n,
    }


def main():
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else 8

    species = pd.read_csv('species_list_17missinggene.csv')
    species['Short_name'] = species['Short_name'].str.strip()
    species['Version'] = species['Version'].str.strip()
    jobs = list(zip(species['Short_name'], species['Version']))

    print(f"Surveying {len(jobs)} species with {threads} processes")
    results = []
    with Pool(threads) as pool:
        for i, row in enumerate(pool.imap_unordered(survey_species, jobs), 1):
            print(f"[{i:>2}/{len(jobs)}] {row['species']:<14} genome={row['genome_file']}  "
                  f"genes={row['genes']}  cds={row['cds_sequences']}  proteins={row['protein_sequences']}")
            missing = [k.split('_')[0] for k in
                       ('genome_file', 'annotation_file', 'cds_file', 'protein_file')
                       if row[k] == 'NOT FOUND']
            if missing:
                print(f"          missing: {', '.join(missing)}")
            results.append(row)

    results.sort(key=lambda r: r['species'])
    summary = pd.DataFrame(results)[COLUMNS]
    summary.to_csv('data_quality_summary.tsv', sep='\t', index=False)

    print()
    print(summary.to_string(index=False))
    print("\nWrote data_quality_summary.tsv")


if __name__ == '__main__':
    main()
