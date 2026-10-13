#!/usr/bin/env python3
"""Strip characters that DIAMOND/OrthoFinder cannot handle from the protein
sequences: '.' placeholders and '*' stop codons.
Headers and record order are unchanged.

Usage: python3 part2b_clean_proteins_v2.py
"""

import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq

INPUT_DIR = 'results/02a_renamed_proteins'
OUTPUT_DIR = 'results/02b_cleaned_proteins'


def main():
    if not os.path.isdir(INPUT_DIR):
        sys.exit(f"ERROR: {INPUT_DIR} not found - run part2a first")

    files = sorted(f for f in os.listdir(INPUT_DIR) if f.endswith('.fa'))
    if not files:
        sys.exit(f"ERROR: no .fa files in {INPUT_DIR}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for filename in files:
        records, modified = [], 0
        for rec in SeqIO.parse(os.path.join(INPUT_DIR, filename), 'fasta'):
            clean = str(rec.seq).replace('.', '').replace('*', '')
            if len(clean) != len(rec.seq):
                modified += 1
            rec.seq = Seq(clean)
            records.append(rec)
        SeqIO.write(records, os.path.join(OUTPUT_DIR, filename), 'fasta')
        print(f"{filename:<15} {len(records)} sequences ({modified} contained '.' or '*')")

    print(f"\nCleaned files written to {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
