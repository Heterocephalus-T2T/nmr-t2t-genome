#!/usr/bin/env python3
"""Rename the part1 proteome files to short species codes for OrthoFinder,
which uses the file name as the species identifier. Only file names change;
protein IDs inside the files are left untouched, since they are needed later
for CDS/GFF mapping. species_codes.tsv records the mapping.

Usage: python3 part2a_rename_samples_v2.py
"""

import os
import shutil
import sys

INPUT_DIR = 'results/01_isoform_selection'
OUTPUT_DIR = 'results/02a_renamed_proteins'

SPECIES_CODES = {
    # bathyergid focal species
    'nmr':          'hetg',  # Heterocephalus glaber
    'batsui':       'bats',  # Bathyergus suillus
    'amr':          'fuka',  # Fukomys anselli
    # rodents
    'mouse':        'mous',  # Mus musculus 
    'wwildmouse':   'wwmo',  # Mus spretus
    'norwayrat':    'nrat',  # Rattus norvegicus 
    'gpig':         'gpig',  # Cavia porcellus
    'degu':         'degu',  # Octodon degus 
    'ltchinchilla': 'chin',  # Chinchilla lanigera
    'pvole':        'pvol',  # Microtus ochrogaster
    'nadeermouse':  'nade',  # Peromyscus maniculatus
    'ersquirrel':   'squi',  # Sciurus vulgaris
    'csmus':        'csmu',  # Acomys cahirinus
    # primates
    'human':        'huma',  # Homo sapiens 
    'bonobo':       'bono',  # Pan paniscus
    'gorilla':      'gori',  # Gorilla gorilla gorilla
    'macaque':      'maca',  # Macaca mulatta
    # other mammals
    'rabbit':       'rabb',  # Oryctolagus cuniculus
    'dog':          'dogs',  # Canis lupus familiaris
    'cat':          'cats',  # Felis catus
    'pig':          'pigs',  # Sus scrofa
    'cattle':       'catt',  # Bos taurus
    'elephant':     'elep',  # Loxodonta africana

}


def main():
    if not os.path.isdir(INPUT_DIR):
        sys.exit(f"ERROR: {INPUT_DIR} not found - run part1 first")

    files = sorted(f for f in os.listdir(INPUT_DIR) if f.endswith('.primary.fa'))
    if not files:
        sys.exit(f"ERROR: no .primary.fa files in {INPUT_DIR}")

    unknown = sorted({f.replace('.primary.fa', '') for f in files} - set(SPECIES_CODES))
    if unknown:
        sys.exit(f"ERROR: species not in SPECIES_CODES: {', '.join(unknown)}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    rows = []
    for filename in files:
        species = filename.replace('.primary.fa', '')
        code = SPECIES_CODES[species]
        src = os.path.join(INPUT_DIR, filename)
        shutil.copy2(src, os.path.join(OUTPUT_DIR, f"{code}.fa"))
        with open(src) as fh:
            n_proteins = sum(1 for line in fh if line.startswith('>'))
        rows.append((species, code, n_proteins, filename, f"{code}.fa"))
        print(f"{species:<20} -> {code}.fa  ({n_proteins} proteins)")

    with open(os.path.join(OUTPUT_DIR, 'species_codes.tsv'), 'w') as out:
        out.write('original_name\tshort_code\tprotein_count\toriginal_file\tnew_file\n')
        for r in rows:
            out.write('\t'.join(str(x) for x in r) + '\n')
    print(f"\n{len(rows)} species, {sum(r[2] for r in rows)} proteins total")
    print(f"Wrote {OUTPUT_DIR}/species_codes.tsv")


if __name__ == '__main__':
    main()
