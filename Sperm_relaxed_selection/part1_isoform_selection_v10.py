#!/usr/bin/env python3
"""Select one representative protein isoform per gene for every species.

Strategy by species class:
  NMR, BatSui            customised single-copy TSV (exact transcript-ID match)
  human                  MANE Select
  NCBI species (bonobo,  pre-filter pseudogene/low-quality annotations, then
  gorilla by format)     one isoform per gene by tier (NP_ > XP_), DIAMOND
                         best-hit score to human, then length
  Ensembl species        APPRIS principal / Ensembl_canonical / basic tag,
                         longest as tie-break
  BRAKER/mikado/simple   GFF gene->transcript grouping; best isoform by
                         DIAMOND best-hit to human, the best-annotated
                         reference proteome ('.'-containing sequences
                         penalized), cleanest/longest as fallback

Project subfolder layout: genomes/, cds/, proteins/, annotations/,
species_list_17missinggene.csv, and the single-copy TSVs for NMR/BatSui.

Output: results/01_isoform_selection/<species>.primary.fa
Log:    logs/01_isoform_selection/01_isoform_selection_v10.log

Usage: python3 part1_isoform_selection_v10.py
"""

import glob
import gzip
import os
import re
import shutil
import subprocess
import sys
import traceback
from collections import defaultdict
from datetime import datetime

import pandas as pd
from Bio import SeqIO

# =============================================================================
# CONFIGURATION
# =============================================================================

DIAMOND_THREADS = 30

COLLABORATOR_TSVS = {
    'nmr':    'mHetGlaPri.single.copy.tsv',
    'batsui': 'BatSui.single.copy.tsv',
}

# Species known to need strict NCBI-style filtering. Others with NP_/XP_
NCBI_SPECIES = {'bonobo'}

# Description/ID patterns that mark pseudogenes or low-quality NCBI models.
NCBI_EXCLUDE_PATTERNS = [
    r'\bpseudogene\b',
    r'-like\b',              # e.g. "GAPDH-like"
    r'\bretrotransposed\b',
    r'\bretro\b',
    r'LOW QUALITY PROTEIN',
    r'\bpartial\b',
    r'\bhypothetical\b',
    # 'uncharacterized' deliberately NOT excluded: real genes carry it
]
NCBI_EXCLUDE_REGEX = re.compile('|'.join(NCBI_EXCLUDE_PATTERNS), re.IGNORECASE)


SPECIES_FILE_PATTERNS = {
    'nmr':    ['mHetGla'],
    'batsui': ['BatSui'],
    'amr':    ['DKFZ_Fans'],
}

PROTEIN_EXTS = ['.fa', '.faa', '.fasta', '.pep.fa', '.pep.all.fa', '.aa']
ANNOT_EXTS = ['.gff3', '.gff', '.gtf', '.gff3.gz', '.gff.gz', '.gtf.gz']
CDS_EXTS = ['.cds.fa', '.cds.all.fa', '.cds.fna', '_cds.fa', '.cds.fasta']

# =============================================================================
# LOGGING
# =============================================================================

class Logger:
    """Print to stdout and append to the log file (reopened per message so
    the log survives a crash mid-run)."""

    def __init__(self, path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.path = path
        with open(path, 'w') as f:
            f.write(f"Isoform selection log - {datetime.now():%Y-%m-%d %H:%M:%S}\n\n")

    def log(self, message):
        print(message)
        with open(self.path, 'a') as f:
            f.write(message + "\n")

logger = None

def log(msg):
    if logger:
        logger.log(msg)
    else:
        print(msg)

# =============================================================================
# PARSING HELPERS
# =============================================================================

def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def parse_gff_attributes(attr_string):
    """Parse GFF3 or GTF."""
    attributes = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if not item:
            continue
        if '=' in item:
            key, value = item.split('=', 1)
        elif ' ' in item:
            key, value = item.split(' ', 1)
        else:
            continue
        attributes[key] = value.strip('"')
    return attributes


def extract_transcript_id_from_protein(record, file_format='ensembl'):
    """Transcript ID for a protein record."""
    if file_format in ['ncbi', 'mikado', 'simple', 'braker']:
        rid = record.id
        base = rid.rsplit('.', 1)[0] if '.' in rid else rid
        return rid, base
    # ensembl: transcript ID lives in the description
    match = re.search(r'transcript:(\S+)', record.description)
    if match:
        tid = match.group(1)
        return tid, tid.split('.')[0]
    return record.id, record.id.split('.')[0]


def detect_file_format(protein_file):
    """Classify a protein FASTA from its first 10 headers:
    'ncbi', 'ensembl', 'mikado', 'braker', 'simple', or 'unknown'."""
    for i, record in enumerate(SeqIO.parse(protein_file, 'fasta')):
        rid = record.id
        if rid.startswith(('NP_', 'XP_', 'WP_')):
            return 'ncbi'
        if 'transcript:' in record.description:
            return 'ensembl'
        if rid.startswith('mikado.'):
            return 'mikado'
        if re.match(r'^g\d+\.t\d+$', rid):
            return 'braker'
        if re.match(r'^\w+_T\d+$', rid):
            return 'simple'
        if i >= 9:
            break
    return 'unknown'


def build_transcript_lookup(transcript_set):
    """Index transcript IDs under both their full and version-stripped forms."""
    lookup = {}
    for tid in transcript_set:
        lookup[tid] = tid
        if '.' in tid:
            lookup.setdefault(tid.rsplit('.', 1)[0], tid)
    return lookup


def get_transcript_info_from_gff(gff_file):
    """transcript_id -> {gene_id, CDS length, tags, appris} from a GFF/GTF."""
    log(f"    Parsing annotation: {os.path.basename(gff_file)}")
    transcripts = defaultdict(lambda: {'gene_id': None, 'length': 0, 'tag': [], 'appris': None})
    count = 0
    with open_text(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            ftype = fields[2]
            if ftype not in ['transcript', 'mRNA', 'CDS']:
                continue
            attrs = parse_gff_attributes(fields[8])
            tid = attrs.get('transcript_id') or attrs.get('ID') or attrs.get('Name')
            if not tid:
                continue
            if ftype in ['transcript', 'mRNA']:
                count += 1
                transcripts[tid]['gene_id'] = attrs.get('gene_id') or attrs.get('Parent')
                if 'tag' in attrs:
                    transcripts[tid]['tag'] = attrs['tag'].split(',')
                transcripts[tid]['appris'] = attrs.get('appris')
            else:  # CDS
                transcripts[tid]['length'] += int(fields[4]) - int(fields[3]) + 1
    log(f"    Found {count} transcripts")
    return dict(transcripts)


def get_ncbi_protein_mapping(gff_file):
    """protein_id -> gene_id via the mRNA layer of an NCBI GFF."""
    log(f"    Parsing NCBI annotation: {os.path.basename(gff_file)}")
    t2g, p2t = {}, {}
    with open_text(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            attrs = parse_gff_attributes(fields[8])
            if fields[2] == 'mRNA':
                tid = attrs.get('transcript_id') or attrs.get('ID')
                parent = attrs.get('Parent', '').replace('gene-', '')
                if tid and parent:
                    t2g[tid] = parent
            elif fields[2] == 'CDS':
                pid = attrs.get('protein_id')
                parent = attrs.get('Parent', '')
                if pid and parent:
                    p2t[pid] = parent
    p2g = {}
    for pid, ref in p2t.items():
        gid = t2g.get(ref) or t2g.get(ref.replace('rna-', ''))
        if gid:
            p2g[pid] = gid
    log(f"    Mapped {len(p2g)} proteins to genes")
    return p2g


def parse_simple_gff_for_gene_mapping(gff_file):
    """transcript_id -> gene_id for CAT/BRAKER-style annotations."""
    log(f"    Parsing GFF for gene-transcript mapping: {os.path.basename(gff_file)}")
    t2g = {}
    with open_text(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] not in ['transcript', 'mRNA']:
                continue
            attrs = parse_gff_attributes(fields[8])
            tid = attrs.get('transcript_id') or attrs.get('ID') or attrs.get('Name')
            gid = attrs.get('gene_id') or attrs.get('Parent')
            if tid and gid:
                t2g[tid] = gid.replace('gene-', '').replace('gene:', '')
    log(f"    Mapped {len(t2g)} transcripts to {len(set(t2g.values()))} genes")
    return t2g

# =============================================================================
# DIAMOND
# =============================================================================

def diamond_best_hit_scores(query_fa, db, out_tsv):
    """qseqid -> best-hit bitscore against db, or None if DIAMOND fails."""
    cmd = ["diamond", "blastp", "-q", query_fa, "-d", db, "-o", out_tsv,
           "-p", str(DIAMOND_THREADS), "--outfmt", "6", "qseqid", "bitscore",
           "--max-target-seqs", "1", "--sensitive"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log("    ERROR: DIAMOND blastp failed")
        log(f"    STDERR: {result.stderr}")
        return None
    scores = {}
    if os.path.exists(out_tsv):
        with open(out_tsv) as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2:
                    scores[parts[0]] = float(parts[1])
        os.remove(out_tsv)
    return scores

# =============================================================================
# NCBI PSEUDOGENE FILTERING
# =============================================================================

def is_ncbi_pseudogene(record):
    """(is_bad, reason) for a record matching the exclusion patterns."""
    match = NCBI_EXCLUDE_REGEX.search(record.description)
    if match:
        return True, f"description contains '{match.group()}'"
    match = NCBI_EXCLUDE_REGEX.search(record.id)
    if match:
        return True, f"ID contains '{match.group()}'"
    return False, None


def ncbi_protein_tier(record):
    """3 = NP_ (curated), 2 = XP_ (predicted), 1 = other."""
    if record.id.startswith('NP_'):
        return 3, "NP_ (curated)"
    if record.id.startswith('XP_'):
        return 2, "XP_ (predicted)"
    return 1, "other"


def filter_ncbi_proteins(protein_file, species):
    """Drop pseudogene/low-quality records; log the breakdown with examples."""
    log(f"    Pre-filtering NCBI proteins for {species}...")
    kept = []
    total = 0
    filtered = defaultdict(int)
    examples = defaultdict(list)
    by_tier = defaultdict(int)

    for record in SeqIO.parse(protein_file, 'fasta'):
        total += 1
        is_bad, reason = is_ncbi_pseudogene(record)
        if is_bad:
            filtered[reason] += 1
            if len(examples[reason]) < 3:
                examples[reason].append(record.id)
        else:
            _, tier_name = ncbi_protein_tier(record)
            by_tier[tier_name] += 1
            kept.append(record)

    log(f"    Kept {len(kept)}/{total} proteins "
        f"({100 * len(kept) / total:.1f}%)" if total else "    Empty protein file")
    for reason, count in sorted(filtered.items(), key=lambda x: -x[1]):
        log(f"      filtered - {reason}: {count} (e.g., {', '.join(examples[reason])})")
    for tier, count in sorted(by_tier.items(), reverse=True):
        log(f"      kept    - {tier}: {count}")
    return kept

# =============================================================================
# SELECTION STRATEGIES
# =============================================================================

def pick_cleanest_longest(isoforms):
    """Fewest '.' placeholder characters, then longest sequence."""
    return max(isoforms, key=lambda r: (-str(r.seq).count('.'), len(r.seq)))


def select_isoform_collaborator(species, protein_file):
    """Exact transcript-ID match against the collaborator single-copy TSV."""
    log("    Strategy: Customised TSV (exact match)")
    tsv = COLLABORATOR_TSVS.get(species)
    if not tsv or not os.path.exists(tsv):
        log(f"    ERROR: Customised TSV not found for {species} "
            f"(expected: {tsv or 'no file defined'})")
        return None

    log(f"    Using TSV: {tsv}")
    df = pd.read_csv(tsv, sep='\t', low_memory=False)
    valid_ids = set(df['tid'].dropna().astype(str))
    log(f"    Found {len(valid_ids)} target transcripts in TSV")

    selected = [r for r in SeqIO.parse(protein_file, 'fasta') if r.id in valid_ids]
    log(f"    Matched {len(selected)} proteins")
    return selected


def select_isoform_mane(protein_file, annotation_file):
    """MANE Select transcripts (human)."""
    log("    Strategy: MANE Select")
    info = get_transcript_info_from_gff(annotation_file)
    mane_ids = {t for t, i in info.items()
                if 'MANE_Select' in i['tag'] or 'MANE Select' in i['tag']}
    log(f"    Found {len(mane_ids)} MANE Select transcripts")
    if not mane_ids:
        return None

    lookup = build_transcript_lookup(mane_ids)
    selected = []
    for r in SeqIO.parse(protein_file, 'fasta'):
        tid, base = extract_transcript_id_from_protein(r, 'ensembl')
        if tid in lookup or base in lookup:
            selected.append(r)
    return selected


def select_isoform_ensembl_canonical(protein_file, annotation_file):
    """One transcript per gene by APPRIS principal > Ensembl_canonical > basic
    tag, longest CDS as tie-break."""
    log("    Strategy: APPRIS principal / Ensembl canonical")
    info = get_transcript_info_from_gff(annotation_file)

    gene_groups = defaultdict(list)
    for tid, data in info.items():
        if not data['gene_id']:
            continue
        prio = 0
        if data['appris'] and 'principal' in str(data['appris']):
            prio = 3
        elif 'Ensembl_canonical' in data['tag']:
            prio = 2
        elif 'basic' in data['tag']:
            prio = 1
        gene_groups[data['gene_id']].append((prio, data['length'], tid))

    best_ids = {max(group, key=lambda x: (x[0], x[1]))[2] for group in gene_groups.values()}
    log(f"    Found {len(best_ids)} canonical transcripts")

    lookup = build_transcript_lookup(best_ids)
    fmt = detect_file_format(protein_file)
    selected = []
    for r in SeqIO.parse(protein_file, 'fasta'):
        tid, base = extract_transcript_id_from_protein(r, fmt)
        if tid in lookup or base in lookup:
            selected.append(r)
    return selected


def select_isoform_ncbi_filtered(species, protein_file, annotation_file, human_db):
    """Pre-filter NCBI proteins, then pick one isoform per gene by
    tier (NP_ > XP_) > DIAMOND bitscore to human > length."""
    log("    Strategy: NCBI filtered + DIAMOND best-hit")

    records = filter_ncbi_proteins(protein_file, species)
    if not records:
        log("    ERROR: no proteins remaining after filtering")
        return None

    temp_filtered = f"temp_{species}_filtered.fa"
    SeqIO.write(records, temp_filtered, "fasta")

    p2g = get_ncbi_protein_mapping(annotation_file) if annotation_file else {}
    if not p2g:
        log("    WARNING: no protein-gene mapping - treating each protein as its own gene")

    scores = diamond_best_hit_scores(temp_filtered, human_db, f"{species}_vs_human.tsv")
    os.remove(temp_filtered)
    if scores is None:
        log("    WARNING: continuing with tier + length only (no BLAST scores)")
        scores = {}
    log(f"    DIAMOND hits: {len(scores)} / {len(records)} proteins")

    gene_groups = defaultdict(list)
    for r in records:
        pid = r.id
        base = pid.rsplit('.', 1)[0] if '.' in pid else pid
        gid = p2g.get(pid) or p2g.get(base) or f"gene_{base}"
        tier, _ = ncbi_protein_tier(r)
        gene_groups[gid].append({'record': r, 'tier': tier,
                                 'score': scores.get(pid, 0), 'length': len(r.seq)})
    log(f"    Gene groups: {len(gene_groups)}")

    selected = []
    reasons = {'tier': 0, 'blast': 0, 'length': 0, 'single': 0}
    for isoforms in gene_groups.values():
        if len(isoforms) == 1:
            selected.append(isoforms[0]['record'])
            reasons['single'] += 1
            continue
        ranked = sorted(isoforms, key=lambda x: (x['tier'], x['score'], x['length']),
                        reverse=True)
        best, second = ranked[0], ranked[1]
        if best['tier'] > second['tier']:
            reasons['tier'] += 1
        elif best['score'] > second['score']:
            reasons['blast'] += 1
        else:
            reasons['length'] += 1
        selected.append(best['record'])

    log(f"    Selection breakdown: single={reasons['single']}, tier={reasons['tier']}, "
        f"blast={reasons['blast']}, length={reasons['length']}")
    log(f"    Final selection: {len(selected)} proteins")
    return selected


def select_isoform_simple_format(protein_file, human_db=None, annotation_file=None):
    """BRAKER/mikado/simple formats: group isoforms per gene via the GFF
    (or the g<N>.t<M> ID pattern) and pick the best by DIAMOND best-hit to
    human, penalizing '.' placeholder characters; cleanest/longest fallback."""
    log("    Strategy: simple format with GFF-based gene grouping")

    records = list(SeqIO.parse(protein_file, 'fasta'))
    log(f"    Total proteins: {len(records)}")

    dot_counts = {r.id: str(r.seq).count('.') for r in records}
    seqs_with_dots = sum(1 for c in dot_counts.values() if c > 0)
    if seqs_with_dots:
        log(f"    WARNING: {seqs_with_dots} sequences contain '.' "
            f"({100 * seqs_with_dots / len(records):.1f}%) - cleaned for DIAMOND "
            f"('.' -> X) and penalized in selection")

    t2g = {}
    if annotation_file and os.path.exists(annotation_file):
        t2g = parse_simple_gff_for_gene_mapping(annotation_file)

    gene_groups = defaultdict(list)
    for r in records:
        gid = t2g.get(r.id)
        if not gid:
            if re.match(r'^g\d+\.t\d+$', r.id):        # BRAKER: g1.t1, g1.t2 -> g1
                gid = r.id.rsplit('.', 1)[0]
            else:                                       # no mapping: treat as unique
                gid = r.id
        gene_groups[gid].append(r)

    multi_isoform = sum(1 for g in gene_groups.values() if len(g) > 1)
    log(f"    Gene groups: {len(gene_groups)} ({multi_isoform} with multiple isoforms)")
    if multi_isoform == 0:
        log("    No isoform selection needed - keeping all proteins")
        return records

    if not human_db or not os.path.exists(human_db):
        log(f"    ERROR: human DIAMOND database not found: {human_db}")
        log("    WARNING: falling back to cleanest/longest isoform")
        return [pick_cleanest_longest(g) for g in gene_groups.values()]

    temp_cleaned = "temp_cleaned_for_diamond.fa"
    with open(temp_cleaned, 'w') as f:
        for r in records:
            f.write(f">{r.id}\n{str(r.seq).replace('.', 'X').rstrip('*')}\n")

    scores = diamond_best_hit_scores(temp_cleaned, human_db, "temp_simple_vs_human.tsv")
    os.remove(temp_cleaned)
    if scores is None:
        log("    WARNING: BLAST failed - falling back to cleanest/longest isoform")
        return [pick_cleanest_longest(g) for g in gene_groups.values()]

    log(f"    DIAMOND hits: {len(scores)} / {len(records)} proteins")
    hit_rate = len(scores) / len(records) if records else 0
    if hit_rate < 0.5 and len(records) > 100:
        log(f"    WARNING: low DIAMOND hit rate ({hit_rate * 100:.1f}%)")
    if not scores:
        log("    WARNING: DIAMOND returned 0 hits - selection falls back to cleanest/longest")

    selected = []
    blast_selected = length_selected = 0
    for isoforms in gene_groups.values():
        if len(isoforms) == 1:
            selected.append(isoforms[0])
            continue
        with_hits = [r for r in isoforms if r.id in scores]
        if with_hits:
            # each '.' placeholder costs ~1 mismatch worth of bitscore
            best = max(with_hits, key=lambda r: scores[r.id] - dot_counts[r.id] * 10)
            blast_selected += 1
        else:
            best = pick_cleanest_longest(isoforms)
            length_selected += 1
        selected.append(best)

    log(f"    Multi-isoform genes resolved by BLAST: {blast_selected}, "
        f"by cleanest/longest: {length_selected}")
    selected_with_dots = sum(1 for r in selected if dot_counts.get(r.id, 0) > 0)
    if selected_with_dots:
        log(f"    WARNING: {selected_with_dots} selected isoforms contain '.'")
    log(f"    Selected {len(selected)} proteins (one per gene)")
    return selected


def keep_all_proteins(protein_file):
    log("    Strategy: keeping ALL proteins (fallback)")
    return list(SeqIO.parse(protein_file, 'fasta'))

# =============================================================================
# INPUT FILE LOCATION
# =============================================================================

def find_input_file(directory, exts, species, version, include_generic=()):
    """First file in directory matching *pattern*, trying the version-derived
    patterns first, then SPECIES_FILE_PATTERNS, then the species name."""
    files = []
    for ext in exts:
        files.extend(glob.glob(os.path.join(directory, f"*{ext}")))
    for ext in include_generic:
        for f in glob.glob(os.path.join(directory, f"*{ext}")):
            if f not in files:
                files.append(f)
    if not files:
        return None

    patterns = []
    if version:
        patterns.append(version)
        if '.' in version:
            parts = version.split('.')
            patterns.append(parts[0])
            patterns.append('.'.join(parts[:2]))
    patterns.extend(SPECIES_FILE_PATTERNS.get(species.lower(), []))
    patterns.extend([species, species.capitalize()])

    for pattern in patterns:
        for f in files:
            if pattern.lower() in os.path.basename(f).lower():
                return f
    return None

# =============================================================================
# MAIN
# =============================================================================

def main():
    global logger
    logger = Logger("logs/01_isoform_selection/01_isoform_selection_v10.log")

    log(f"Isoform selection - NCBI species with strict filtering: {NCBI_SPECIES}")
    log(f"NCBI exclusion patterns: {NCBI_EXCLUDE_PATTERNS}\n")

    if not shutil.which("diamond"):
        log("ERROR: 'diamond' not found on PATH")
        sys.exit(1)
    if not os.path.exists('species_list_17missinggene.csv'):
        log("ERROR: species_list_17missinggene.csv not found")
        sys.exit(1)

    species_df = pd.read_csv('species_list_17missinggene.csv')
    os.makedirs('results/01_isoform_selection', exist_ok=True)

    # --- human DIAMOND reference database ---
    log("--- Initializing human reference ---")
    human_row = species_df[species_df['Short_name'].str.strip() == 'human'].iloc[0]
    human_prot = find_input_file('proteins', PROTEIN_EXTS, 'human',
                                 human_row['Version'].strip())
    if not human_prot:
        log("ERROR: human protein file not found")
        sys.exit(1)
    log(f"  Human protein file: {os.path.basename(human_prot)}")

    human_db = "human_mane_db.dmnd"
    result = subprocess.run(["diamond", "makedb", "--in", human_prot, "-d", "human_mane_db"],
                            capture_output=True, text=True)
    if result.returncode != 0 or not os.path.exists(human_db) or os.path.getsize(human_db) == 0:
        log("ERROR: failed to build human DIAMOND database")
        log(f"STDERR: {result.stderr}")
        sys.exit(1)
    log(f"  Human DIAMOND database created ({os.path.getsize(human_db) / 1024 / 1024:.1f} MB)\n")

    # --- process every species ---
    summary = []
    for _, row in species_df.iterrows():
        species = row['Short_name'].strip().lower()
        version = row['Version'].strip()

        log(f"Processing: {species} (version: {version})")

        p_file = find_input_file('proteins', PROTEIN_EXTS, species, version)
        a_file = find_input_file('annotations', ANNOT_EXTS, species, version)
        c_file = find_input_file('cds', CDS_EXTS, species, version,
                                 include_generic=('.fa', '.fna'))  # logged for provenance

        if not p_file:
            log(f"  Protein file not found for '{species}' "
                f"(patterns: {SPECIES_FILE_PATTERNS.get(species, [species])}) - SKIPPING\n")
            summary.append({'species': species, 'count': 0, 'format': '-',
                            'strategy': 'SKIPPED (no protein file)'})
            continue

        log(f"  Protein file:    {os.path.basename(p_file)}")
        log(f"  Annotation file: {os.path.basename(a_file) if a_file else 'NOT FOUND'}")
        log(f"  CDS file:        {os.path.basename(c_file) if c_file else 'NOT FOUND'}")

        fmt = detect_file_format(p_file)
        log(f"  Detected format: {fmt}")

        try:
            selected, strategy = None, ""

            if species in ['nmr', 'batsui']:
                selected = select_isoform_collaborator(species, p_file)
                strategy = "Customised TSV"
                if selected is None:
                    log("    Falling back to simple format strategy...")
                    selected = select_isoform_simple_format(p_file, human_db)
                    strategy = "Simple format (TSV fallback)"

            elif species == 'human':
                if a_file:
                    selected = select_isoform_mane(p_file, a_file)
                    strategy = "MANE Select"

            elif species in NCBI_SPECIES or fmt == 'ncbi':
                log("  NCBI species - applying strict pseudogene filtering")
                selected = select_isoform_ncbi_filtered(species, p_file, a_file, human_db)
                strategy = "NCBI Filtered + Best-Hit"

            elif fmt == 'ensembl':
                if a_file:
                    selected = select_isoform_ensembl_canonical(p_file, a_file)
                    strategy = "Ensembl Canonical"
                else:
                    log("  WARNING: Ensembl format but no annotation - keeping all")
                    selected = keep_all_proteins(p_file)
                    strategy = "ALL (no annotation)"

            elif fmt in ['simple', 'braker', 'mikado', 'unknown']:
                selected = select_isoform_simple_format(p_file, human_db, a_file)
                strategy = f"Simple format ({fmt})"

            if not selected:
                selected = keep_all_proteins(p_file)
                strategy = "ALL (Fallback)"

            out = os.path.join('results/01_isoform_selection', f"{species}.primary.fa")
            SeqIO.write(selected, out, "fasta")
            log(f"  Selected {len(selected)} isoforms ({strategy}) -> {out}\n")
            summary.append({'species': species, 'count': len(selected),
                            'format': fmt, 'strategy': strategy})

        except Exception as e:
            log(f"  ERROR while processing {species}: {e}")
            log(traceback.format_exc())
            summary.append({'species': species, 'count': 0, 'format': fmt,
                            'strategy': f'FAILED ({e})'})

    if os.path.exists(human_db):
        os.remove(human_db)

    # --- summary ---
    log("SUMMARY")
    log(f"{'Species':<15} {'Count':<10} {'Format':<12} {'Strategy'}")
    log("-" * 70)
    for r in summary:
        log(f"{r['species']:<15} {r['count']:<10} {r['format']:<12} {r['strategy']}")

    problems = [r for r in summary if r['strategy'].startswith(('SKIPPED', 'FAILED'))]
    if problems:
        log(f"\nWARNING: {len(problems)} species produced no output: "
            + ", ".join(r['species'] for r in problems))
    log(f"\nDone. Log: {logger.path}")


if __name__ == "__main__":
    main()
