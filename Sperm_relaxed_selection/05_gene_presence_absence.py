#!/usr/bin/env python3
"""
05_gene_presence_absence.py

Generates a species-by-gene presence/absence matrix from the final
ClipKIT-trimmed alignments that were used as RELAX input.

This provides a definitive record of which species retained each gene
through all pipeline filtering steps (OrthoFinder -> tree pruning ->
CDS extraction -> MACSE alignment -> ClipKIT trimming -> gap filtering).

Also integrates degradation flags (frameshifts, internal stop codons)
from the focal alignment analysis if available.

Output:
  - presence_absence_matrix.tsv   Gene x Species matrix (1/0)
  - presence_absence_detailed.tsv Long-format with sequence lengths
  - presence_absence_summary.txt  Console-formatted report

Usage:
  python3 05_gene_presence_absence.py
  python3 05_gene_presence_absence_v4.py --alignment_dir results/08_clipkit/sperm/clipkit_trimmed
  python3 05_gene_presence_absence.py --degradation results/10_focal_alignments/degradation_flags.tsv
"""

import argparse
import os
import re
import sys
from collections import OrderedDict

# =============================================================================
# CONFIGURATION
# =============================================================================

TARGET_GENES = [
    "ADIG", "ATP1A4", "CABS1", "DNHD1", "DPY19L2",
    "GARIN3", "GARIN1A", "FSIP2", "MROH2B", "NDC80",   # were FAM71B, FAM71F2
    "NPM2", "RIMBP3", "SPATA33", "SPEM3", "TMPRSS12",
    "TTLL8", "ZAN",
]

# All 23 species in phylogenetic display order
SPECIES_ORDER = [
    "NMR", "BATS", "AMR",           # Bathyergidae
    "GPIG", "CHIN", "DEGU",         # Hystricomorpha
    "ACAH",                          # Acomys
    "MOUS", "RAT", "WWMO",          # Muridae
    "PVOL", "NADE",                  # Other Myomorpha
    "SQUI",                          # Sciuridae
    "RABB",                          # Lagomorpha
    "DOG", "CAT",                    # Carnivora
    "PIG", "BOS",                    # Artiodactyla
    "HUMA", "BONO", "GORI", "MACA", # Primates
    "ELEP",                          # Afrotheria
]

SPECIES_LABELS = {
    "NMR":  "HetGla",
    "BATS": "BatSui",
    "AMR":  "FukAns",
    "GPIG": "Guinea pig",
    "CHIN": "Chinchilla",
    "DEGU": "Degu",
    "ACAH": "Spiny mouse",
    "MOUS": "Mouse",
    "RAT":  "Rat",
    "WWMO": "W. wild mouse",
    "PVOL": "Prairie vole",
    "NADE": "Deer mouse",
    "SQUI": "Red squirrel",
    "RABB": "Rabbit",
    "DOG":  "Dog",
    "CAT":  "Cat",
    "PIG":  "Pig",
    "BOS":  "Cattle",
    "HUMA": "Human",
    "BONO": "Bonobo",
    "GORI": "Gorilla",
    "MACA": "Macaque",
    "ELEP": "Elephant",
}

SPECIES_GROUPS = {
    "NMR": "Bathyergidae", "BATS": "Bathyergidae", "AMR": "Bathyergidae",
    "GPIG": "Hystricomorpha", "CHIN": "Hystricomorpha", "DEGU": "Hystricomorpha",
    "ACAH": "Acomys",
    "MOUS": "Muridae", "RAT": "Muridae", "WWMO": "Muridae",
    "PVOL": "Myomorpha", "NADE": "Myomorpha",
    "SQUI": "Sciuridae", "RABB": "Lagomorpha",
    "DOG": "Carnivora", "CAT": "Carnivora",
    "PIG": "Artiodactyla", "BOS": "Artiodactyla",
    "HUMA": "Primates", "BONO": "Primates", "GORI": "Primates", "MACA": "Primates",
    "ELEP": "Afrotheria",
}


# =============================================================================
# FUNCTIONS
# =============================================================================

def extract_species_from_header(header):
    """Extract species code from alignment sequence header.

    Handles main pipeline formats (Ensembl/mikado IDs):
      >ENSMUST00000123456        (Ensembl transcript)
      >ENSBTAT00000147632_1      (Ensembl with copy suffix)
      >mikado.chr1.23.t1         (mikado for NMR)
      >mikado.ptg000123.23.t1    (mikado for BATS/BatSui)
      >g123.t1                   (FukAns/AMR)
      >lcl|NC_012345.1_cds_...   (NCBI for Bonobo)

    Also handles control pipeline format:
      >NMR_ADIG                  (species_gene)
    """
    seq_id = header.split()[0].lstrip('>')

    # --- Ensembl prefix mapping (transcript ID prefix -> species code) ---
    ENSEMBL_PREFIX_MAP = {
        'ENSMUST': 'MOUS',   # Mus musculus
        'ENSRNOT': 'RAT',    # Rattus norvegicus
        'ENST0':   'HUMA',   # Homo sapiens
        'ENSCPOT': 'GPIG',   # Cavia porcellus
        'ENSODET': 'DEGU',   # Octodon degus
        'ENSCLAT': 'CHIN',   # Chinchilla lanigera
        'ENSACHT': 'ACAH',   # Acomys cahirinus
        'ENSMOCT': 'PVOL',   # Microtus ochrogaster
        'ENSMSPT': 'WWMO',   # Mus spretus
        'ENSSVLT': 'SQUI',   # Sciurus vulgaris
        'ENSMMUT': 'MACA',   # Macaca mulatta
        'ENSPEMT': 'NADE',   # Peromyscus maniculatus
        'ENSBTAT': 'BOS',    # Bos taurus
        'ENSSSCT': 'PIG',    # Sus scrofa
        'ENSCAFT': 'DOG',    # Canis lupus familiaris
        'ENSFCAT': 'CAT',    # Felis catus (Ensembl standard)
        'ENSFCTT': 'CAT',    # Felis catus (alternate prefix)
        'ENSLAFT': 'ELEP',   # Loxodonta africana
        'ENSOCUT': 'RABB',   # Oryctolagus cuniculus
        'ENSGGOT': 'GORI',   # Gorilla gorilla
        'ENSPPAT': 'BONO',   # Pan paniscus (if Ensembl)
    }

    # Check Ensembl prefixes (longest match first to avoid ambiguity)
    for prefix, species in sorted(ENSEMBL_PREFIX_MAP.items(),
                                   key=lambda x: -len(x[0])):
        if seq_id.startswith(prefix):
            return species

    # Mikado IDs for bathyergid species
    # Handles both dot and underscore formats:
    #   mikado.chr1.23.t1  or  mikado_chr15G2025_1
    if seq_id.startswith('mikado_chr') or seq_id.startswith('mikado.chr') or \
       seq_id.startswith('mikado_Chr') or seq_id.startswith('mikado.Chr'):
        return 'NMR'    # HetGla (naked mole-rat)
    if seq_id.startswith('mikado_ptg') or seq_id.startswith('mikado.ptg') or \
       seq_id.startswith('mikado_Ptg') or seq_id.startswith('mikado.Ptg'):
        return 'BATS'   # BatSui (cape dune mole-rat, solitary)

    # FukAns (Ansell's mole-rat, eusocial) uses g\d+.t\d+ or g\d+_t\d+ pattern
    if re.match(r'^g\d+[._]t\d+', seq_id):
        return 'AMR'

    # NCBI format for Bonobo - handles both pipe and underscore separators:
    #   lcl|NC_012345...  or  lcl_NC_012345...
    if seq_id.startswith('lcl|NC_') or seq_id.startswith('lcl_NC_') or \
       seq_id.startswith('lcl|NW_') or seq_id.startswith('lcl_NW_'):
        return 'BONO'

    # Control pipeline format: SPECIES_GENE
    parts = seq_id.split('_')
    if parts[0] in SPECIES_LABELS:
        return parts[0]

    return None


def parse_alignment(filepath):
    """Parse a FASTA alignment, returning {species: sequence_length} for
    the ungapped sequence length, plus the total alignment length."""
    species_data = {}
    current_id = None
    current_seq = []
    aln_len = None

    with open(filepath) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id is not None:
                    seq = ''.join(current_seq)
                    if aln_len is None:
                        aln_len = len(seq)
                    sp = extract_species_from_header(current_id)
                    if sp:
                        ungapped = len(seq.replace('-', '').replace('!', ''))
                        gap_pct = ((len(seq) - ungapped) / len(seq) * 100
                                   if len(seq) > 0 else 0)
                        species_data[sp] = {
                            'seq_len': ungapped,
                            'aln_len': len(seq),
                            'gap_pct': gap_pct,
                        }
                current_id = line
                current_seq = []
            else:
                current_seq.append(line)

        # Last sequence
        if current_id is not None:
            seq = ''.join(current_seq)
            if aln_len is None:
                aln_len = len(seq)
            sp = extract_species_from_header(current_id)
            if sp:
                ungapped = len(seq.replace('-', '').replace('!', ''))
                gap_pct = ((len(seq) - ungapped) / len(seq) * 100
                           if len(seq) > 0 else 0)
                species_data[sp] = {
                    'seq_len': ungapped,
                    'aln_len': len(seq),
                    'gap_pct': gap_pct,
                }

    return species_data, aln_len


# legacy display symbol -> symbol actually used in the v5 filenames
FILE_ALIASES = {"FAM71B": "GARIN3", "FAM71F2": "GARIN1A"}


def find_gene_files(alignment_dir):
    """Find alignment files for each target gene.

    Handles naming patterns:
      GENE_NT.clipkit.fa
      GENE_NT_clipkit.fasta
      GENE_trimmed.fa
    """
    gene_files = {}
    all_files = os.listdir(alignment_dir)

    for gene in TARGET_GENES:
        # Try common patterns (files carry the current HGNC symbol)
        stem = FILE_ALIASES.get(gene, gene)
        candidates = [f for f in all_files
                      if f.startswith(stem + '_') and
                      any(f.endswith(ext) for ext in ['.fa', '.fasta', '.fna'])]
        if candidates:
            # Prefer NT (nucleotide) over AA
            nt_files = [f for f in candidates if 'NT' in f or 'nt' in f]
            if nt_files:
                gene_files[gene] = os.path.join(alignment_dir, nt_files[0])
            else:
                gene_files[gene] = os.path.join(alignment_dir, candidates[0])

    return gene_files


def load_degradation_flags(deg_file):
    """Load degradation flags TSV if available.

    Returns dict: {(gene, species): {'frameshifts': N, 'stops': N}}
    """
    flags = {}
    if not deg_file or not os.path.exists(deg_file):
        return flags

    with open(deg_file) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                gene = parts[0]
                species = parts[1]
                n_fs = int(parts[3]) if parts[3].isdigit() else 0
                n_stops = int(parts[4]) if parts[4].isdigit() else 0
                flags[(gene, species)] = {
                    'frameshifts': n_fs,
                    'stops': n_stops,
                }

    return flags


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate gene x species presence/absence matrix"
    )
    parser.add_argument("--alignment_dir", default="results/08_clipkit/sperm/clipkit_trimmed",
                        help="Directory with final trimmed alignments")
    parser.add_argument("--output_dir", default="results/12_presence_absence",
                        help="Output directory")
    parser.add_argument("--degradation", default=None,
                        help="Path to degradation_flags.tsv (optional)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 75)
    print("GENE x SPECIES PRESENCE/ABSENCE MATRIX")
    print("=" * 75)
    print(f"Alignment dir:  {args.alignment_dir}")
    print(f"Output dir:     {args.output_dir}")
    print(f"Target genes:   {len(TARGET_GENES)}")
    print(f"Species:        {len(SPECIES_ORDER)}")
    if args.degradation:
        print(f"Degradation:    {args.degradation}")
    print("=" * 75)

    # Find alignment files
    gene_files = find_gene_files(args.alignment_dir)
    print(f"\nFound alignments for {len(gene_files)}/{len(TARGET_GENES)} genes")

    missing_genes = [g for g in TARGET_GENES if g not in gene_files]
    if missing_genes:
        print(f"  Missing: {', '.join(missing_genes)}")

    # Load degradation flags
    deg_flags = load_degradation_flags(args.degradation)
    if deg_flags:
        print(f"Loaded {len(deg_flags)} degradation flags")

    # Parse all alignments
    print("\nParsing alignments...")
    all_data = {}  # {gene: {species: {seq_len, aln_len, gap_pct}}}
    unrecognized_ids = set()

    for gene in TARGET_GENES:
        if gene not in gene_files:
            all_data[gene] = {}
            continue

        species_data, aln_len = parse_alignment(gene_files[gene])
        all_data[gene] = species_data

        # Check for unrecognized sequences
        with open(gene_files[gene]) as f:
            for line in f:
                if line.startswith('>'):
                    seq_id = line.strip().lstrip('>').split()[0]
                    sp = extract_species_from_header(line.strip())
                    if sp is None:
                        unrecognized_ids.add(seq_id)

        n_sp = len(species_data)
        print(f"  {gene:<12s} {n_sp:2d} species  aln_len={aln_len} bp")

    if unrecognized_ids:
        print(f"\n  WARNING: {len(unrecognized_ids)} unrecognized sequence IDs:")
        for uid in sorted(unrecognized_ids):
            print(f"    {uid}")
        print("  These sequences were not mapped to any species.")

    # =========================================================================
    # OUTPUT 1: Presence/absence matrix (wide format)
    # =========================================================================
    matrix_file = os.path.join(args.output_dir, "presence_absence_matrix.tsv")
    with open(matrix_file, 'w') as f:
        # Header
        f.write("Gene\t" + "\t".join(SPECIES_ORDER) + "\tTotal\n")

        for gene in TARGET_GENES:
            row = [gene]
            count = 0
            for sp in SPECIES_ORDER:
                present = 1 if sp in all_data.get(gene, {}) else 0
                row.append(str(present))
                count += present
            row.append(str(count))
            f.write("\t".join(row) + "\n")

        # Species totals row
        row = ["TOTAL"]
        for sp in SPECIES_ORDER:
            total = sum(1 for gene in TARGET_GENES if sp in all_data.get(gene, {}))
            row.append(str(total))
        row.append("")
        f.write("\t".join(row) + "\n")

    print(f"\nSaved: {matrix_file}")

    # =========================================================================
    # OUTPUT 2: Detailed long-format table
    # =========================================================================
    detail_file = os.path.join(args.output_dir, "presence_absence_detailed.tsv")
    with open(detail_file, 'w') as f:
        f.write("Gene\tSpecies\tLabel\tGroup\tPresent\tSeq_len_bp\t"
                "Aln_len_bp\tGap_pct\tFrameshifts\tInternal_stops\tStatus\n")

        for gene in TARGET_GENES:
            for sp in SPECIES_ORDER:
                label = SPECIES_LABELS.get(sp, sp)
                group = SPECIES_GROUPS.get(sp, "Other")
                sp_data = all_data.get(gene, {}).get(sp)

                if sp_data:
                    present = 1
                    seq_len = sp_data['seq_len']
                    aln_len = sp_data['aln_len']
                    gap_pct = f"{sp_data['gap_pct']:.1f}"
                else:
                    present = 0
                    seq_len = 0
                    aln_len = 0
                    gap_pct = ""

                # Degradation info
                deg = deg_flags.get((gene, sp), {})
                n_fs = deg.get('frameshifts', '')
                n_stops = deg.get('stops', '')

                # Status classification
                if not present:
                    status = "absent"
                elif n_fs or n_stops:
                    status = "degraded"
                else:
                    status = "intact"

                f.write(f"{gene}\t{sp}\t{label}\t{group}\t{present}\t"
                        f"{seq_len}\t{aln_len}\t{gap_pct}\t"
                        f"{n_fs}\t{n_stops}\t{status}\n")

    print(f"Saved: {detail_file}")

    # =========================================================================
    # OUTPUT 3: Heatmap figure
    # =========================================================================
    print("\nGenerating heatmap figure...")

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np

        # Build matrix (genes x species)
        matrix = []
        gene_labels = []
        for gene in TARGET_GENES:
            row = []
            for sp in SPECIES_ORDER:
                if sp in all_data.get(gene, {}):
                    # Check degradation
                    deg = deg_flags.get((gene, sp), {})
                    if deg.get('frameshifts', 0) > 0 or deg.get('stops', 0) > 0:
                        row.append(0.5)   # Degraded
                    else:
                        row.append(1.0)   # Present
                else:
                    row.append(0.0)       # Absent
            matrix.append(row)
            gene_labels.append(gene)

        matrix = np.array(matrix)

        # Sort genes by number of species present (fewest first = most losses at top)
        gene_counts_arr = matrix.sum(axis=1)
        sort_idx = np.argsort(gene_counts_arr)
        matrix = matrix[sort_idx]
        gene_labels = [gene_labels[i] for i in sort_idx]

        # Custom colormap: absent=coral, degraded=gold, present=steelblue
        from matplotlib.colors import ListedColormap, BoundaryNorm
        cmap = ListedColormap(['#F4847B', '#F5C542', '#4A90D9'])
        bounds = [-0.25, 0.25, 0.75, 1.25]
        norm = BoundaryNorm(bounds, cmap.N)

        # Figure dimensions
        n_genes = len(gene_labels)
        n_species = len(SPECIES_ORDER)
        fig_w = max(10, n_species * 0.45 + 3)
        fig_h = max(5, n_genes * 0.4 + 2.5)

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        im = ax.imshow(matrix, aspect='auto', cmap=cmap, norm=norm,
                        interpolation='nearest')

        # Species labels (x-axis)
        sp_labels = [SPECIES_LABELS.get(sp, sp) for sp in SPECIES_ORDER]
        ax.set_xticks(range(n_species))
        ax.set_xticklabels(sp_labels, rotation=55, ha='right', fontsize=11)

        # Gene labels (y-axis)
        ax.set_yticks(range(n_genes))
        ax.set_yticklabels(gene_labels, fontsize=9, fontweight='bold')

        # Color x-axis labels by taxonomic group
        group_colors = {
            'Bathyergidae': '#D62728', 'Hystricomorpha': '#FF7F0E',
            'Acomys': '#8C564B', 'Muridae': '#2CA02C',
            'Myomorpha': '#17BECF', 'Sciuridae': '#9467BD',
            'Lagomorpha': '#BCBD22', 'Carnivora': '#E377C2',
            'Artiodactyla': '#7F7F7F', 'Primates': '#1F77B4',
            'Afrotheria': '#D4A017',
        }
        for i, sp in enumerate(SPECIES_ORDER):
            grp = SPECIES_GROUPS.get(sp, 'Other')
            color = group_colors.get(grp, 'black')
            ax.get_xticklabels()[i].set_color(color)
            # Bold for Bathyergidae
            if grp == 'Bathyergidae':
                ax.get_xticklabels()[i].set_fontweight('bold')

        # Group separator lines (vertical)
        prev_group = None
        for i, sp in enumerate(SPECIES_ORDER):
            grp = SPECIES_GROUPS.get(sp, 'Other')
            if prev_group is not None and grp != prev_group:
                ax.axvline(x=i - 0.5, color='white', linewidth=2)
            prev_group = grp

        # Grid lines
        for i in range(n_genes + 1):
            ax.axhline(y=i - 0.5, color='white', linewidth=0.5)
        for i in range(n_species + 1):
            ax.axvline(x=i - 0.5, color='white', linewidth=0.5)

        # Title
        n_present = int(matrix[matrix >= 0.5].sum())
        n_absent = int((matrix < 0.25).sum())
        ax.set_title(f'Gene Presence/Absence Across {n_species} Mammalian Species',
                    
                     fontsize=12, fontweight='bold', pad=15)

        # Legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#4A90D9', edgecolor='white', label='Present'),
            Patch(facecolor='#F4847B', edgecolor='white', label='Absent'),
        ]
        if deg_flags:
            legend_elements.insert(1,
                Patch(facecolor='#F5C542', edgecolor='white', label='Degraded'))

        ax.legend(handles=legend_elements, loc='upper left',
                  bbox_to_anchor=(1.01, 1.0), ncol=1, fontsize=9,
                  frameon=True, edgecolor='grey')

        plt.tight_layout()

        heatmap_file = os.path.join(args.output_dir, 'presence_absence_heatmap')
        plt.savefig(heatmap_file + '.png', dpi=600, bbox_inches='tight',
                    facecolor='white')
        plt.savefig(heatmap_file + '.pdf', bbox_inches='tight',
                    facecolor='white')
        plt.close()
        print(f"Saved: {heatmap_file}.png / .pdf")

    except ImportError:
        print("  matplotlib not available, skipping heatmap")

    # =========================================================================
    # OUTPUT 4: Console report + text file
    # =========================================================================
    lines = []
    lines.append("=" * 90)
    lines.append("GENE x SPECIES PRESENCE/ABSENCE MATRIX")
    lines.append(f"Source: {args.alignment_dir}")
    lines.append("=" * 90)
    lines.append("")

    # Header row
    sp_short = [sp[:4] for sp in SPECIES_ORDER]
    hdr = f"{'Gene':<12s} |"
    for sp in SPECIES_ORDER:
        hdr += f" {sp:>4s}"
    hdr += " | N"
    lines.append(hdr)
    lines.append("-" * len(hdr))

    gene_counts = {}
    species_counts = {sp: 0 for sp in SPECIES_ORDER}

    for gene in TARGET_GENES:
        row = f"{gene:<12s} |"
        count = 0
        for sp in SPECIES_ORDER:
            present = sp in all_data.get(gene, {})
            if present:
                # Check degradation
                deg = deg_flags.get((gene, sp), {})
                if deg.get('frameshifts', 0) > 0 or deg.get('stops', 0) > 0:
                    row += "   D"  # Degraded
                else:
                    row += "   +"  # Present and intact
                count += 1
                species_counts[sp] += 1
            else:
                row += "   -"  # Absent
        row += f" | {count}"
        gene_counts[gene] = count
        lines.append(row)

    # Totals row
    lines.append("-" * len(hdr))
    row = f"{'TOTAL':<12s} |"
    for sp in SPECIES_ORDER:
        row += f" {species_counts[sp]:>4d}"
    row += " |"
    lines.append(row)

    lines.append("")
    lines.append("Legend:  + = present   - = absent   D = present but degraded")
    lines.append("")

    # Summary statistics
    lines.append("=" * 60)
    lines.append("SUMMARY STATISTICS")
    lines.append("=" * 60)

    n_present = sum(gene_counts.values())
    n_total = len(TARGET_GENES) * len(SPECIES_ORDER)
    lines.append(f"Total cells:     {n_total}")
    lines.append(f"Present:         {n_present} ({100*n_present/n_total:.1f}%)")
    lines.append(f"Absent:          {n_total - n_present} ({100*(n_total-n_present)/n_total:.1f}%)")
    lines.append("")

    # Genes with fewest species
    lines.append("Genes with fewest species:")
    for gene, count in sorted(gene_counts.items(), key=lambda x: x[1]):
        missing_sp = [sp for sp in SPECIES_ORDER if sp not in all_data.get(gene, {})]
        if missing_sp:
            lines.append(f"  {gene:<12s} {count:2d}/{len(SPECIES_ORDER)} species  "
                         f"missing: {', '.join(missing_sp)}")

    lines.append("")

    # Species with fewest genes
    lines.append("Species with fewest genes:")
    for sp, count in sorted(species_counts.items(), key=lambda x: x[1]):
        if count < len(TARGET_GENES):
            missing_genes = [g for g in TARGET_GENES if sp not in all_data.get(g, {})]
            lines.append(f"  {sp:<5s} ({SPECIES_LABELS.get(sp, sp):20s}) "
                         f"{count:2d}/{len(TARGET_GENES)} genes  "
                         f"missing: {', '.join(missing_genes)}")

    lines.append("")

    # Bathyergidae focus
    lines.append("=" * 60)
    lines.append("BATHYERGIDAE FOCUS")
    lines.append("=" * 60)
    bath_species = ["NMR", "BATS", "AMR"]
    for sp in bath_species:
        present_genes = [g for g in TARGET_GENES if sp in all_data.get(g, {})]
        absent_genes = [g for g in TARGET_GENES if sp not in all_data.get(g, {})]
        lines.append(f"  {sp} ({SPECIES_LABELS[sp]}):")
        lines.append(f"    Present: {len(present_genes)}/{len(TARGET_GENES)}")
        if absent_genes:
            lines.append(f"    Absent:  {', '.join(absent_genes)}")
        # Degraded
        degraded = [g for g in present_genes
                    if deg_flags.get((g, sp), {}).get('frameshifts', 0) > 0
                    or deg_flags.get((g, sp), {}).get('stops', 0) > 0]
        if degraded:
            lines.append(f"    Degraded: {', '.join(degraded)}")
    lines.append("")

    # Print to console
    report_text = "\n".join(lines)
    print("\n" + report_text)

    # Save report
    report_file = os.path.join(args.output_dir, "presence_absence_summary.txt")
    with open(report_file, 'w') as f:
        f.write(report_text)
    print(f"Saved: {report_file}")

    print(f"\n{'=' * 60}")
    print("DONE")
    print(f"{'=' * 60}")
    print(f"Output: {args.output_dir}/")
    print(f"  presence_absence_matrix.tsv    Wide-format matrix")
    print(f"  presence_absence_detailed.tsv  Long-format with lengths")
    print(f"  presence_absence_heatmap.png   Heatmap figure")
    print(f"  presence_absence_summary.txt   Text report")


if __name__ == "__main__":
    main()
