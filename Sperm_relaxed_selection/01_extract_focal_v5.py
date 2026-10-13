#!/usr/bin/env python3
"""
01_extract_focal_v5.py

Extract focal + outgroup species from MACSE codon-aware NT alignments for
protein-level disruption visualization (unified pipeline).

Reads MACSE raw alignments (results/.../macse_aligned/{GENE}_NT.aln) rather
than ClipKIT-trimmed ones, because MACSE preserves the internal stop codons,
frameshifts (marked '!'), and degraded regions that the disruption plot is
meant to SHOW -- ClipKIT trimming + stop-removal would erase exactly that signal.

v5 notes (v4 notes retained below):
  - MACSE headers are BARE SPECIES LABELS (NMR, AMR, BATS, HUMA, ...), so
    species selection is a simple label filter -- no Ensembl/mikado ID prefix
    matching needed.
  - Records ABSENCE explicitly: if a chosen species is not in a gene's
    alignment (e.g. NMR lost in CABS1), it is written to a companion
    presence/absence table so the R plot can draw the pink "gene absent" bar.

Output per gene:
  {output_dir}/{GENE}_focal_NT.aln     the sub-alignment (chosen species, ordered)
  {output_dir}/{GENE}_presence.tsv     species | present(0/1)   (for absence bars)

Usage:
  python3 01_extract_focal_v5.py \
    --input_dir results/07_macse/sperm/macse_aligned \
    --output_dir results/focal_alignments_17genes \
    --genes ADIG,ATP1A4,CABS1,DNHD1,DPY19L2,FAM71B,FAM71F2,FSIP2,MROH2B,NDC80,NPM2,RIMBP3,SPATA33,SPEM3,TMPRSS12,TTLL8,ZAN
"""

import argparse
import os
import sys
from collections import OrderedDict

# =============================================================================
# SPECIES SELECTION + PHYLOGENETIC ORDER
# =============================================================================
# The species to display (focal mole-rats + key outgroups + human), in the
# order they should appear top-to-bottom in the plot. Only species actually
# present in a given alignment are extracted; the rest are marked absent.
DISPLAY_ORDER = [
    # Bathyergidae (focal)
    "NMR",    # HetGla - naked mole-rat (eusocial)
    "BATS",   # BatSui - Cape dune mole-rat (solitary)
    "AMR",    # FukAns - Ansell's mole-rat (cooperative)
    # Hystricomorpha outgroups
    "GPIG",   # guinea pig
    "CHIN",   # chinchilla
    "DEGU",   # degu
    # Murid outgroups
    "MOUS",   # mouse
    "RAT",    # rat
    # Primate reference (well-conserved anchor)
    "HUMA",   # human
    # broader-mammal outgroup (signal phylogenetic breadth of the 23-sp analysis)
    "ELEP",   # elephant (Afrotheria - most distant major branch from rodents)
]

# Full display labels + taxonomic group (group drives label color in R)
SPECIES_INFO = OrderedDict([
    ("NMR",  {"label": "Naked mole-rat (HetGla)",   "group": "Bathyergidae"}),
    ("BATS", {"label": "Cape dune mole-rat (BatSui)","group": "Bathyergidae"}),
    ("AMR",  {"label": "Ansell's mole-rat (FukAns)", "group": "Bathyergidae"}),
    ("GPIG", {"label": "Guinea pig",                 "group": "Hystricomorpha"}),
    ("CHIN", {"label": "Chinchilla",                 "group": "Hystricomorpha"}),
    ("DEGU", {"label": "Degu",                        "group": "Hystricomorpha"}),
    ("MOUS", {"label": "Mouse",                       "group": "Muridae"}),
    ("RAT",  {"label": "Rat",                         "group": "Muridae"}),
    ("HUMA", {"label": "Human",                       "group": "Primates"}),
    ("BONO", {"label": "Bonobo",                      "group": "Primates"}),
    ("ELEP", {"label": "Elephant",                     "group": "Afrotheria"}),
])

# WHY A SPECIES IS ABSENT -- these are NOT the same claim and must not share a
# colour in the figure.
#
#   unannotated   the gene is not annotated in that species' proteome and its
#                 orthogroup contains no protein from it. tBLASTn recovers the
#                 locus from the genome for several of these (SPEM3, CABS1 in
#                 gpig/chin/degu at ~51-59% identity), so the gene very likely
#                 EXISTS and is simply unannotated. Resolving it properly is a
#                 whole-genome-alignment problem (TOGA), out of scope here.
#                 Label: "unannotated / not detected"  -- never "GENE ABSENT".
#
#   not_single_copy  NMR and BatSui sequences come from curated single-copy
#                 TSVs. A gene absent from that call (multi-copy, tandem
#                 duplicated, or uncalled) yields no sequence. This is a
#                 filtering decision, not evidence of gene loss.
#                 Label: "excluded (not single-copy)"
#
#   real_absence  no hit in the genome either. Only elephant SPEM3/ZAN met this
#                 in our tBLASTn screen. Label: "not detected in genome"
ABSENCE_REASON = {
    "NMR":  "not_single_copy",
    "BATS": "not_single_copy",
    "AMR":  "unannotated",       # BRAKER models, no curated single-copy call
    "GPIG": "unannotated", "CHIN": "unannotated", "DEGU": "unannotated",
    "MOUS": "unannotated", "RAT": "unannotated", "HUMA": "unannotated",
    "BONO": "unannotated", "ELEP": "unannotated",
}

TARGET_GENES_DEFAULT = [
    "ADIG", "ATP1A4", "CABS1", "DNHD1", "DPY19L2",
    "FAM71B", "FAM71F2", "FSIP2", "MROH2B", "NDC80",
    "NPM2", "RIMBP3", "SPATA33", "SPEM3", "TMPRSS12",
    "TTLL8", "ZAN",
]

def parse_fasta(path):
    """Return OrderedDict header->sequence."""
    seqs, cur, buf = OrderedDict(), None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            seqs[cur] = "".join(buf)
    return seqs

# Transcript-ID prefix -> species label. The 17-gene pipeline MACSE alignments
# use raw transcript IDs (mikado.chr..., ENST..., g123.t1), NOT clean labels.
import re
ID_PREFIX_MAP = [
    ("mikado.chr", "NMR"),    # HetGla naked mole-rat
    ("mikado.ptg", "BATS"),   # BatSui Cape dune mole-rat
    ("ENSCPOT",    "GPIG"),   # guinea pig
    ("ENSCLAT",    "CHIN"),   # chinchilla
    ("ENSODET",    "DEGU"),   # degu
    ("ENSMUST",    "MOUS"),   # mouse
    ("ENSRNOT",    "RAT"),    # rat
    ("ENST0",      "HUMA"),   # human (ENST0... ; before generic ENS)
    ("ENSBTAT",    "BOS"), ("ENSSSCT", "PIG"), ("ENSCAFT", "DOG"),
    ("ENSFCTT",    "CAT"), ("ENSLAFT", "ELEP"), ("ENSGGOT", "GORI"),
    ("ENSMMUT",    "MACA"), ("ENSOCUT", "RABB"), ("ENSSVLT", "SQUI"),
    ("ENSACHT",    "ACAH"), ("ENSMSPT", "WWMO"), ("ENSPEMT", "NADE"),
    ("ENSMOCT",    "PVOL"),
]

# Bonobo is NCBI RefSeq. In the 17-gene MACSE alignments its header is the full
# CDS locus tag, e.g.  lcl|NC_073270.2_cds_XP_054960461.2_74547
# The old fallback (tok.split("_")[0].upper()) turned that into "LCL|NC" and the
# sequence was silently discarded. Bonobo is now present in ADIG / RIMBP3 /
# DPY19L2 after the paralog-blocklist fix, so this must resolve.
BONOBO_RE = re.compile(r"(?:^|_)(?:cds_)?[NX]P_\d+")

# collected during a run so unmapped headers are reported, never dropped silently
UNMAPPED = []

def species_of(header):
    """Map a header to a species label.

    Handles three header conventions:
      atlas macse_aligned  bare labels           >NMR  >HUMA
      17-gene macse        transcript IDs        mikado.chr21G1260.1, ENST00000..., g11462.t1
      bonobo (RefSeq)      CDS locus tag         lcl|NC_073270.2_cds_XP_054960461.2_74547

    Returns None if the header cannot be resolved, so the caller can warn
    instead of silently inventing a species code.
    """
    tok = header.split()[0]

    # Package-A copy labels (HUMA_c2, ...) resolve to the base species; the
    # duplicate-warning path below then reports them instead of UNMAPPED noise
    m = re.match(r"^([A-Za-z]{2,5})_c\d+$", tok)
    if m:
        return m.group(1).upper()

    # already a clean species label (atlas style)
    if tok.upper() in SPECIES_INFO or (tok.isalpha() and 2 <= len(tok) <= 5):
        return tok.upper()

    if re.match(r"^g\d+\.t\d+", tok):
        return "AMR"                     # FukAns (BRAKER gene models)

    for prefix, sp in ID_PREFIX_MAP:
        if tok.startswith(prefix):
            return sp

    if BONOBO_RE.search(tok):
        return "BONO"                    # Pan paniscus, NCBI RefSeq

    UNMAPPED.append(tok)
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_dir",
                    default="results/07_macse/sperm/macse_aligned")
    ap.add_argument("--output_dir",
                    default="results/focal_alignments_17genes")
    ap.add_argument("--genes", default=",".join(TARGET_GENES_DEFAULT),
                    help="comma-separated gene list")
    ap.add_argument("--selected_from", default=None,
                    help="a part-9 scheme's pruned/ dir; when given, each "
                         "species' displayed sequence is the copy ortholog "
                         "selection kept (default: first copy, as before)")
    args = ap.parse_args()

    genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 70)
    print("01_extract_focal_v5: focal + outgroup extraction from MACSE alignments")
    print("=" * 70)
    print(f"  input : {args.input_dir}")
    print(f"  output: {args.output_dir}")
    print(f"  genes : {len(genes)}")
    print(f"  species shown: {', '.join(DISPLAY_ORDER)}")
    print("=" * 70)

    summary = []
    for gene in genes:
        _stem = {"FAM71B": "GARIN3", "FAM71F2": "GARIN1A"}.get(gene, gene)
        aln_path = os.path.join(args.input_dir, f"{_stem}_NT.aln")
        if not os.path.exists(aln_path):
            print(f"  {gene}: NO MACSE alignment ({aln_path}) -> skipped")
            summary.append((gene, "no_alignment", 0, 0))
            continue

        seqs = parse_fasta(aln_path)

        # which copy index per species did part 9 keep? (only with --selected_from)
        wanted_idx = {}
        if args.selected_from:
            pruned = os.path.join(args.selected_from, f"{_stem}_NT.fa")
            if os.path.exists(pruned):
                for uid in parse_fasta(pruned):
                    m = re.match(rf"^([A-Za-z]{{2,5}})_{_stem}(?:_c(\d+))?$", uid)
                    if m:
                        sp, k = m.group(1).upper(), int(m.group(2) or 1)
                        if sp in wanted_idx:
                            print(f"    NOTE {gene}: {sp} has multiple selected "
                                  f"copies (test-exempt) -- displaying the first")
                        else:
                            wanted_idx[sp] = k
            else:
                print(f"    WARNING {gene}: {pruned} not found -- "
                      f"falling back to first-copy display")

        by_species = {}
        seen_count = {}
        aln_len = 0
        for hdr, seq in seqs.items():
            sp = species_of(hdr)
            if sp is None:
                continue                 # reported at the end, not dropped silently
            seen_count[sp] = seen_count.get(sp, 0) + 1
            if sp in wanted_idx:
                if seen_count[sp] != wanted_idx[sp]:
                    continue             # not the selected copy
            elif sp in by_species:
                print(f"    WARNING {gene}: two sequences map to {sp} "
                      f"(second header: {hdr[:50]}) -- keeping the first")
                continue
            by_species[sp] = seq
            aln_len = max(aln_len, len(seq))

        # build the ordered sub-alignment + presence table
        out_records = []          # (species, sequence)  present species only
        presence_rows = []        # (species, present 0/1, reason)
        n_present = 0
        for sp in DISPLAY_ORDER:
            seq = by_species.get(sp)
            informative = seq is not None and set(seq.upper()) - set("-!?N") != set()
            if informative:
                out_records.append((sp, seq))
                presence_rows.append((sp, 1, ""))
                n_present += 1
            else:
                presence_rows.append((sp, 0, ABSENCE_REASON.get(sp, "unannotated")))

        if n_present < 2:
            print(f"  {gene}: only {n_present} of the chosen species present "
                  f"-> plot may be sparse (kept anyway)")

        # write sub-alignment (present species, in DISPLAY_ORDER)
        out_aln = os.path.join(args.output_dir, f"{gene}_focal_NT.aln")
        with open(out_aln, "w") as f:
            for sp, seq in out_records:
                f.write(f">{sp}\n")
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i + 60] + "\n")

        # write presence/absence table
        out_pres = os.path.join(args.output_dir, f"{gene}_presence.tsv")
        with open(out_pres, "w") as f:
            f.write("species\tpresent\treason\n")
            for sp, p, why in presence_rows:
                f.write(f"{sp}\t{p}\t{why}\n")

        n_absent = len(DISPLAY_ORDER) - n_present
        print(f"  {gene}: {n_present} present, {n_absent} absent, aln_len={aln_len}")
        summary.append((gene, "ok", n_present, n_absent))

    # write a master species-info table for the R script (labels + groups + order)
    info_path = os.path.join(args.output_dir, "species_info.tsv")
    with open(info_path, "w") as f:
        f.write("species\torder\tlabel\tgroup\n")
        for i, sp in enumerate(DISPLAY_ORDER):
            info = SPECIES_INFO[sp]
            f.write(f"{sp}\t{i}\t{info['label']}\t{info['group']}\n")

    print("=" * 70)
    print(f"  wrote species_info.tsv (labels/groups/order for the R plotter)")
    print(f"  {sum(1 for s in summary if s[1]=='ok')}/{len(genes)} genes extracted")
    if UNMAPPED:
        uniq = sorted(set(UNMAPPED))
        print(f"\n  *** {len(UNMAPPED)} headers could not be mapped to a species "
              f"({len(uniq)} distinct) ***")
        for u in uniq[:10]:
            print(f"      {u}")
        if len(uniq) > 10:
            print(f"      ... and {len(uniq)-10} more")
        print("  These were skipped. If any belong to a displayed species, the")
        print("  figure is WRONG -- add the prefix to ID_PREFIX_MAP.")
    else:
        print("  all headers mapped to a species")
    print("=" * 70)

if __name__ == "__main__":
    main()
