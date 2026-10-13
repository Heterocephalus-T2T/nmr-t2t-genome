#!/usr/bin/env python3
"""
01_extract_controls_v1.py

Prepare control-gene MACSE alignments for protein-disruption plotting, the same
way 01_extract_focal_v5.py does for the sperm genes -- but control genes have
GENE-SPECIFIC test lineages, so each gene highlights the correct species.

Control design (from control_gene_validation_v21.py):
  POSITIVE controls, test = subterranean bathyergids (expect damage/relaxation):
    CRYAA CRYAB CRYBB1 RHO GRK1   -> NMR, BATS, AMR
    OPN1SW                        -> NMR, BATS
  POSITIVE controls, test = catarrhine primates (expect damage/relaxation):
    OR1A1 OR2W1                   -> HUMA, MACA, BONO, GORI
  NEGATIVE controls, test = primates (expect INTACT, conserved everywhere):
    ATP5F1A RPL13 RPS13 OPN1LW VCP GPI REEP5 VPS29 SNRPD3 -> HUMA, MACA, BONO, GORI

The disruption plot then shows:
  - positive controls: frameshifts/stops in the test lineage = visual confirmation
    of relaxed selection
  - negative controls: a clean conserved alignment = visual confirmation the gene
    is intact (no false damage)

Reads MACSE alignments (retains '!' frameshifts and internal stops that ClipKIT
strips). Writes, per gene:
  {out}/{GENE}_focal_NT.aln    sub-alignment (test lineage first, then references)
  {out}/{GENE}_presence.tsv    species | present(0/1) | reason
  {out}/species_info.tsv       label/group/order for the R plotter (per gene set)

Because the test lineage differs per gene, species_info is written PER GENE
(suffix by gene) so the R script orders each correctly.

Usage:
  python3 01_extract_controls_v1.py \
      --input_dir results/07_macse/validation/macse_aligned
      --output_dir results/control_gene_validation_v22/focal_alignments
"""

import argparse, os, re
from collections import OrderedDict

# ---- per-gene test lineage + type (parsed from control_gene_validation_v21) ----
CONTROL_GENES = OrderedDict([
    ("CRYAA",  {"type": "positive", "test": ["NMR", "BATS", "AMR"]}),
    ("CRYAB",  {"type": "positive", "test": ["NMR", "BATS", "AMR"]}),
    ("CRYBB1", {"type": "positive", "test": ["NMR", "BATS", "AMR"]}),
    ("RHO",    {"type": "positive", "test": ["NMR", "BATS", "AMR"]}),
    ("GRK1",   {"type": "positive", "test": ["NMR", "BATS", "AMR"]}),
    ("OPN1SW", {"type": "positive", "test": ["NMR", "BATS"]}),
    ("OR1A1",  {"type": "positive", "test": ["HUMA", "MACA", "BONO", "GORI"]}),
    ("OR2W1",  {"type": "positive", "test": ["HUMA", "MACA", "BONO", "GORI"]}),
    ("ATP5F1A",{"type": "negative", "test": ["HUMA", "MACA", "BONO", "GORI"]}),
    ("RPL13",  {"type": "negative", "test": ["HUMA", "MACA", "BONO", "GORI"]}),
    ("RPS13",  {"type": "negative", "test": ["HUMA", "MACA", "BONO", "GORI"]}),
    ("OPN1LW", {"type": "negative", "test": ["HUMA", "MACA", "BONO", "GORI"]}),
    ("VCP",    {"type": "negative", "test": ["HUMA", "MACA", "BONO", "GORI"]}),
    ("GPI",    {"type": "negative", "test": ["HUMA", "MACA", "BONO", "GORI"]}),
    ("REEP5",  {"type": "negative", "test": ["HUMA", "MACA", "BONO", "GORI"]}),
    ("VPS29",  {"type": "negative", "test": ["HUMA", "MACA", "BONO", "GORI"]}),
    ("SNRPD3", {"type": "negative", "test": ["HUMA", "MACA", "BONO", "GORI"]}),
])

# Reference species shown beneath the test lineage (broad, conserved outgroups).
# Kept consistent so every panel has the same visual reference frame.
REFERENCE_ORDER = ["HUMA", "BONO", "GORI", "MACA", "MOUS", "RAT",
                   "GPIG", "CHIN", "DEGU", "DOG", "CAT", "PIG", "BOS", "ELEP"]

SPECIES_LABEL = {
    "NMR": "Naked mole-rat (HetGla)", "BATS": "Cape dune mole-rat (BatSui)",
    "AMR": "Ansell's mole-rat (FukAns)", "HUMA": "Human", "BONO": "Bonobo",
    "GORI": "Gorilla", "MACA": "Macaque", "MOUS": "Mouse", "RAT": "Rat",
    "GPIG": "Guinea pig", "CHIN": "Chinchilla", "DEGU": "Degu", "DOG": "Dog",
    "CAT": "Cat", "PIG": "Pig", "BOS": "Cattle", "ELEP": "Elephant",
}
SPECIES_GROUP = {
    "NMR": "Bathyergidae", "BATS": "Bathyergidae", "AMR": "Bathyergidae",
    "HUMA": "Primates", "BONO": "Primates", "GORI": "Primates", "MACA": "Primates",
    "MOUS": "Muridae", "RAT": "Muridae", "GPIG": "Hystricomorpha",
    "CHIN": "Hystricomorpha", "DEGU": "Hystricomorpha", "DOG": "Laurasiatheria",
    "CAT": "Laurasiatheria", "PIG": "Laurasiatheria", "BOS": "Laurasiatheria",
    "ELEP": "Afrotheria",
}

BONOBO_RE = re.compile(r"(?:^|_)(?:cds_)?[NX]P_\d+")
ID_PREFIX_MAP = [
    ("mikado.chr", "NMR"), ("mikado.ptg", "BATS"),
    ("ENSCPOT", "GPIG"), ("ENSCLAT", "CHIN"), ("ENSODET", "DEGU"),
    ("ENSMUST", "MOUS"), ("ENSRNOT", "RAT"), ("ENST0", "HUMA"),
    ("ENSBTAT", "BOS"), ("ENSSSCT", "PIG"), ("ENSCAFT", "DOG"),
    ("ENSFCTT", "CAT"), ("ENSLAFT", "ELEP"), ("ENSGGOT", "GORI"),
    ("ENSMMUT", "MACA"),
]
UNMAPPED = []


def species_of(header):
    tok = header.split()[0]
    up = tok.upper()
    # bare species code
    if up in SPECIES_LABEL or (tok.isalpha() and 2 <= len(tok) <= 5):
        return up
    # SPECIES_GENE form (e.g. NMR_CRYAA, HUMA_VCP): first field is the species
    first = up.split("_")[0]
    if first in SPECIES_LABEL:
        return first
    if re.match(r"^g\d+\.t\d+", tok):
        return "AMR"
    for pre, sp in ID_PREFIX_MAP:
        if tok.startswith(pre):
            return sp
    if BONOBO_RE.search(tok):
        return "BONO"
    UNMAPPED.append(tok)
    return None


def parse_fasta(path):
    seqs, cur, buf = OrderedDict(), None, []
    for line in open(path):
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


def find_aln(input_dir, gene):
    for pat in (f"{gene}_NT.aln", f"{gene}_NT.fa", f"{gene}.aln",
                f"{gene}_codon.aln"):
        p = os.path.join(input_dir, pat)
        if os.path.exists(p):
            return p
    # OG-suffixed names
    import glob
    hits = glob.glob(os.path.join(input_dir, f"{gene}_*NT*.*")) + \
           glob.glob(os.path.join(input_dir, f"{gene}_OG*"))
    hits = [h for h in hits if h.endswith((".aln", ".fa", ".fasta"))]
    return hits[0] if hits else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_dir", required=True,
                    help="dir with control MACSE alignments (retaining '!')")
    ap.add_argument("--output_dir",
                    default="results/07_macse/validation/protein_disruption_figures")
    ap.add_argument("--selected_from", default=None,
                    help="part-9 pruned/ dir (e.g. results/09_relax/validation/"
                         "general/pruned): display the selected copy per species")
    args = ap.parse_args()
    # safety: never write control figures into the sperm-gene focal dir
    if "focal_alignments_17genes" in args.output_dir or \
       args.output_dir.rstrip("/").endswith("10_extended_sperm_genes/focal_alignments"):
        raise SystemExit(f"ERROR: {args.output_dir} is a SPERM-gene output dir. "
                         "Choose a control-specific --output_dir.")
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 70)
    print("01_extract_controls_v1: control-gene alignment extraction")
    print(f"  input : {args.input_dir}")
    print(f"  output: {args.output_dir}")
    print("=" * 70)

    info_rows = []  # (gene, species, order, label, group, is_test)
    done = 0
    for gene, meta in CONTROL_GENES.items():
        aln = find_aln(args.input_dir, gene)
        if not aln:
            print(f"  {gene}: no alignment found -> skipped")
            continue
        seqs = parse_fasta(aln)

        wanted_idx = {}
        if args.selected_from:
            pruned = os.path.join(args.selected_from, f"{gene}_NT.fa")
            if os.path.exists(pruned):
                for uid in parse_fasta(pruned):
                    m = re.match(rf"^([A-Za-z]{{2,5}})_{gene}(?:_c(\d+))?$", uid)
                    if m:
                        sp, k = m.group(1).upper(), int(m.group(2) or 1)
                        wanted_idx.setdefault(sp, k)
            else:
                print(f"  {gene}: {pruned} not found -- first-copy display")

        by_sp, seen = {}, {}
        for hdr, s in seqs.items():
            sp = species_of(hdr)
            if sp is None:
                continue
            seen[sp] = seen.get(sp, 0) + 1
            if sp in wanted_idx:
                if seen[sp] != wanted_idx[sp]:
                    continue
            elif sp in by_sp:
                print(f"    WARNING {gene}: two sequences map to {sp} "
                      f"(second: {hdr[:50]}) -- keeping the first")
                continue
            by_sp[sp] = s
        if not by_sp:
            print(f"  {gene}: no sequences mapped to a species -> skipped")
            continue

        # display order: test lineage first (in defined order), then references
        test = [s for s in meta["test"] if s in by_sp]
        refs = [s for s in REFERENCE_ORDER if s in by_sp and s not in meta["test"]]
        order = test + refs

        # write sub-alignment
        out_aln = os.path.join(args.output_dir, f"{gene}_focal_NT.aln")
        with open(out_aln, "w") as f:
            for sp in order:
                f.write(f">{sp}\n")
                seq = by_sp[sp]
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i + 60] + "\n")

        # presence table (test + reference species we chose to display)
        display = meta["test"] + [s for s in REFERENCE_ORDER if s not in meta["test"]]
        out_pres = os.path.join(args.output_dir, f"{gene}_presence.tsv")
        with open(out_pres, "w") as f:
            f.write("species\tpresent\treason\n")
            for sp in display:
                if sp in by_sp and set(by_sp[sp].upper()) - set("-!?N"):
                    f.write(f"{sp}\t1\t\n")
                else:
                    # absent: test-species bathyergids from single-copy filter -> "absent";
                    # everything else -> unannotated
                    reason = "not_single_copy" if sp in ("NMR", "BATS") else "unannotated"
                    f.write(f"{sp}\t0\t{reason}\n")

        # per-gene species_info (order + test flag drive the R plot)
        for i, sp in enumerate(display):
            info_rows.append((gene, sp, i, SPECIES_LABEL.get(sp, sp),
                              SPECIES_GROUP.get(sp, "Other"),
                              1 if sp in meta["test"] else 0))

        n_test = len(test)
        n_ref = len(refs)
        tag = meta["type"]
        print(f"  {gene:8s} [{tag}] test={n_test}/{len(meta['test'])} present, "
              f"refs={n_ref}, len={len(next(iter(by_sp.values())))}")
        done += 1

    # write combined species_info (gene-keyed so R can subset per panel)
    info_path = os.path.join(args.output_dir, "species_info.tsv")
    with open(info_path, "w") as f:
        f.write("gene\tspecies\torder\tlabel\tgroup\tis_test\n")
        for row in info_rows:
            f.write("\t".join(str(x) for x in row) + "\n")

    print("=" * 70)
    print(f"  {done}/{len(CONTROL_GENES)} genes extracted -> {args.output_dir}")
    if UNMAPPED:
        uniq = sorted(set(UNMAPPED))
        print(f"  WARNING: {len(UNMAPPED)} unmapped headers ({len(uniq)} distinct):")
        for u in uniq[:8]:
            print(f"    {u}")
        print("  If any belong to a displayed species, add its prefix to ID_PREFIX_MAP.")
    else:
        print("  all headers mapped to a species")
    print("=" * 70)


if __name__ == "__main__":
    main()
