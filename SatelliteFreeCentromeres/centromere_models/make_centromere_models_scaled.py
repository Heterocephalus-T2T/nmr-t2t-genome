#!/usr/bin/env python3
"""
Third version: the centromere_models schematic (same glyphs and layout), with
arrays, gene copies, L1 blocks and CENP-A domains drawn to scale in a 5-Mb window
centred on CENP-A. Satellite texture and nucleosome glyphs are not to scale.
Arms that lie beyond the window are drawn as short stubs after a // break.

Sources:
  Human  - active alpha-satellite array 2.5 Mb and CENP-A/CDR 240 kb: Gao et al.
           2026 Nature (means). Satellite region either side of the active array:
           CHM13 censat v2.1 medians over 19 non-acrocentric chromosomes, p side
           2.24 Mb, q side 1.78 Mb (data/chm13_centromere_region_sizes.tsv).
  Mouse  - telomere -> minor satellite (core, CENP-A) -> major satellite; minor
           satellite 0.48-2.3 Mb, major satellite 2.2-23 Mb: Chaudhry et al. 2025
           Nat Commun. ~1 Mb and the ~0.8-Mb CENP-A extent are chosen within those
           ranges (approximate).
  NMR    - Ptgr1 cluster: the 11 real copies on chr13 (data/chr13_ptgr1_copies.bed),
           positions preserved. CENP-A domain 220 kb = median of the 22 SFC domains.
           L1 array: the real L1 blocks of the chrX CENP-A domain (L1 merged within
           1 kb, blocks >= 3 kb; data/chrX_cenpa_domain.l1_blocks.bed). The spacing
           between the two features is a layout choice.

Outputs: centromere_models_scaled.svg / .png (PNG rendered from the SVG).
Do not hand-edit the SVG; change this script and regenerate.
"""
import os
import subprocess

import make_centromere_schematic as M
from make_centromere_schematic import C, SANS, Draw, X, array, dna, telomere, nucleosome, under

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

# leave room either side of the 5-Mb strip for arm stubs
M.STRIP_X, M.STRIP_W = 430, 860
M.ROW0 = 18                 # rows start at the top now that the title is gone
KB = 1000 / 5000            # strip units per kb (5 Mb window = 1000 units)
W, H = M.W, 568


def u(kb):
    return kb * KB


def read_bed(name):
    with open(os.path.join(DATA, name)) as f:
        return [(int(p[1]), int(p[2])) for p in (line.split() for line in f)]


KB_PER_GLYPH = 120          # one CENP-A glyph per 120 kb of domain, same in every row
GLYPH_R = 10


def cenpa_domain(d, kb0, kb1, label):
    """CENP-A glyphs across a to-scale domain; glyph count is proportional to domain length."""
    a, b = u(kb0), u(kb1)
    n = max(1, round((kb1 - kb0) / KB_PER_GLYPH))
    step = (b - a) / n
    for i in range(n):                      # glyphs sit on top of the track so the DNA stays visible
        nucleosome(d, X(a + step * (i + 0.5)), -GLYPH_R + 1, r=GLYPH_R)
    d.text(X((a + b) / 2), -2 * GLYPH_R - 10, label, 14, C["cenpa"], anchor="middle", weight="bold")


def gene_copy(d, kb0, kb1):
    x0, x1, m = X(u(kb0)), X(u(kb1)), M.BAR / 2
    head = min(8, (x1 - x0) * 0.45)
    d.path(f"M{x0:.2f},{m - 10} L{x1 - head:.2f},{m - 10} L{x1:.2f},{m} L{x1 - head:.2f},{m + 10} "
           f"L{x0:.2f},{m + 10} Z", fill=C["gene"], stroke=C["gene_s"], sw=1)


def l1_block(d, kb0, kb1):
    d.rect(X(u(kb0)), M.BAR / 2 - 10, X(u(kb1)) - X(u(kb0)), 20, C["l1"], C["l1_s"], 0.8)


def arm_stub(d, side, note):
    """// break at the window edge, then a short stub of arm DNA (not to scale)."""
    m = M.BAR / 2
    if side == "left":
        bx, s0, s1 = X(0) - 9, X(0) - 62, X(0) - 18
    else:
        bx, s0, s1 = X(1000) + 9, X(1000) + 18, X(1000) + 62
    for dx in (-3, 3):
        d.line(bx + dx - 4, m + 12, bx + dx + 4, m - 12, C["ink"], 1.6)
    k = M.STRIP_W / 1000
    dna(d, (s0 - M.STRIP_X) / k, (s1 - M.STRIP_X) / k)
    d.text((s0 + s1) / 2, -6, "arm", 12, C["muted"], anchor="middle")
    d.text((s0 + s1) / 2, M.BAR + 26, note, 11.5, C["muted"], anchor="middle")


def human(g):
    array(g, 0, u(1250), "url(#peri)", C["peri"])
    array(g, u(1250), u(3750), "url(#sat)", C["sat"])
    array(g, u(3750), u(5000), "url(#peri)", C["peri"])
    cenpa_domain(g, 2380, 2620, "CENP-A · 240 kb")
    under(g, u(625), "pericentric satellite", C["peri"])
    under(g, u(2500), "α-satellite array · 2.5 Mb", C["sat"])
    under(g, u(4375), "pericentric satellite", C["peri"])
    arm_stub(g, "left", "+1.0 Mb")        # 2.24 Mb p side - 1.25 Mb shown
    arm_stub(g, "right", "+0.5 Mb")       # 1.78 Mb q side - 1.25 Mb shown


def mouse(g):
    g.text(X(u(1945)), M.BAR / 2 + 5, "chromosome end", 12.5, C["faint"], anchor="end")
    telomere(g, u(1975), u(2000))
    array(g, u(2000), u(3000), "url(#sat)", C["sat"])
    array(g, u(3000), u(5000), "url(#peri)", C["peri"])
    cenpa_domain(g, 2100, 2900, "CENP-A · ~0.8 Mb (approx.)")
    under(g, u(2500), "minor satellite · ~1 Mb", C["sat"])
    under(g, u(4000), "major satellite · 2–23 Mb in total", C["peri"])
    arm_stub(g, "right", "+0.2–21 Mb")


def nmr(g):
    dna(g, 0, 1000)
    copies = read_bed("chr13_ptgr1_copies.bed")
    c0, c1 = copies[0][0], copies[-1][1]
    shift = 2500 * 1000 - (c0 + c1) // 2              # centre the real cluster on CENP-A
    for a, b in copies:
        gene_copy(g, (a + shift) / 1000, (b + shift) / 1000)
    blocks = read_bed("chrX_cenpa_domain.l1_blocks.bed")
    b0 = blocks[0][0]
    for a, b in blocks:                                # real chrX L1 blocks, placed at 3.3 Mb
        l1_block(g, 3300 + (a - b0) / 1000, 3300 + (b - b0) / 1000)
    cenpa_domain(g, 2390, 2610, "CENP-A · 220 kb")
    under(g, u(2500), "Gene", C["gene"])
    under(g, u(3300 + (blocks[-1][1] - b0) / 2000), "L1 (LINE) array", C["l1"])
    under(g, u(1100), "no satellite DNA", C["sfc"])
    arm_stub(g, "left", "")
    arm_stub(g, "right", "")


def scale_bar(d, y):
    x0, x1 = X(0), X(u(1000))
    d.line(x0, y, x1, y, C["ink"], 2)
    for x in (x0, x1):
        d.line(x, y - 5, x, y + 5, C["ink"], 2)
    d.text(x1 + 10, y + 5, "1 Mb", 13, C["ink"], weight="bold")


def build():
    # no title, subtitle or footnotes: those live in the figure legend (caption.md)
    d = Draw()
    d.rect(0, 0, W, H, C["bg"])
    M.row(d, 0, "Human", "all 24 chromosomes", M.icon_meta, human)
    M.row(d, 1, "Mouse", "20 of 21 chromosomes", M.icon_telo, mouse)
    M.row(d, 2, "Naked mole-rat", "23 of 31 chromosomes, incl. X and Y", M.icon_meta, nmr, highlight=True)
    scale_bar(d, M.ROW0 + 3 * M.ROW_H + 16)
    M.legend(d, M.ROW0 + 3 * M.ROW_H + 66)
    return d.out()


if __name__ == "__main__":
    svg = os.path.join(HERE, "centromere_models_scaled.svg")
    png = os.path.join(HERE, "centromere_models_scaled.png")
    with open(svg, "w", encoding="utf-8") as f:
        f.write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
                f'font-family="{SANS}">{M.defs()}{build()}</svg>')
    print("wrote", os.path.basename(svg))
    subprocess.run(["rsvg-convert", "-w", str(W * 2), "-h", str(H * 2), svg, "-o", png], check=True)
    print("wrote", os.path.basename(png))
