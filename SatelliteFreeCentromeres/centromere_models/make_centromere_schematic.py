#!/usr/bin/env python3
"""
Supplementary figure: the sequence neighbourhood of naked mole-rat satellite-free
centromeres (SFCs) compared with the satellite-based centromeres of human and mouse.

One idea per figure: CENP-A (red) is centred at the same position in every row,
so the eye compares only what surrounds it. Schematic, not to scale.

Outputs, into this folder:
  centromere_models.svg   - the figure (legend embedded in the SVG)
  centromere_models.png   - rendered FROM the SVG with rsvg-convert

Sources for the drawn layouts:
  Human  - CENP-A inside an alpha-satellite HOR array (171-bp monomers), flanked by
           pericentric satellites. Altemose et al. 2022; Gao et al. 2026.
  Mouse  - telocentric: telomere, minor satellite (120 bp) with CENP-A, then major
           satellite (234 bp) pericentromere. Liu et al. 2024; Chaudhry et al. 2025.
  NMR    - 22 of 31 centromeres (incl. chrX, chrY) have no centromeric satellite
           (Fig. 2B). Neighbourhoods per manuscript: olfactory-receptor clusters
           (chr2, 8, 16, 19, 20), Cyp450 (chr10, 28), tandem gene duplications
           (chr3, 9, 12, 13, 26, Y) and LINE arrays (chr1, 6, 15, 17, 18, 22, 24,
           25, X); checked against mikado gene models in cenpa_domain_composition.tsv.

Do not hand-edit the SVG; change this script and regenerate.
"""
import math
import os
import subprocess
from html import escape

C = dict(
    ink="#1A1D21", muted="#5C6169", faint="#8A8F97", hair="#DADDE2", bg="#FFFFFF",
    chrom="#AEB3BA", chrom_s="#4B5058",
    cenpa="#D1263B", cenpa_s="#8E1224",
    sat="#7B3FA8", sat2="#9C68C6", sat3="#5A2A85",
    peri="#3F4A5C",
    telo="#D4A017", telo_s="#8F6B00",
    gene="#1F6FC5", gene_s="#12457D",
    l1="#E07000", l1_s="#924800",
    sfc="#139460", sfc_tint="#E9F6EF",
    dna="#2B2F36",
)
SANS = "Helvetica, Arial, sans-serif"

W, H = 1400, 660
STRIP_X, STRIP_W = 376, 984          # CENP-A centred at the strip midpoint
BAR = 28                              # height of a repeat array
ROW_H = 150
ROW0 = 118


def esc(s):
    return escape(str(s), quote=True)


class Draw:
    def __init__(self):
        self.p = []

    def rect(self, x, y, w, h, fill, stroke=None, sw=1, rx=0):
        a = f'x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" fill="{fill}"'
        if rx:
            a += f' rx="{rx}"'
        if stroke:
            a += f' stroke="{stroke}" stroke-width="{sw}"'
        self.p.append(f"<rect {a}/>")

    def line(self, x1, y1, x2, y2, stroke, sw=1, cap=None):
        a = f'x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{stroke}" stroke-width="{sw}"'
        if cap:
            a += f' stroke-linecap="{cap}"'
        self.p.append(f"<line {a}/>")

    def path(self, d, fill="none", stroke=None, sw=1, op=None):
        a = f'd="{d}" fill="{fill}" stroke-linejoin="round" stroke-linecap="round"'
        if stroke:
            a += f' stroke="{stroke}" stroke-width="{sw}"'
        if op is not None:
            a += f' opacity="{op}"'
        self.p.append(f"<path {a}/>")

    def circle(self, cx, cy, r, fill, stroke=None, sw=1):
        a = f'cx="{cx:.2f}" cy="{cy:.2f}" r="{r}" fill="{fill}"'
        if stroke:
            a += f' stroke="{stroke}" stroke-width="{sw}"'
        self.p.append(f"<circle {a}/>")

    def ellipse(self, cx, cy, rx, ry, fill, stroke=None, sw=1):
        self.p.append(f'<ellipse cx="{cx:.2f}" cy="{cy:.2f}" rx="{rx}" ry="{ry}" fill="{fill}" '
                      f'stroke="{stroke}" stroke-width="{sw}"/>')

    def text(self, x, y, s, size, fill, anchor="start", weight="normal"):
        self.p.append(f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" fill="{fill}" '
                      f'style="font-family:{SANS};font-size:{size}px;font-weight:{weight};">{esc(s)}</text>')

    def out(self):
        return "".join(self.p)


def defs():
    h = BAR
    chev = [(0, 0), (9, 0), (13, h / 2), (9, h), (0, h), (4, h / 2)]

    def poly(dx, fill):
        pts = " ".join(f"{x + dx:g},{y:g}" for x, y in chev)
        return f'<polygon points="{pts}" fill="{fill}" stroke="{fill}" stroke-width="0.6"/>'

    return (
        "<defs>"
        f'<pattern id="hor" width="39" height="{h}" patternUnits="userSpaceOnUse">'
        f'<rect width="39" height="{h}" fill="#FFFFFF"/>'
        f'{poly(0, C["sat"])}{poly(13, C["sat2"])}{poly(26, C["sat3"])}</pattern>'
        f'<pattern id="sat" width="13" height="{h}" patternUnits="userSpaceOnUse">'
        f'<rect width="13" height="{h}" fill="#FFFFFF"/>{poly(0, C["sat"])}</pattern>'
        '<pattern id="peri" width="8" height="8" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'
        f'<rect width="8" height="8" fill="{C["peri"]}" opacity="0.35"/>'
        f'<rect width="3.5" height="8" fill="{C["peri"]}"/></pattern>'
        "</defs>"
    )


# ---------------------------------------------------------------- features (row-local y: bar top = 0)
def X(u):
    """strip units (0..1000) -> page x"""
    return STRIP_X + u * STRIP_W / 1000


def array(d, u0, u1, fill, stroke):
    d.rect(X(u0), 0, X(u1) - X(u0), BAR, fill, stroke, 1.2)


def dna(d, u0, u1):
    x0, x1, mid = X(u0), X(u1), BAR / 2
    for phase in (0, math.pi):
        pts, x = [], x0
        while x <= x1:
            pts.append(f"{x:.1f},{mid + 4 * math.sin((x - x0) / 16 * 2 * math.pi + phase):.1f}")
            x += 2
        d.path("M" + " L".join(pts), stroke=C["dna"], sw=1.5)


def gene(d, u0, u1):
    x0, x1, m = X(u0), X(u1), BAR / 2
    d.path(f"M{x0},{m - 11} L{x1 - 10},{m - 11} L{x1},{m} L{x1 - 10},{m + 11} L{x0},{m + 11} Z",
           fill=C["gene"], stroke=C["gene_s"], sw=1.3)


def l1(d, u0, u1):
    d.rect(X(u0), BAR / 2 - 10, X(u1) - X(u0), 20, C["l1"], C["l1_s"], 1.3, rx=4)


def telomere(d, u0, u1):
    d.rect(X(u0), 2, X(u1) - X(u0), BAR - 4, C["telo"], C["telo_s"], 1.2, rx=3)


def nucleosome(d, x, y, r=14):
    d.circle(x, y, r, C["cenpa"], C["cenpa_s"], 1.6)
    d.path(f"M{x - r + 2:.1f},{y + 3:.1f} Q{x:.1f},{y + 10:.1f} {x + r - 2:.1f},{y + 3:.1f}",
           stroke="#FFFFFF", sw=2, op=0.85)
    d.path(f"M{x - r + 4:.1f},{y - 5:.1f} Q{x:.1f},{y + 1:.1f} {x + r - 4:.1f},{y - 5:.1f}",
           stroke="#FFFFFF", sw=1.4, op=0.55)


def cenpa(d):
    """CENP-A domain, always centred on u=500."""
    for u in (420, 460, 500, 540, 580):
        nucleosome(d, X(u), BAR / 2)
    d.text(X(500), -20, "CENP-A", 15, C["cenpa"], anchor="middle", weight="bold")


def under(d, u, s, color, weight="bold", size=13.5):
    d.text(X(u), BAR + 26, s, size, color, anchor="middle", weight=weight)


# ---------------------------------------------------------------- chromosome icons
def icon_meta(d, cx, cy):
    for width, col in ((17, C["chrom_s"]), (13, C["chrom"])):
        for sx in (-1, 1):
            for sy in (-1, 1):
                d.line(cx + sx * 4, cy, cx + sx * 15, cy + sy * 40, col, width, cap="round")
    d.ellipse(cx, cy, 8.5, 6, C["cenpa"], C["cenpa_s"], 1.2)


def icon_telo(d, cx, cy):
    top = cy - 38
    for width, col in ((17, C["chrom_s"]), (13, C["chrom"])):
        for sx in (-1, 1):
            d.line(cx + sx * 5, top + 6, cx + sx * 11, cy + 42, col, width, cap="round")
    d.ellipse(cx, top + 5, 9, 6, C["cenpa"], C["cenpa_s"], 1.2)
    d.ellipse(cx, top - 4, 7, 2.8, C["telo"], C["telo_s"], 1)


# ---------------------------------------------------------------- rows
def row(d, i, name, sub, icon, draw, highlight=False):
    top = ROW0 + i * ROW_H
    bar_y = top + 62
    if highlight:
        d.rect(18, top, W - 36, ROW_H - 12, C["sfc_tint"], C["sfc"], 2.2, rx=12)
    d.text(40, bar_y + 6, name, 19, C["ink"], weight="bold")
    d.text(40, bar_y + 28, sub, 13.5, C["sfc"] if highlight else C["muted"],
           weight="bold" if highlight else "normal")
    icon(d, 326, bar_y + 14)
    g = Draw()
    draw(g)
    d.p.append(f'<g transform="translate(0,{bar_y})">{g.out()}</g>')


def human(g):
    array(g, 0, 150, "url(#peri)", C["peri"])
    array(g, 150, 850, "url(#sat)", C["sat"])
    array(g, 850, 1000, "url(#peri)", C["peri"])
    cenpa(g)
    under(g, 75, "pericentric satellite", C["peri"])
    under(g, 500, "α-satellite array", C["sat"])
    under(g, 925, "pericentric satellite", C["peri"])


def mouse(g):
    g.text(X(262), BAR / 2 + 5, "chromosome end", 12.5, C["faint"], anchor="end")
    telomere(g, 270, 300)
    array(g, 300, 680, "url(#sat)", C["sat"])
    array(g, 680, 1000, "url(#peri)", C["peri"])
    cenpa(g)
    under(g, 285, "telomere", C["telo_s"], weight="normal", size=12.5)
    under(g, 500, "minor satellite", C["sat"])
    under(g, 840, "major satellite", C["peri"])


def nmr(g):
    # Manuscript: SFCs sit in gene clusters (OR, Cyp450, tandem duplicates such as
    # Ptgr1 x12, Camk2d, Adad1, Cyp2b) or in LINE arrays; no satellite array.
    dna(g, 0, 1000)
    gene(g, 60, 160)
    gene(g, 180, 280)
    gene(g, 300, 400)
    gene(g, 440, 560)
    l1(g, 630, 700)
    l1(g, 720, 790)
    l1(g, 810, 880)
    cenpa(g)
    under(g, 280, "tandem gene copies", C["gene"])
    under(g, 500, "no satellite DNA", C["sfc"])
    under(g, 755, "L1 (LINE) array", C["l1"])


# ---------------------------------------------------------------- legend
def legend(d, y):
    items = [
        (lambda x: nucleosome(d, x + 12, y - 5, r=11), 26, "CENP-A (marks the centromere)"),
        (lambda x: d.rect(x, y - 16, 34, 20, "url(#sat)", C["sat"], 1), 34, "Centromeric satellite repeat"),
        (lambda x: d.rect(x, y - 16, 34, 20, "url(#peri)", C["peri"], 1), 34, "Pericentric satellite repeat"),
        (lambda x: d.path(f"M{x},{y - 14} L{x + 22},{y - 14} L{x + 30},{y - 6} L{x + 22},{y + 2} L{x},{y + 2} Z",
                          fill=C["gene"], stroke=C["gene_s"], sw=1.2), 30, "Gene"),
        (lambda x: d.rect(x, y - 13, 28, 14, C["l1"], C["l1_s"], 1.2, rx=3), 28, "L1 (LINE-1) sequence"),
        (lambda x: d.rect(x, y - 15, 16, 18, C["telo"], C["telo_s"], 1.2, rx=3), 16, "Telomere"),
    ]
    widths = [w + 10 + len(t) * 7.4 for _, w, t in items]
    gap = 34
    x = (W - (sum(widths) + gap * (len(items) - 1))) / 2
    for (draw, w, t), tw in zip(items, widths):
        draw(x)
        d.text(x + w + 10, y, t, 13.5, C["ink"])
        x += tw + gap


def build():
    d = Draw()
    d.rect(0, 0, W, H, C["bg"])
    d.text(40, 44, "Most naked mole-rat centromeres have no satellite DNA", 26, C["ink"], weight="bold")
    d.text(40, 72, "Sequence around the centromere (CENP-A, red). Human and mouse centromeres sit inside "
                   "long satellite repeats; most naked mole-rat centromeres sit in gene clusters or L1 (LINE) arrays.",
           15, C["muted"])
    row(d, 0, "Human", "all 24 chromosomes", icon_meta, human)
    row(d, 1, "Mouse", "all 21 chromosomes", icon_telo, mouse)
    row(d, 2, "Naked mole-rat", "22 of 31 chromosomes, incl. X and Y", icon_meta, nmr, highlight=True)
    legend(d, ROW0 + 3 * ROW_H + 26)
    d.text(W / 2, H - 16, "Schematic, not to scale.", 12, C["faint"], anchor="middle")
    return d.out()


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    svg = os.path.join(here, "centromere_models.svg")
    png = os.path.join(here, "centromere_models.png")
    with open(svg, "w", encoding="utf-8") as f:
        f.write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
                f'font-family="{SANS}">{defs()}{build()}</svg>')
    print("wrote", os.path.basename(svg))
    subprocess.run(["rsvg-convert", "-w", str(W * 2), "-h", str(H * 2), svg, "-o", png], check=True)
    print("wrote", os.path.basename(png))
