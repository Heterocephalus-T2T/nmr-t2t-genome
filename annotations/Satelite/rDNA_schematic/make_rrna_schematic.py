#!/usr/bin/env python3
"""
Generate the rRNA large-subunit schematic (Fig. 2D) as vector SVG.
Outputs, into this folder:
  rRNA_large_subunit.svg    - the large-subunit panel (the Fig. 2D schematic)
  rRNA_small_subunit.svg    - the 5S / TRC10 panel
  rRNA_schematic_combined.svg - both panels + legend (landscape)

All coordinates are drawn from the mHetGla1.1.pri annotation; the single 45S
repeat unit is to scale in base pairs. Numbers follow the manuscript:
TRC7 = 681 bp, TRC10 = 2,905 bp, chr29 rRNA satellite ~2 Mb.
"""
from html import escape

# ---- print palette (light) ----
C = dict(ink="#1A1D21", muted="#5C6169", faint="#8A8F97", hair="#E7E2D7",
         hair2="#D8D2C4", igs="#DED8CA", igs_s="#BEB6A4", s18="#2E8B8B",
         s58="#E0A32E", s28="#2C5C8F", s5="#3F9A5C", sat="#8E57A0",
         pol="#334049", inf="#A79E8C", accent="#1F6F6F", bg="#FFFFFF",
         k4me3="#C24030", k27ac="#E28A2A", k4me1="#C7A544")
SANS = "Helvetica, Arial, sans-serif"
MONO = "'Courier New', Courier, monospace"

# ---- one 45S repeat unit (genomic, - strand; drawn 5'->3' by flipping) ----
U_END, U_START = 16832922, 16806306
U_LEN = U_END - U_START            # 26616

def esc(s): return escape(str(s), quote=True)

class Draw:
    def __init__(self): self.p = []
    def rect(self, x, y, w, h, fill, stroke=None, sw=1, rx=0, dash=None, op=None):
        a = f'x="{x:.2f}" y="{y:.2f}" width="{max(w,0):.2f}" height="{h:.2f}"'
        if rx: a += f' rx="{rx}"'
        a += f' fill="{fill}"'
        if stroke: a += f' stroke="{stroke}" stroke-width="{sw}"'
        if dash: a += f' stroke-dasharray="{dash}"'
        if op is not None: a += f' opacity="{op}"'
        self.p.append(f'<rect {a}/>')
    def line(self, x1, y1, x2, y2, stroke, sw=1, dash=None, marker=None):
        a = f'x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{stroke}" stroke-width="{sw}"'
        if dash: a += f' stroke-dasharray="{dash}"'
        if marker: a += f' marker-end="url(#{marker})"'
        self.p.append(f'<line {a}/>')
    def circle(self, cx, cy, r, fill):
        self.p.append(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{r}" fill="{fill}"/>')
    def text(self, x, y, s, size, fill, anchor="start", weight="normal",
             italic=False, mono=False):
        st = (f'font-family:{MONO};' if mono else f'font-family:{SANS};')
        st += f'font-size:{size}px;font-weight:{weight};'
        if italic: st += 'font-style:italic;'
        self.p.append(f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" '
                      f'fill="{fill}" style="{st}">{esc(s)}</text>')
    def g(self, body, tx=0, ty=0):
        self.p.append(f'<g transform="translate({tx},{ty})">{body}</g>')

def defs():
    return (
        '<defs>'
        '<pattern id="sat" width="6" height="6" patternUnits="userSpaceOnUse" '
        f'patternTransform="rotate(45)"><rect width="6" height="6" fill="{C["sat"]}" '
        f'opacity="0.30"/><rect width="2.4" height="6" fill="{C["sat"]}"/></pattern>'
        f'<marker id="ph" markerWidth="9" markerHeight="9" refX="6" refY="3" '
        f'orient="auto"><path d="M0,0 L6,3 L0,6 Z" fill="{C["pol"]}"/></marker>'
        '</defs>')

# ---------------- LARGE SUBUNIT panel (900 x 300 internal) ----------------
def panel_large():
    d = Draw()
    # array strip (schematic)
    sx0, sx1, sy, sh = 20, 880, 34, 26
    d.rect(sx0, sy+7, sx1-sx0, sh-14, C["igs"], C["igs_s"], 1, 2)
    n = 6; uW = (sx1-sx0)/n
    for i in range(n):
        ux = sx0+i*uW; gx = ux+10; gW = uW*0.34
        d.rect(gx,            sy, gW*0.32, sh, C["s18"], rx=1.5)
        d.rect(gx+gW*0.32,    sy, gW*0.05, sh, C["s58"])
        d.rect(gx+gW*0.37,    sy, gW*0.63, sh, C["s28"], rx=1.5)
        d.rect(gx+gW+uW*0.20, sy+4, uW*0.20, sh-8, "url(#sat)", C["sat"], 0.8)
    d.text(sx0, sy-10, "tandem 45S units + TRC7 satellite (~2 Mb, chr29)", 11.5, C["ink"], weight="bold")
    d.text(sx1, sy-10, "detected by FISH →", 11, C["sat"], anchor="end", weight="bold")

    # to-scale unit
    x0, x1 = 60, 868; W = x1-x0; sc = W/U_LEN
    X   = lambda p: x0+(U_END-p)*sc
    seg = lambda a, b: (X(b), (b-a)*sc)
    mid = lambda a, b: X(b)+(b-a)*sc/2
    yG, hG, yTx = 132, 30, 104

    # connectors from strip unit #2 to the zoom
    d.line(sx0+1*uW+10, sy+sh, x0,  yTx-8, C["hair2"], 1, "3 3")
    d.line(sx0+2*uW,    sy+sh, x1,  yTx-8, C["hair2"], 1, "3 3")

    d.rect(x0, yG+6, W, hG-12, C["igs"], C["igs_s"], 1, 2)
    for a, b in [(16823617,16824737),(16820696,16823467)]:      # ITS measured
        x, w = seg(a, b); d.rect(x, yG+4, w, hG-8, C["igs"], C["igs_s"], 0.8)
    for a, b in [(16826603,16831923),(16815835,16816835)]:      # ETS inferred
        x, w = seg(a, b); d.rect(x, yG+4, w, hG-8, C["inf"], C["inf"], 1, dash="4 3", op=0.28)
    # TRC7 is one continuous tandem satellite. The annotated monomers are broken
    # only by degenerate copies RepeatMasker mislabels as simple-repeat (start-to-
    # start spacing is ~800 bp, with occasional ~1540 bp = one missed monomer), so
    # draw a single block rather than three islands.
    t7a, t7b = 16810636, 16813673
    x7, w7 = seg(t7a, t7b); d.rect(x7, yG+2, w7, hG-4, "url(#sat)", C["sat"], 0.9)
    d.text((X(t7b)+X(t7a))/2, yG+hG+14, "TRC7 satellite (~681 bp monomer, tandem)",
           9.5, C["sat"], anchor="middle", weight="bold")
    for a, b, c in [(16824737,16826603,"s18"),(16823467,16823617,"s58"),(16816835,16820696,"s28")]:
        x, w = seg(a, b); d.rect(x, yG, max(w,2), hG, C[c], "rgba(0,0,0,.15)", 0.6, 2.5)
    for n2, a, b, c in [("18S",16824737,16826603,"s18"),("28S",16816835,16820696,"s28")]:
        d.text(mid(a,b), yG-9, n2, 13, C[c], anchor="middle", weight="bold")
    m58 = mid(16823467,16823617)
    d.line(m58, yG-5, m58, yG-1, C["s58"], 1)
    d.text(m58, yG-9, "5.8S", 10.5, C["s58"], anchor="middle", weight="bold")
    for n2, a, b, c, it in [("5′ETS",16826603,16831923,"inf",True),
                            ("ITS1",16823617,16824737,"muted",False),
                            ("ITS2",16820696,16823467,"muted",False),
                            ("3′ETS",16815835,16816835,"inf",True)]:
        if (b-a)*sc < 24: continue
        d.text(mid(a,b), yG+hG+14, n2, 9.5, C[c], anchor="middle", italic=it)

    tS, tE = X(16826603), X(16816835)
    d.line(tS, yTx, tE-4, yTx, C["pol"], 2, marker="ph")
    d.circle(tS, yTx, 3, C["pol"])
    d.text(tS+6, yTx-7, "Pol I", 10, C["pol"], weight="bold")
    d.text((tS+tE)/2, yTx-7, "45S pre-rRNA", 10, C["muted"], anchor="middle", italic=True)

    # ---- chromatin marks overlaid on the unit (liver ChIP-seq; K4me3 also heart) ----
    # Manuscript: peaks overlaid on gene/repeat tracks, asking which feature each falls on.
    lanes = [("H3K4me3", "k4me3", [(16831707,16832151),(16809632,16810217)]),   # promoter
             ("H3K27ac", "k27ac", [(16809902,16810085)]),                        # active enhancer
             ("H3K4me1", "k4me1", [(16831137,16832026),(16809328,16810106)])]    # enhancer
    yc = 182
    for name, col, peaks in lanes:
        d.line(x0, yc+3.5, x1, yc+3.5, C["hair"], 1)
        d.text(x0-5, yc+7, name, 8.5, C[col], anchor="end", weight="bold")
        for a, b in peaks:
            xx, ww = seg(a, b); d.rect(xx, yc, max(ww,3), 7, C[col], rx=1.5)
        yc += 10
    # functional zones: promoter (5′) vs distal-IGS enhancer (3′)
    def zone(a, b, lab, sub, col):
        L, R = X(b), X(a); yb = 216
        d.line(L, yb, R, yb, C[col], 1.2); d.line(L, yb-4, L, yb, C[col], 1.2); d.line(R, yb-4, R, yb, C[col], 1.2)
        d.text((L+R)/2, yb+12, lab, 10, C[col], anchor="middle", weight="bold")
        d.text((L+R)/2, yb+22, sub, 8.5, C["muted"], anchor="middle")
    zone(16831137, 16832151, "Pol I promoter", "K4me3 · K4me1", "k4me3")
    zone(16809328, 16810217, "distal-IGS enhancer", "K27ac · K4me1 · K4me3", "k27ac")

    yR = 250
    d.line(x0, yR, x1, yR, C["hair2"], 1)
    dd = 0
    while dd <= U_LEN:
        x = x0+dd*sc; d.line(x, yR, x, yR+5, C["faint"], 1)
        d.text(x, yR+15, f'{dd//1000}{" kb" if dd==0 else ""}', 9.5, C["faint"], anchor="middle", mono=True)
        dd += 5000
    d.text((x0+x1)/2, yR+31, "one 45S repeat unit ≈ 26 kb, to scale (5′→3′)",
           10.5, C["muted"], anchor="middle", weight="bold")
    d.text((x0+x1)/2, yR+44, "ChIP-seq: H3K4me3 liver + heart; H3K27ac & H3K4me1 liver (heart H3K4me1 excluded — noisy, 152k peaks genome-wide)",
           8.5, C["faint"], anchor="middle")
    d.text((x0+x1)/2, yR+56, "CENPA absent from the array (0 peaks; chr29 centromere ~3.0–3.3 Mb, ~14 Mb away)",
           8.5, C["faint"], anchor="middle")
    return "".join(d.p)

# ---------------- SMALL / 5S panel (900 x 120 internal) ----------------
def panel_small():
    d = Draw()
    ax0, ax1, ay, ah = 20, 430, 40, 26
    d.rect(ax0, ay, ax1-ax0, ah, "url(#sat)", C["sat"], 1, 2)
    for i in range(13):
        x = ax0+8+i*((ax1-ax0-16)/12); d.rect(x-1.5, ay, 3, ah, C["s5"])
    d.text(ax0, ay-10, "single 5S satellite (~196 kb, chr6)", 11.5, C["ink"], weight="bold")
    d.text((ax0+ax1)/2, ay+ah+18, "~69 tandem TRC10 monomers", 10, C["muted"], anchor="middle")
    d.line(460, 24, 460, 96, C["hair"], 1)
    bx0, bx1, by, bh, mlen = 500, 880, 40, 26, 2905; ms = (bx1-bx0)/mlen
    d.rect(bx0, by, bx1-bx0, bh, "url(#sat)", C["sat"], 1, 2)
    g5x = bx0+(2905-114-170)*ms; g5w = max(114*ms, 4)
    d.rect(g5x, by, g5w, bh, C["s5"], "rgba(0,0,0,.15)", 0.6, 1.5)
    d.text(g5x+g5w/2, by-8, "5S", 11.5, C["s5"], anchor="middle", weight="bold")
    d.text((bx0+bx1)/2, by-24, "one TRC10 monomer — to scale", 11, C["ink"], anchor="middle", weight="bold")
    yR = by+bh+8; d.line(bx0, yR, bx1, yR, C["hair2"], 1)
    dd = 0
    while dd <= 2905:
        x = bx0+dd*ms; d.line(x, yR, x, yR+5, C["faint"], 1)
        d.text(x, yR+16, f'{dd//1000}{" kb" if dd==0 else ""}', 9.5, C["faint"], anchor="middle", mono=True)
        dd += 1000
    d.text((bx0+bx1)/2, yR+30, "114-bp 5S gene · remainder is TRC10 satellite (2,905 bp)",
           9.5, C["muted"], anchor="middle")
    return "".join(d.p)

def legend_row(y):
    d = Draw()
    items = [(C["s18"],"18S",False),(C["s58"],"5.8S",False),(C["s28"],"28S",False),
             (C["s5"],"5S",False),(C["igs"],"IGS / ITS1 / ITS2",False),
             (C["inf"],"5′/3′ETS (inferred)",False),("url(#sat)","TRC7 / TRC10 satellite",False),
             (C["pol"],"Pol I pre-rRNA",True),
             (C["k4me3"],"H3K4me3 (promoter)",False),(C["k27ac"],"H3K27ac (active enh.)",False),
             (C["k4me1"],"H3K4me1 (enhancer)",False)]
    x = 20
    for fill, lab, isline in items:
        w = 34 + len(lab)*7.0
        if x + w > 890:                       # wrap to next line
            x = 20; y += 20
        if isline: d.line(x, y-3, x+20, y-3, fill, 3)
        else: d.rect(x, y-11, 20, 12, fill, "rgba(0,0,0,.12)", 1, 3)
        d.text(x+26, y, lab, 12, C["muted"])
        x += w
    return "".join(d.p)

def svg_doc(width, height, body):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
            f'width="{width}" height="{height}" font-family="{SANS}">'
            f'{defs()}<rect width="{width}" height="{height}" fill="{C["bg"]}"/>'
            f'{body}</svg>')

def group(body, tx=0, ty=0):
    return f'<g transform="translate({tx},{ty})">{body}</g>'

def title(t, sub, w):
    d = Draw()
    d.text(20, 22, t, 16, C["ink"], weight="bold")
    d.text(20, 40, sub, 12, C["muted"])
    return "".join(d.p)

if __name__ == "__main__":
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    def write(name, s):
        with open(os.path.join(here, name), "w") as f: f.write(s)
        print("wrote", name)

    T = ("Full annotation of the rRNA large subunit",
         "Heterocephalus glaber (mHetGla1.1.pri) — Fig 2C companion to the rRNA (TRC7) FISH")

    # large only
    body = title(*T, 900) + group(panel_large(), 0, 50) + group(legend_row(0), 0, 372)
    write("rRNA_large_subunit.svg", svg_doc(900, 424, body))

    # small only
    body = group(panel_small(), 0, 30)
    write("rRNA_small_subunit.svg", svg_doc(900, 160, body))

    # combined
    body = (title(*T, 900)
            + group(panel_large(), 0, 50)
            + group(panel_small(), 0, 372)
            + group(legend_row(0), 0, 486))
    write("rRNA_schematic_combined.svg", svg_doc(900, 540, body))
