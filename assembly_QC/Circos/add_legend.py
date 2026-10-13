#!/usr/bin/env python3
"""Add a track legend to the circos_v2 figure.

The legend lives IN the SVG, and the PNG is rendered FROM that SVG
(never composited onto the raster). This script is idempotent - it strips any legend
it added before, so it can be re-run after every `circos` render.

Pipeline:
    ~/software/circos-0.69-9/bin/circos -conf conf/circos_v2.conf -outputdir . -outputfile circos_v2.png
    python3 add_legend.py            # edits circos_v2.svg in place, re-renders circos_v2.png

Swatch colors are the *actual* RGB circos emitted (verified against the SVG), not the
named-color inputs, so the legend matches the vector exactly. If you change a track colour
in circos_v2.conf, update the matching swatch RGB here too.
"""
import os, re, subprocess

RUN = os.path.expanduser('~/Desktop/Circos/circos_run')
SVG = f'{RUN}/circos_v2.svg'
PNG = f'{RUN}/circos_v2.png'
RSVG = '/opt/homebrew/bin/rsvg-convert'

PLOT_W = 3000          # circos canvas (square)
PANEL  = 900           # added right-hand strip for the legend
NEW_W  = PLOT_W + PANEL
X0     = PLOT_W + 55   # left margin of legend content
SW     = 60            # swatch width/height
TX     = X0 + SW + 28  # text x
FONT   = 'CMUBright-Roman, DejaVu Sans, sans-serif'

def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

def row_solid(y, rgb, label, opacity=1.0, stroke=None):
    st = stroke if stroke else rgb
    return (f'<rect x="{X0}" y="{y}" width="{SW}" height="{SW}" rx="6" '
            f'style="fill:rgb({rgb});fill-opacity:{opacity};stroke:rgb({st});stroke-width:2"/>'
            f'<text x="{TX}" y="{y+SW*0.72:.0f}" font-size="34" font-family="{FONT}" '
            f'style="fill:rgb(20,20,20)">{esc(label)}</text>')

def row_gradient(y, gid, label):
    return (f'<rect x="{X0}" y="{y}" width="{SW}" height="{SW}" rx="6" '
            f'style="fill:url(#{gid});stroke:rgb(120,120,120);stroke-width:1"/>'
            f'<text x="{TX}" y="{y+SW*0.72:.0f}" font-size="34" font-family="{FONT}" '
            f'style="fill:rgb(20,20,20)">{esc(label)}</text>')

def row_line(y, rgb, width, label):
    cy = y + SW/2
    return (f'<line x1="{X0}" y1="{cy:.0f}" x2="{X0+SW}" y2="{cy:.0f}" '
            f'style="stroke:rgb({rgb});stroke-width:{width};stroke-linecap:round"/>'
            f'<text x="{TX}" y="{y+SW*0.72:.0f}" font-size="34" font-family="{FONT}" '
            f'style="fill:rgb(20,20,20)">{esc(label)}</text>')

def row_curve(y, rgb, label):
    # a shallow bezier to read as a "link"
    x1, x2 = X0, X0+SW
    cy = y + SW/2
    d = f'M{x1},{cy+14:.0f} Q{(x1+x2)/2:.0f},{cy-24:.0f} {x2},{cy+14:.0f}'
    return (f'<path d="{d}" style="fill:none;stroke:rgb({rgb});stroke-width:5;'
            f'stroke-opacity:0.85;stroke-linecap:round"/>'
            f'<text x="{TX}" y="{y+SW*0.72:.0f}" font-size="34" font-family="{FONT}" '
            f'style="fill:rgb(20,20,20)">{esc(label)}</text>')

def header(y, text, size=40):
    return (f'<text x="{X0}" y="{y}" font-size="{size}" font-family="{FONT}" '
            f'style="fill:rgb(0,0,0);font-weight:bold">{esc(text)}</text>')

def subhead(y, text):
    return (f'<text x="{X0}" y="{y}" font-size="34" font-family="{FONT}" '
            f'style="fill:rgb(70,70,70);font-weight:bold">{esc(text)}</text>')

def build_legend():
    parts = ['<g id="legend-v2">']
    # panel background
    parts.append(f'<rect x="{PLOT_W}" y="0" width="{PANEL}" height="{PLOT_W}" style="fill:rgb(255,255,255)"/>')
    parts.append(f'<line x1="{PLOT_W}" y1="0" x2="{PLOT_W}" y2="{PLOT_W}" style="stroke:rgb(210,210,210);stroke-width:2"/>')
    # gradient def for the gene-density heatmap swatch (blues-9: low -> high)
    parts.append(
        '<defs><linearGradient id="lg-genedens" x1="0" y1="0" x2="1" y2="0">'
        '<stop offset="0" stop-color="rgb(247,251,255)"/>'
        '<stop offset="1" stop-color="rgb(8,48,107)"/>'
        '</linearGradient></defs>')

    parts.append(header(150, 'circos_v2 — track legend', 46))
    parts.append(f'<text x="{X0}" y="205" font-size="32" font-family="{FONT}" '
                 f'style="fill:rgb(90,90,90);font-style:italic">rings, outer → inner</text>')

    y = 300
    STEP = 108
    rows = [
        ('solid',  '218,165,32',  'Telomere caps (TRC2/TRC13)', 1.0, '160,120,10'),
        ('solid',  '224,60,132',  'New sequence (vs mHetGlaV3)', 1.0, None),
        ('solid',  '251,106,74',  'Centromere region', 0.33, '251,106,74'),
        ('grad',   'lg-genedens', 'Gene density (low → high)'),
        ('line',   '0,128,128', 5, 'TE / repeat density'),
        ('line',   '0,0,0', 5,    'GC content'),
        ('solid',  '150,150,150', 'Heterozygosity (histogram)', 0.5, '150,150,150'),
    ]
    for r in rows:
        if r[0] == 'solid':
            parts.append(row_solid(y, r[1], r[2], r[3], r[4]))
        elif r[0] == 'grad':
            parts.append(row_gradient(y, r[1], r[2]))
        elif r[0] == 'line':
            parts.append(row_line(y, r[1], r[2], r[3]))
        y += STEP

    y += 18
    parts.append(subhead(y, 'Satellite families (by role)')); y += 68
    for rgb, lab in [('136,23,152', 'Centromeric (TRC1/3/5)'),
                     ('150,75,20',  'Peritelomeric (TRC2/13)'),
                     ('34,139,34',  'rRNA large subunit (TRC7)'),
                     ('190,30,140', '5S rDNA (TRC10)')]:
        parts.append(row_solid(y, rgb, lab)); y += STEP

    y += 18
    parts.append(row_solid(y, '251,106,74', 'CENPA peaks (histogram)', 0.5, '251,106,74')); y += STEP

    y += 18
    parts.append(subhead(y, 'Segmental duplications (≥100 kb)')); y += 68
    parts.append(row_curve(y, '107,174,214', 'Intra-chromosomal')); y += STEP
    parts.append(row_curve(y, '253,141,60',  'Inter-chromosomal')); y += STEP

    parts.append('</g>')
    return '\n'.join(parts)

def main():
    svg = open(SVG).read()
    # idempotent: strip a legend we added before
    svg = re.sub(r'<g id="legend-v2">.*?</g>\s*(?=</svg>)', '', svg, flags=re.S)
    # widen canvas (always to NEW_W so re-runs are stable)
    svg = re.sub(r'(<svg\s+width=")\d+(px")', rf'\g<1>{NEW_W}\g<2>', svg, count=1)
    # inject legend before </svg>
    svg = svg.replace('</svg>', build_legend() + '\n</svg>')
    open(SVG, 'w').write(svg)
    print(f'legend injected -> {SVG} (canvas {NEW_W}x{PLOT_W})')
    subprocess.run([RSVG, '-w', str(NEW_W), '-h', str(PLOT_W), SVG, '-o', PNG], check=True)
    print(f'rendered -> {PNG}')

if __name__ == '__main__':
    main()
