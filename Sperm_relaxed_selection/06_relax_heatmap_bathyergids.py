#!/usr/bin/env python3
"""
06_relax_heatmap_bathyergids.py

Creates a publication-quality heatmap of RELAX K values for the three
bathyergid species (HetGla, BatSui, FukAns) tested against all others.

Reads RELAX JSON outputs directly. Shows:
  - Color scale: blue (relaxed, K<1) through white (K=1) to red (intensified, K>1)
  - Significance markers: * p<0.05, ** p<0.01, *** p<0.001
  - Artifact annotations: K~0 and K>5 shown with text labels
  - Gene absence marked with X

Usage:
  python3 06_relax_heatmap_bathyergids.py
  python3 06_relax_heatmap_bathyergids.py --output figures/relax_heatmap_bath.pdf
"""

import os
import sys
import json
import glob
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch

# =============================================================================
# CONFIGURATION
# =============================================================================

RELAX_DIR = 'results/09_relax/targets17'

# Schemes to show (in display order left to right). v4 directory names:
SCHEMES = ['nmr_vs_nonmolerat', 'bats_vs_nonmolerat', 'amr_vs_nonmolerat',
           'social_vs_nonmolerat', 'allmolerat_vs_nonmolerat']

SCHEME_LABELS = {
    'nmr_vs_nonmolerat':        'HetGla\n(eusocial)',
    'bats_vs_nonmolerat':       'BatSui\n(solitary)',
    'amr_vs_nonmolerat':        'FukAns\n(cooperative)',
    'social_vs_nonmolerat':     'HetGla +\nFukAns',
    'allmolerat_vs_nonmolerat': 'All\nBathyergidae',
}

# All 17 target genes
TARGET_GENES = [
    "ADIG", "ATP1A4", "CABS1", "DNHD1", "DPY19L2",
    "FAM71B", "FAM71F2", "FSIP2", "MROH2B", "NDC80",
    "NPM2", "RIMBP3", "SPATA33", "SPEM3", "TMPRSS12",
    "TTLL8", "ZAN",
]

# Thresholds (K<=0.01 lower boundary matches the settled atlas thresholds; a
# boundary MLE invalidates the chi-square LRT, so these cells are flagged, not
# plotted on the K colour scale)
K_ZERO_THRESHOLD = 0.01
K_HIGH_THRESHOLD = 5.0
ALPHA = 0.05

# files on disk carry current HGNC symbols; TARGET_GENES uses the legacy ones
DISPLAY_NAME = {'GARIN3': 'FAM71B', 'GARIN1A': 'FAM71F2'}

# files on disk carry current HGNC symbols; TARGET_GENES uses the legacy ones
DISPLAY_NAME = {'GARIN3': 'FAM71B', 'GARIN1A': 'FAM71F2'}


# =============================================================================
# FUNCTIONS
# =============================================================================

def bh_fdr(p):
    # Benjamini-Hochberg q-values
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = p[o] * n / np.arange(1, n + 1)
    q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.minimum(q, 1.0)


def bh_fdr(p):
    # Benjamini-Hochberg q-values
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = p[o] * n / np.arange(1, n + 1)
    q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.minimum(q, 1.0)


def parse_relax_json(json_file):
    """Parse a single RELAX JSON file."""
    try:
        with open(json_file) as f:
            data = json.load(f)
        tr = data.get('test results', {})
        k = tr.get('relaxation or intensification parameter')
        p = tr.get('p-value')
        return k, p
    except Exception:
        return None, None


def load_scheme(scheme, relax_dir=RELAX_DIR):
    """Load all RELAX results for a scheme. Returns {gene: (K, p)}."""
    scheme_dir = os.path.join(relax_dir, scheme)
    results = {}

    if not os.path.exists(scheme_dir):
        return results

    for jf in glob.glob(os.path.join(scheme_dir, 'json', '*_RELAX.json')):
        gene = os.path.basename(jf).replace('_RELAX.json', '')
        gene = DISPLAY_NAME.get(gene, gene)
        k, p = parse_relax_json(jf)
        results[gene] = (k, p)

    return results


def classify(k, p):
    """Classify a result. Returns category string.
    Boundary thresholds (K<=0.01, K>=5) match the atlas rebuild convention."""
    if k is None or p is None:
        return 'absent'
    if k <= K_ZERO_THRESHOLD:
        return 'artifact_zero'
    if k >= K_HIGH_THRESHOLD:
        return 'artifact_high'
    if p < ALPHA and k < 1:
        return 'relaxed'
    if p < ALPHA and k > 1:
        return 'intensified'
    return 'not_significant'


def sig_stars(p):
    """Return significance stars."""
    if p is None:
        return ''
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return ''


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="RELAX heatmap for bathyergid species"
    )
    parser.add_argument("--output",
                        default="results/11_figures/targets17/relax_heatmap_bathyergids",
                        help="Output file path (without extension)")
    parser.add_argument("--relax_dir", default=RELAX_DIR,
                        help="RELAX results directory")
    parser.add_argument("--fdr", action="store_true",
                        help="add a dagger to cells surviving Benjamini-Hochberg "
                             "correction (q<0.05, computed within scheme over "
                             "testable cells); off by default, so the figure is "
                             "unchanged unless asked for")
    args = parser.parse_args()

    relax_dir = args.relax_dir

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)

    # Load data
    print("Loading RELAX results...")
    all_data = {}
    for scheme in SCHEMES:
        all_data[scheme] = load_scheme(scheme, relax_dir)
        print(f"  {scheme}: {len(all_data[scheme])} genes")

    # Build matrices
    n_genes = len(TARGET_GENES)
    n_schemes = len(SCHEMES)

    k_matrix = np.full((n_genes, n_schemes), np.nan)
    p_matrix = np.full((n_genes, n_schemes), np.nan)
    cat_matrix = [['absent'] * n_schemes for _ in range(n_genes)]

    for j, scheme in enumerate(SCHEMES):
        for i, gene in enumerate(TARGET_GENES):
            if gene in all_data[scheme]:
                k, p = all_data[scheme][gene]
                if k is not None:
                    k_matrix[i, j] = k
                if p is not None:
                    p_matrix[i, j] = p
                cat_matrix[i][j] = classify(k, p)

    # BH q-values within each scheme, over testable cells only (boundary
    # estimates are excluded from the denominator as well as from the test)
    q_matrix = np.full((n_genes, n_schemes), np.nan)
    for j in range(n_schemes):
        idx = [i for i in range(n_genes)
               if cat_matrix[i][j] in ('relaxed', 'intensified', 'not_significant')]
        if idx:
            q_matrix[idx, j] = bh_fdr(p_matrix[idx, j])

    # BH q-values within each scheme, over testable cells only (boundary
    # estimates are excluded from the denominator as well as from the test)
    q_matrix = np.full((n_genes, n_schemes), np.nan)
    for j in range(n_schemes):
        idx = [i for i in range(n_genes)
               if cat_matrix[i][j] in ('relaxed', 'intensified', 'not_significant')]
        if idx:
            q_matrix[idx, j] = bh_fdr(p_matrix[idx, j])

    # Sort genes: by number of schemes with evidence (relaxed or artifact),
    # then by mean K (lowest first)
    def sort_key(idx):
        gene = TARGET_GENES[idx]
        n_evidence = sum(1 for j in range(n_schemes)
                         if cat_matrix[idx][j] in ('relaxed', 'artifact_zero',
                                                     'artifact_high'))
        n_absent = sum(1 for j in range(n_schemes)
                       if cat_matrix[idx][j] == 'absent')
        # Mean K for sorting (use 0.001 for artifacts, 10 for absent)
        k_vals = []
        for j in range(n_schemes):
            k = k_matrix[idx, j]
            if np.isnan(k):
                k_vals.append(10)  # absent goes to bottom
            elif k < K_ZERO_THRESHOLD:
                k_vals.append(0.001)
            elif k > K_HIGH_THRESHOLD:
                k_vals.append(0.001)  # artifacts = extreme relaxation
            else:
                k_vals.append(k)
        mean_k = np.mean(k_vals)
        # Sort: most evidence first, then lowest K, then fewest absent
        return (-n_evidence, mean_k, -n_absent)

    sort_indices = sorted(range(n_genes), key=sort_key)
    gene_order = [TARGET_GENES[i] for i in sort_indices]
    k_sorted = k_matrix[sort_indices, :]
    p_sorted = p_matrix[sort_indices, :]
    q_sorted = q_matrix[sort_indices, :]
    cat_sorted = [cat_matrix[i] for i in sort_indices]

    # =========================================================================
    # PLOT
    # =========================================================================

    fig_h = max(5, n_genes * 0.38 + 2.5)
    fig_w = 7.5
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Color matrix: use log2(K) for display, clip to [-3, 3]
    display_matrix = np.full_like(k_sorted, np.nan)
    for i in range(n_genes):
        for j in range(n_schemes):
            k = k_sorted[i, j]
            cat = cat_sorted[i][j]
            if cat == 'absent':
                display_matrix[i, j] = np.nan
            elif cat == 'artifact_zero':
                display_matrix[i, j] = np.nan  # will be overdrawn with purple
            elif cat == 'artifact_high':
                display_matrix[i, j] = np.nan  # will be overdrawn with orange
            elif not np.isnan(k) and k > 0:
                display_matrix[i, j] = np.clip(np.log2(k), -3, 3)
            else:
                display_matrix[i, j] = 0

    # Colormap: blue (relaxed) - white (neutral) - red (intensified)
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad(color='#F0F0F0')  # light grey for absent

    norm = TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)

    im = ax.imshow(display_matrix, aspect='auto', cmap=cmap, norm=norm,
                    interpolation='nearest')

    # Overdraw artifact cells with distinct colors
    from matplotlib.patches import Rectangle
    for i in range(len(gene_order)):
        for j in range(n_schemes):
            cat = cat_sorted[i][j]
            if cat == 'artifact_high':
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                             facecolor='#FFA500', edgecolor='white',
                             linewidth=1, zorder=2))
            elif cat == 'artifact_zero':
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1,
                             facecolor='#9932CC', edgecolor='white',
                             linewidth=1, zorder=2))

    # Cell annotations
    for i in range(len(gene_order)):
        for j in range(n_schemes):
            cat = cat_sorted[i][j]
            k = k_sorted[i, j]
            p = p_sorted[i, j]

            if cat == 'absent':
                ax.text(j, i, 'X', ha='center', va='center',
                        fontsize=11, fontweight='bold', color='#999999')

            elif cat == 'artifact_zero':
                ax.text(j, i, 'K~0', ha='center', va='center',
                        fontsize=7, fontweight='bold', color='white')

            elif cat == 'artifact_high':
                ax.text(j, i, f'K={k:.0f}', ha='center', va='center',
                        fontsize=7, fontweight='bold', color='white')

            elif not np.isnan(p):
                stars = sig_stars(p)
                if args.fdr and not np.isnan(q_sorted[i, j]) \
                        and q_sorted[i, j] < ALPHA:
                    stars += '\u2020'
                if stars:
                    # Determine text color for contrast
                    log2k = np.log2(k) if k > 0 else 0
                    text_color = 'white' if abs(log2k) > 1.5 else 'black'
                    ax.text(j, i, stars, ha='center', va='center',
                            fontsize=10, fontweight='bold', color=text_color)

    # Grid lines
    for i in range(len(gene_order) + 1):
        ax.axhline(y=i - 0.5, color='white', linewidth=1)
    for j in range(n_schemes + 1):
        ax.axvline(x=j - 0.5, color='white', linewidth=1)

    # Thicker separator between individual species and group schemes
    ax.axvline(x=2.5, color='black', linewidth=2.0, alpha=0.6)

    # Axes
    ax.set_xticks(range(n_schemes))
    ax.set_xticklabels([SCHEME_LABELS[s] for s in SCHEMES],
                        fontsize=9, fontweight='bold')
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')

    ax.set_yticks(range(len(gene_order)))
    ax.set_yticklabels(gene_order, fontsize=10, fontweight='bold',
                        fontstyle='italic')

    # Column group labels above species labels
    fig.text(0.37, 0.89, 'Individual species', ha='center', fontsize=8,
             fontstyle='italic', color='grey')
    fig.text(0.72, 0.89, 'Combined', ha='center', fontsize=8,
             fontstyle='italic', color='grey')

    # Title
    ax.set_title('RELAX Selection Analysis\nBathyergidae vs All Other Species',
                 fontsize=12, fontweight='bold', pad=55)

    # Colorbar - placed below the heatmap
    cbar_ax = fig.add_axes([0.15, 0.04, 0.55, 0.02])  # [left, bottom, width, height]
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_ticks([-3, -2, -1, 0, 1, 2, 3])
    cbar.set_ticklabels(['0.12', '0.25', '0.5', '1', '2', '4', '8'])
    cbar.set_label('K value', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # Legend - placed to the right of colorbar
    legend_elements = [
        Patch(facecolor='#2166AC', edgecolor='grey', label='Relaxed (K<1, p<0.05)'),
        Patch(facecolor='#B2182B', edgecolor='grey', label='Intensified (K>1, p<0.05)'),
        Patch(facecolor='#9932CC', edgecolor='grey', label='K~0 artifact'),
        Patch(facecolor='#FFA500', edgecolor='grey', label='K>5 artifact'),
    ]

    fig.legend(handles=legend_elements, loc='lower right',
               bbox_to_anchor=(0.98, 0.005), fontsize=7.5,
               ncol=1, frameon=True, edgecolor='grey',
               fancybox=False)

    # Marker legend as text
    _marker_txt = ('*** p<0.001   ** p<0.01   * p<0.05   X = gene absent   '
                   'K~0 / K=N = optimiser artifact')
    if args.fdr:
        _marker_txt = ('*** p<0.001   ** p<0.01   * p<0.05   \u2020 BH q<0.05   '
                       'X = gene not tested   K~0 / K=N = optimiser artifact')
    fig.text(0.15, 0.13, _marker_txt,
             fontsize=8, color='grey', va='top')

    plt.subplots_adjust(bottom=0.16, top=0.82, left=0.18, right=0.95)

    # Save
    plt.savefig(args.output + '.png', dpi=600, bbox_inches='tight',
                facecolor='white')
    plt.savefig(args.output + '.pdf', bbox_inches='tight',
                facecolor='white')
    plt.close()

    print(f"Saved: {args.output}.png / .pdf")

    # Print summary
    print("\nSummary:")
    for j, scheme in enumerate(SCHEMES):
        n_relax = sum(1 for i in range(len(gene_order))
                      if cat_sorted[i][j] == 'relaxed')
        n_art = sum(1 for i in range(len(gene_order))
                    if cat_sorted[i][j] in ('artifact_zero', 'artifact_high'))
        n_absent = sum(1 for i in range(len(gene_order))
                       if cat_sorted[i][j] == 'absent')
        print(f"  {SCHEME_LABELS[scheme].replace(chr(10), ' ')}: "
              f"{n_relax} relaxed, {n_art} artifacts, {n_absent} absent")


if __name__ == "__main__":
    main()
