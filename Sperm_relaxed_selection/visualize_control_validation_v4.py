#!/usr/bin/env python3
"""
visualize_control_validation_v4.py

Comprehensive visualization of control gene validation results (v21).

Updates from v3:
- Matches v21 gene set: removed GULO, UOX, TAS1R2, CMAH (never in OrthoFinder)
- Final set: 8 positive + 9 negative = 17 control genes

Updates from v2:
- Reflects v19 control gene set: PRM1/PRM2 removed, 5 Eisenberg & Levanon 2013
  housekeeping genes added (VCP, GPI, REEP5, VPS29, SNRPD3)
- 8 positive controls (CRYAA, CRYAB, CRYBB1, RHO, GRK1, OPN1SW, OR1A1, OR2W1)
- 9 negative controls (ATP5F1A, RPL13, RPS13, OPN1LW, VCP, GPI, REEP5, VPS29,
  SNRPD3)

Creates:
1. Paired bar plot - General vs Specific modes side-by-side per gene
2. Forest plot - K values with expected ranges, positive vs negative grouping
3. Per-mode bar plots with pass/fail/artifact classification
4. Text summary report

Usage:
    python3 visualize_control_validation_v4.py [results.tsv]

Default input: results/control_gene_validation_v21/validation_results.tsv
"""

import os
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
from pathlib import Path

# =============================================================================
# CONFIGURATION
# =============================================================================

DEFAULT_RESULTS = "results/09_relax/validation/general/relax_results.tsv"
K_ZERO_THRESHOLD = 0.0001
K_UNRELIABLE_THRESHOLD = 5.0

# --- v19 positive controls (8) ---
ALL_POSITIVE_GENES = {
    "CRYAA":  {"desc": "Lens crystallin αA (subterranean bathyergids)"},
    "CRYAB":  {"desc": "Lens crystallin αB (subterranean bathyergids)"},
    "CRYBB1": {"desc": "Lens crystallin βB1 (subterranean bathyergids)"},
    "RHO":    {"desc": "Rhodopsin (subterranean)"},
    "GRK1":   {"desc": "Rhodopsin kinase (subterranean)"},
    "OPN1SW": {"desc": "Blue opsin (subterranean)"},
    "OR1A1":  {"desc": "Olfactory receptor (catarrhine primates)"},
    "OR2W1":  {"desc": "Olfactory receptor (catarrhine primates)"},
}

# --- v19 negative controls (9) ---
# PRM1/PRM2 removed (Wyckoff et al. 2000 shows positive selection)
# Added VCP, GPI, REEP5, VPS29, SNRPD3 (Eisenberg & Levanon 2013)
ALL_NEGATIVE_GENES = {
    "ATP5F1A": {"desc": "ATP synthase F1α — housekeeping"},
    "RPL13":   {"desc": "Ribosomal protein L13 — housekeeping"},
    "RPS13":   {"desc": "Ribosomal protein S13 — housekeeping"},
    "OPN1LW":  {"desc": "Red/green opsin — conserved"},
    "VCP":     {"desc": "Valosin-containing protein — housekeeping"},
    "GPI":     {"desc": "Glucose-6-phosphate isomerase — housekeeping"},
    "REEP5":   {"desc": "Receptor accessory protein 5 — housekeeping"},
    "VPS29":   {"desc": "Vacuolar protein sorting 29 — housekeeping"},
    "SNRPD3":  {"desc": "Small nuclear ribonucleoprotein D3 — housekeeping"},
}

# Display order
GENE_ORDER_POS = ["CRYAA", "CRYAB", "CRYBB1", "RHO", "GRK1", "OPN1SW",
                   "OR1A1", "OR2W1"]
GENE_ORDER_NEG = ["ATP5F1A", "RPL13", "RPS13", "OPN1LW",
                   "VCP", "GPI", "REEP5", "VPS29", "SNRPD3"]

# Colors
COLORS = {
    'pos_pass':   '#2ecc71',   # Green
    'pos_fail':   '#e74c3c',   # Red
    'neg_pass':   '#3498db',   # Blue
    'neg_fail':   '#e67e22',   # Orange
    'artifact':   '#95a5a6',   # Grey
    'error':      '#d5d8dc',   # Light grey
    'missing':    '#f0f0f0',   # Very light grey
    'positive':   '#27ae60',   # Green (for y-axis labels)
    'negative':   '#7f8c8d',   # Grey
}


# =============================================================================
# DATA LOADING AND CLASSIFICATION
# =============================================================================

def load_results(filepath):
    """Load validation results TSV."""
    df = pd.read_csv(filepath, sep='\t')
    # v5 pipeline layout: Gene/K/p_value/LRT/Scheme/Set -> gene/mode
    if 'Gene' in df.columns and 'gene' not in df.columns:
        df = df.rename(columns={'Gene': 'gene', 'Scheme': 'mode'})
    if 'mode' not in df.columns:
        df['mode'] = 'general'
    # Extract gene name from gene_og if needed
    if 'gene' not in df.columns and 'gene_og' in df.columns:
        df['gene'] = df['gene_og'].str.extract(r'^([A-Z0-9]+)_OG')
    # Ensure control_type is populated
    if 'control_type' not in df.columns or df['control_type'].isna().all():
        pos_genes = set(ALL_POSITIVE_GENES.keys())
        neg_genes = set(ALL_NEGATIVE_GENES.keys())
        df['control_type'] = df['gene'].apply(
            lambda g: 'positive' if g in pos_genes
            else ('negative' if g in neg_genes else 'unknown'))
    return df


def classify_result(row):
    """
    Classify a result as pass/fail/artifact/error.

    Returns one of:
        'pass', 'fail', 'error', 'k_zero_artifact', 'k_high_artifact'
    """
    k = row['K']
    # Use FDR q-value for the significance call when present, so control
    # classification matches the atlas (which reports q<0.05). Fall back to
    # raw p only if q is absent.
    sig_val = row['q_value'] if ('q_value' in row and pd.notna(row.get('q_value'))) else row['p_value']
    ctype = row['control_type']

    if pd.isna(k) or pd.isna(sig_val):
        return 'error'
    if k < K_ZERO_THRESHOLD:
        return 'k_zero_artifact'
    if k > K_UNRELIABLE_THRESHOLD:
        return 'k_high_artifact'

    if ctype == 'positive':
        return 'pass' if (k < 1 and sig_val < 0.05) else 'fail'
    else:  # negative
        return 'pass' if (sig_val >= 0.05 or (0.5 <= k <= 2.0)) else 'fail'


def get_color(row, classification):
    """Get bar color based on control type and classification."""
    if classification in ('k_zero_artifact', 'k_high_artifact'):
        return COLORS['artifact']
    if classification == 'error':
        return COLORS['error']
    if row['control_type'] == 'positive':
        return COLORS['pos_pass'] if classification == 'pass' else COLORS['pos_fail']
    else:
        return COLORS['neg_pass'] if classification == 'pass' else COLORS['neg_fail']


def significance_label(p):
    """Return significance asterisks."""
    if pd.isna(p):
        return ''
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return ''


# =============================================================================
# FIGURE 1: Paired bar plot (General vs Specific per gene)
# =============================================================================

def create_paired_barplot(df, output_dir):
    """Side-by-side bars for general vs specific mode, grouped by gene."""

    genes_in_data = set(df['gene'].unique())

    # Build ordered gene list: positives then negatives (only those present)
    gene_order = [g for g in GENE_ORDER_POS if g in genes_in_data]
    gene_order += [g for g in GENE_ORDER_NEG if g in genes_in_data]

    n_genes = len(gene_order)
    if n_genes == 0:
        print("  No genes to plot for paired barplot")
        return

    fig, ax = plt.subplots(figsize=(max(14, n_genes * 0.85), 7))

    bar_width = 0.35
    x = np.arange(n_genes)

    for i, mode in enumerate(['general', 'specific']):
        mode_data = df[df['mode'] == mode]

        ks = []
        colors = []
        hatches = []
        edge_styles = []

        for gene in gene_order:
            row_data = mode_data[mode_data['gene'] == gene]
            if len(row_data) == 0:
                ks.append(0)
                colors.append(COLORS['missing'])
                hatches.append('//')
                edge_styles.append('--')
            else:
                row = row_data.iloc[0]
                k = row['K'] if not pd.isna(row['K']) else 0
                classification = classify_result(row)

                # Cap K for display
                if k > K_UNRELIABLE_THRESHOLD:
                    k = K_UNRELIABLE_THRESHOLD

                ks.append(k)
                colors.append(get_color(row, classification))
                if classification in ('k_zero_artifact', 'k_high_artifact'):
                    hatches.append('xx')
                elif classification == 'error':
                    hatches.append('...')
                else:
                    hatches.append('')
                edge_styles.append('-')

        offset = -bar_width / 2 + i * bar_width
        bars = ax.bar(x + offset, ks, bar_width, color=colors,
                      edgecolor='black', linewidth=0.5)

        # Apply hatching
        for bar, hatch in zip(bars, hatches):
            bar.set_hatch(hatch)

        # Add K value labels
        for j, (bar, gene) in enumerate(zip(bars, gene_order)):
            row_data = mode_data[mode_data['gene'] == gene]
            if len(row_data) > 0:
                row = row_data.iloc[0]
                k = row['K']
                p = row['p_value']
                if pd.isna(k):
                    # Mark errors
                    ax.annotate('ERR', xy=(bar.get_x() + bar.get_width() / 2, 0.02),
                                ha='center', va='bottom', fontsize=5.5,
                                rotation=45, color='red', fontweight='bold')
                else:
                    sig = significance_label(p)
                    height = min(k, K_UNRELIABLE_THRESHOLD)
                    if k > K_UNRELIABLE_THRESHOLD:
                        label = f'>{K_UNRELIABLE_THRESHOLD:.0f}{sig}'
                    elif k < 0.01:
                        label = f'{k:.4f}{sig}'
                    else:
                        label = f'{k:.2f}{sig}'
                    ax.annotate(label,
                                xy=(bar.get_x() + bar.get_width() / 2, height),
                                xytext=(0, 3), textcoords='offset points',
                                ha='center', va='bottom', fontsize=6, rotation=45)

    # Reference line
    ax.axhline(y=1.0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
    ax.axhspan(0, 1, alpha=0.05, color='green')

    # Separator between positive and negative controls
    n_pos = len([g for g in gene_order if g in ALL_POSITIVE_GENES])
    if n_pos < n_genes:
        ax.axvline(x=n_pos - 0.5, color='grey', linestyle=':', linewidth=1, alpha=0.6)

    # X-axis labels
    ax.set_xticks(x)
    ax.set_xticklabels(gene_order, rotation=45, ha='right', fontsize=10, fontweight='bold')
    ax.set_ylabel('K (relaxation/intensification parameter)', fontsize=12)
    ax.set_title('Control Gene Validation: General vs Specific Modes',
                 fontsize=13, fontweight='bold')

    # Cap y-axis
    max_k = df['K'].dropna()
    max_k = max_k[max_k <= K_UNRELIABLE_THRESHOLD]
    if len(max_k) > 0:
        ax.set_ylim(0, min(max_k.max() * 1.3, 3.0))
    else:
        ax.set_ylim(0, 2.0)

    # Section labels (placed AFTER ylim is set, using axes transform)
    if n_pos > 0 and n_pos < n_genes:
        # Convert gene positions to axes coordinates for label placement
        # Positive controls center
        pos_center_data = (n_pos / 2) - 0.5
        neg_center_data = (n_pos + n_genes) / 2 - 0.5
        # Use transform to mix data x with axes y
        trans = ax.get_xaxis_transform()
        label_bbox = dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7, edgecolor='none')
        ax.text(pos_center_data, 0.96, 'POSITIVE CONTROLS',
                ha='center', fontsize=9, fontstyle='italic', color='grey',
                transform=trans, bbox=label_bbox, zorder=10)
        ax.text(neg_center_data, 0.96, 'NEGATIVE CONTROLS',
                ha='center', fontsize=9, fontstyle='italic', color='grey',
                transform=trans, bbox=label_bbox, zorder=10)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COLORS['pos_pass'], edgecolor='k', label='Positive PASS'),
        mpatches.Patch(facecolor=COLORS['pos_fail'], edgecolor='k', label='Positive FAIL'),
        mpatches.Patch(facecolor=COLORS['neg_pass'], edgecolor='k', label='Negative PASS'),
        mpatches.Patch(facecolor=COLORS['neg_fail'], edgecolor='k', label='Negative FAIL'),
        mpatches.Patch(facecolor=COLORS['artifact'], edgecolor='k', hatch='xx',
                       label='K artifact (≈0 or >5)'),
        mpatches.Patch(facecolor=COLORS['error'], edgecolor='k', hatch='...',
                       label='RELAX error (N/A)'),
        mlines.Line2D([], [], color='black', linestyle='--', label='K=1 (neutral)'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=7.5, ncol=2,
              framealpha=0.9, edgecolor='grey')

    # Mode indicator
    # Manual mode legend (bar colours are mixed, so indicate with text)
    ax.text(0.01, 0.98,
            'Left bar = General mode\nRight bar = Specific mode\n\n'
            'K < 1: Relaxed\nK > 1: Intensified\n'
            '* p<0.05  ** p<0.01  *** p<0.001',
            transform=ax.transAxes, fontsize=7, va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.yaxis.grid(True, linestyle=':', alpha=0.4)
    ax.set_axisbelow(True)

    plt.tight_layout()
    outfile = os.path.join(output_dir, 'validation_paired_barplot.png')
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.savefig(outfile.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"  Saved: {outfile}")
    plt.close()


# =============================================================================
# FIGURE 2: Forest plot
# =============================================================================

def create_forest_plot(df, output_dir):
    """Forest plot showing K values with expected ranges, positive vs negative grouping."""

    genes_in_data = set(df['gene'].unique())

    # Build gene list (only genes present in data)
    plot_genes = []
    for g in GENE_ORDER_POS:
        if g in genes_in_data:
            desc = ALL_POSITIVE_GENES[g]['desc']
            plot_genes.append((g, 'positive', desc))
    for g in GENE_ORDER_NEG:
        if g in genes_in_data:
            desc = ALL_NEGATIVE_GENES[g]['desc']
            plot_genes.append((g, 'negative', desc))

    n_genes = len(plot_genes)
    if n_genes == 0:
        print("  No genes for forest plot")
        return

    fig, ax = plt.subplots(figsize=(11, max(6, n_genes * 0.45)))

    y_positions = np.arange(n_genes)[::-1]  # top to bottom

    # Expected ranges
    expected_pos = (0, 1.0)    # Positive: K < 1
    expected_neg = (0.5, 2.0)  # Negative: K ≈ 1

    # Draw expected range bands
    for i, (gene, ctype, desc) in enumerate(plot_genes):
        y = y_positions[i]
        if ctype == 'positive':
            ax.barh(y, expected_pos[1] - expected_pos[0], left=expected_pos[0],
                    height=0.7, color='#d5f5e3', alpha=0.5, zorder=0)
        else:
            ax.barh(y, expected_neg[1] - expected_neg[0], left=expected_neg[0],
                    height=0.7, color='#d6eaf8', alpha=0.5, zorder=0)

    # Plot K values for each mode
    markers = {'general': 'o', 'specific': 's'}
    offsets = {'general': 0.12, 'specific': -0.12}
    mode_colors = {'general': '#2c3e50', 'specific': '#c0392b'}

    for mode in ['general', 'specific']:
        mode_data = df[df['mode'] == mode]

        for i, (gene, ctype, desc) in enumerate(plot_genes):
            y = y_positions[i] + offsets[mode]
            row_data = mode_data[mode_data['gene'] == gene]

            if len(row_data) == 0:
                continue

            row = row_data.iloc[0]
            k = row['K']
            p = row['p_value']

            if pd.isna(k):
                # Mark errors with an X
                ax.plot(0.5, y, marker='x', markersize=8,
                        color=mode_colors[mode], markeredgewidth=2,
                        alpha=0.5, zorder=5)
                ax.text(0.7, y, 'ERR', fontsize=6, va='center',
                        color=mode_colors[mode], alpha=0.5)
                continue

            classification = classify_result(row)

            # Determine marker properties
            if classification in ('k_zero_artifact', 'k_high_artifact'):
                mfc = COLORS['artifact']
                marker_size = 6
                alpha = 0.5
                k_display = min(max(k, 0.001), 4.5)
            else:
                mfc = mode_colors[mode]
                marker_size = 8
                alpha = 1.0
                k_display = min(k, 4.5)

            # Filled if significant, open if not
            mec = mfc
            if not pd.isna(p) and p >= 0.05:
                mfc = 'white'  # Open marker for non-significant

            ax.plot(k_display, y, marker=markers[mode], markersize=marker_size,
                    markerfacecolor=mfc, markeredgecolor=mec, markeredgewidth=1.5,
                    alpha=alpha, zorder=5)

            # Arrow for off-scale K values
            if k > 4.5:
                ax.annotate('', xy=(4.7, y), xytext=(4.5, y),
                           arrowprops=dict(arrowstyle='->', color=mec, lw=1.5))
                ax.text(4.75, y, f'K={k:.1f}', fontsize=7, va='center', color=mec)

    # Reference line at K=1
    ax.axvline(x=1.0, color='black', linestyle='--', linewidth=1.5, alpha=0.6)

    # Y-axis labels
    y_labels = []
    y_label_colors = []
    for gene, ctype, desc in plot_genes:
        if ctype == 'positive':
            label = f"{gene}"
            color = COLORS['positive']
        else:
            label = f"{gene}"
            color = COLORS['negative']
        y_labels.append(label)
        y_label_colors.append(color)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=10, fontweight='bold')
    for tick_label, color in zip(ax.get_yticklabels(), y_label_colors):
        tick_label.set_color(color)

    # Add descriptions on right side
    for i, (gene, ctype, desc) in enumerate(plot_genes):
        y = y_positions[i]
        ax.text(5.0, y, desc, fontsize=6.5, va='center', color='grey', style='italic')

    ax.set_xlabel('K (relaxation/intensification parameter)', fontsize=12)
    ax.set_title('Control Gene Validation — Forest Plot', fontsize=13, fontweight='bold')
    ax.set_xlim(-0.1, 9.0)

    # Separator between positive and negative controls
    prev_ctype = None
    for i, (gene, ctype, desc) in enumerate(plot_genes):
        if prev_ctype is not None and ctype != prev_ctype:
            y_mid = (y_positions[i] + y_positions[i - 1]) / 2
            ax.axhline(y=y_mid, color='grey', linestyle=':', linewidth=0.8, alpha=0.6)
        prev_ctype = ctype

    # Legend
    legend_elements = [
        mlines.Line2D([], [], marker='o', color='w', markerfacecolor='#2c3e50',
                       markeredgecolor='#2c3e50', markersize=8,
                       label='General mode (p<0.05)'),
        mlines.Line2D([], [], marker='o', color='w', markerfacecolor='white',
                       markeredgecolor='#2c3e50', markersize=8,
                       label='General mode (p≥0.05)'),
        mlines.Line2D([], [], marker='s', color='w', markerfacecolor='#c0392b',
                       markeredgecolor='#c0392b', markersize=8,
                       label='Specific mode (p<0.05)'),
        mlines.Line2D([], [], marker='s', color='w', markerfacecolor='white',
                       markeredgecolor='#c0392b', markersize=8,
                       label='Specific mode (p≥0.05)'),
        mlines.Line2D([], [], marker='x', color='#2c3e50', markersize=8,
                       markeredgewidth=2, linestyle='None',
                       label='RELAX error (N/A)'),
        mlines.Line2D([], [], marker='o', color='w', markerfacecolor=COLORS['artifact'],
                       markeredgecolor=COLORS['artifact'], markersize=6,
                       label='K artifact', alpha=0.5),
        mpatches.Patch(facecolor='#d5f5e3', alpha=0.5, label='Expected range (positive)'),
        mpatches.Patch(facecolor='#d6eaf8', alpha=0.5, label='Expected range (negative)'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=7.5,
              framealpha=0.9, edgecolor='grey')

    ax.xaxis.grid(True, linestyle=':', alpha=0.3)
    ax.set_axisbelow(True)

    plt.tight_layout()
    outfile = os.path.join(output_dir, 'validation_forest_plot.png')
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.savefig(outfile.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"  Saved: {outfile}")
    plt.close()


# =============================================================================
# FIGURE 3: Per-mode bar plots
# =============================================================================

def create_barplot(df, mode, output_dir):
    """Bar plot for a single mode with pass/fail/artifact classification."""

    data = df[df['mode'] == mode].copy()
    if len(data) == 0:
        print(f"  No data for mode: {mode}")
        return

    # Order genes
    gene_order = [g for g in GENE_ORDER_POS + GENE_ORDER_NEG if g in data['gene'].values]
    data = data.set_index('gene').loc[gene_order].reset_index()

    fig, ax = plt.subplots(figsize=(max(14, len(gene_order) * 0.8), 6))

    x = np.arange(len(data))

    colors = []
    hatches = []
    for _, row in data.iterrows():
        classification = classify_result(row)
        colors.append(get_color(row, classification))
        if classification in ('k_zero_artifact', 'k_high_artifact'):
            hatches.append('xx')
        elif classification == 'error':
            hatches.append('...')
        else:
            hatches.append('')

    # Cap K values for display; errors (NaN) -> 0
    k_display = data['K'].copy()
    k_display = k_display.clip(upper=K_UNRELIABLE_THRESHOLD).fillna(0)

    bars = ax.bar(x, k_display, color=colors, edgecolor='black', linewidth=0.5)
    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)

    ax.axhline(y=1.0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
    ax.axhspan(0, 1, alpha=0.05, color='green')

    # Separator between positive and negative
    n_pos = len([g for g in gene_order if g in ALL_POSITIVE_GENES])
    if n_pos < len(gene_order):
        ax.axvline(x=n_pos - 0.5, color='grey', linestyle=':', linewidth=1)

    # Value labels
    for bar, (_, row) in zip(bars, data.iterrows()):
        k = row['K']
        p = row['p_value']
        if pd.isna(k):
            ax.annotate('ERR', xy=(bar.get_x() + bar.get_width() / 2, 0.02),
                        ha='center', va='bottom', fontsize=7, color='red',
                        fontweight='bold')
            continue

        sig = significance_label(p)
        height = min(k, K_UNRELIABLE_THRESHOLD)
        if k > K_UNRELIABLE_THRESHOLD:
            label = f'>{K_UNRELIABLE_THRESHOLD:.0f}{sig}'
        elif k < 0.01:
            label = f'{k:.4f}{sig}'
        elif k < 1:
            label = f'{k:.3f}{sig}'
        else:
            label = f'{k:.2f}{sig}'
        ax.annotate(label, xy=(bar.get_x() + bar.get_width() / 2, height),
                   xytext=(0, 3), textcoords='offset points',
                   ha='center', va='bottom', fontsize=7.5, rotation=0)

    ax.set_xticks(x)
    ax.set_xticklabels([row['gene'] for _, row in data.iterrows()],
                       rotation=45, ha='right', fontsize=10, fontweight='bold')
    ax.set_ylabel('K parameter', fontsize=12)

    # --- Pass rate calculation (excluding errors and artifacts) ---
    pos = data[data['control_type'] == 'positive']
    neg = data[data['control_type'] == 'negative']

    def count_clean(subset):
        """Count clean results (not error, not artifact)."""
        clean = subset[subset['K'].notna() &
                       (subset['K'] >= K_ZERO_THRESHOLD) &
                       (subset['K'] <= K_UNRELIABLE_THRESHOLD)]
        return clean

    pos_clean = count_clean(pos)
    neg_clean = count_clean(neg)

    pos_pass = len(pos_clean[(pos_clean['K'] < 1) & (pos_clean['p_value'] < 0.05)])
    neg_pass = len(neg_clean[(neg_clean['p_value'] >= 0.05) |
                             ((neg_clean['K'] >= 0.5) & (neg_clean['K'] <= 2.0))])

    n_errors = len(pos[pos['K'].isna()]) + len(neg[neg['K'].isna()])
    n_artifacts = (len(pos[pos['K'].notna() & (pos['K'] < K_ZERO_THRESHOLD)]) +
                   len(pos[pos['K'].notna() & (pos['K'] > K_UNRELIABLE_THRESHOLD)]) +
                   len(neg[neg['K'].notna() & (neg['K'] < K_ZERO_THRESHOLD)]) +
                   len(neg[neg['K'].notna() & (neg['K'] > K_UNRELIABLE_THRESHOLD)]))

    mode_name = "GENERAL (Test vs ALL)" if mode == "general" else "SPECIFIC (Test vs reference)"
    ax.set_title(f'Control Gene Validation — {mode_name}\n',
                 fontsize=11, fontweight='bold')

    max_k_val = k_display.max()
    ax.set_ylim(0, max(min(max_k_val * 1.3, 3.0), 1.8))

    legend_elements = [
        mpatches.Patch(facecolor=COLORS['pos_pass'], edgecolor='k', label='Positive PASS'),
        mpatches.Patch(facecolor=COLORS['pos_fail'], edgecolor='k', label='Positive FAIL'),
        mpatches.Patch(facecolor=COLORS['neg_pass'], edgecolor='k', label='Negative PASS'),
        mpatches.Patch(facecolor=COLORS['neg_fail'], edgecolor='k', label='Negative FAIL'),
        mpatches.Patch(facecolor=COLORS['artifact'], edgecolor='k', hatch='xx',
                       label='K artifact (≈0 or >5)'),
        mpatches.Patch(facecolor=COLORS['error'], edgecolor='k', hatch='...',
                       label='RELAX error'),
        mlines.Line2D([], [], color='black', linestyle='--', label='K=1 (neutral)'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8, ncol=2)

    ax.text(0.01, 0.98,
            'K < 1: Relaxed selection\nK > 1: Intensified selection\n'
            '* p<0.05  ** p<0.01  *** p<0.001',
            transform=ax.transAxes, fontsize=7, va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.yaxis.grid(True, linestyle=':', alpha=0.4)
    ax.set_axisbelow(True)

    plt.tight_layout()
    outfile = os.path.join(output_dir, f'validation_barplot_{mode}.png')
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.savefig(outfile.replace('.png', '.pdf'), bbox_inches='tight')
    print(f"  Saved: {outfile}")
    plt.close()


# =============================================================================
# TEXT REPORT
# =============================================================================

def generate_report(df, output_dir):
    """Generate comprehensive text report."""

    lines = []
    lines.append("=" * 80)
    lines.append("CONTROL GENE VALIDATION REPORT")
    lines.append("=" * 80)
    lines.append("")
    lines.append("8 positive controls + 9 negative controls = 17 genes")
    lines.append("")

    # Missing genes
    found_genes = set(df['gene'].unique())
    all_expected = set(ALL_POSITIVE_GENES.keys()) | set(ALL_NEGATIVE_GENES.keys())
    missing = all_expected - found_genes

    lines.append(f"Expected: {len(all_expected)} genes "
                 f"({len(ALL_POSITIVE_GENES)} positive, {len(ALL_NEGATIVE_GENES)} negative)")
    lines.append(f"Found:    {len(found_genes & all_expected)} genes")
    if missing:
        lines.append(f"Missing:  {len(missing)} genes — {', '.join(sorted(missing))}")
        lines.append("  (Not found in OrthoFinder results)")
    lines.append("")

    # Per-mode results
    for mode in ['general', 'specific']:
        mode_data = df[df['mode'] == mode]
        if len(mode_data) == 0:
            continue

        mode_name = ("GENERAL (Test vs ALL others)" if mode == "general"
                     else "SPECIFIC (Test vs defined reference)")
        lines.append("=" * 80)
        lines.append(f"MODE: {mode_name}")
        lines.append("=" * 80)
        lines.append("")

        # --- Positive controls ---
        lines.append("POSITIVE CONTROLS (expect K < 1, p < 0.05):")
        lines.append("-" * 78)
        lines.append(f"{'Gene':<12} {'K':<12} {'p-value':<14} "
                     f"{'Spp':<6} {'Status':<22} {'Note'}")
        lines.append("-" * 78)

        n_pos_pass = 0
        n_pos_clean = 0
        n_pos_error = 0
        n_pos_artifact = 0

        for gene in GENE_ORDER_POS:
            row_data = mode_data[mode_data['gene'] == gene]

            if len(row_data) == 0:
                lines.append(f"{gene:<12} {'—':<12} {'—':<14} "
                             f"{'—':<6} {'NOT FOUND':<22}")
                continue

            row = row_data.iloc[0]
            k = row['K']
            p = row['p_value']

            # Get species count if available
            spp_col = [c for c in ['n_species', 'species_count'] if c in row.index]
            spp = ''
            if spp_col and not pd.isna(row[spp_col[0]]):
                spp = str(int(row[spp_col[0]]))

            classification = classify_result(row)
            k_str = f"{k:.4f}" if not pd.isna(k) else "N/A"
            p_str = f"{p:.2e}" if not pd.isna(p) else "N/A"

            if classification == 'error':
                status = "? ERROR"
                n_pos_error += 1
                note = "RELAX failed"
            elif classification == 'k_zero_artifact':
                status = "⚠ K≈0 ARTIFACT"
                n_pos_artifact += 1
                note = "Numerical boundary"
            elif classification == 'k_high_artifact':
                status = f"⚠ K>{K_UNRELIABLE_THRESHOLD:.0f} ARTIFACT"
                n_pos_artifact += 1
                note = "Unreliable estimate"
            elif classification == 'pass':
                status = "✓ PASS"
                n_pos_pass += 1
                n_pos_clean += 1
                note = ""
            else:
                status = "✗ FAIL"
                n_pos_clean += 1
                if not pd.isna(p) and p >= 0.05:
                    note = "Not significant"
                elif not pd.isna(k) and k >= 1:
                    note = "K≥1 (not relaxed)"
                else:
                    note = ""

            lines.append(f"{gene:<12} {k_str:<12} {p_str:<14} "
                         f"{spp:<6} {status:<22} {note}")

        lines.append("")

        # --- Negative controls ---
        lines.append("NEGATIVE CONTROLS (expect K ≈ 1, p > 0.05):")
        lines.append("-" * 78)
        lines.append(f"{'Gene':<12} {'Source':<14} {'K':<12} {'p-value':<14} "
                     f"{'Spp':<6} {'Status':<22} {'Note'}")
        lines.append("-" * 78)

        n_neg_pass = 0
        n_neg_clean = 0
        n_neg_error = 0
        n_neg_artifact = 0

        # Source labels for reporting
        neg_sources = {
            "ATP5F1A": "de Jonge 2007",
            "RPL13": "de Jonge 2007",
            "RPS13": "de Jonge 2007",
            "OPN1LW": "Zhao+ 2009",
            "VCP": "E&L 2013",
            "GPI": "E&L 2013",
            "REEP5": "E&L 2013",
            "VPS29": "E&L 2013",
            "SNRPD3": "E&L 2013",
        }

        for gene in GENE_ORDER_NEG:
            row_data = mode_data[mode_data['gene'] == gene]
            source = neg_sources.get(gene, "—")

            if len(row_data) == 0:
                lines.append(f"{gene:<12} {source:<14} {'—':<12} {'—':<14} "
                             f"{'—':<6} {'NOT FOUND':<22}")
                continue

            row = row_data.iloc[0]
            k = row['K']
            p = row['p_value']

            spp_col = [c for c in ['n_species', 'species_count'] if c in row.index]
            spp = ''
            if spp_col and not pd.isna(row[spp_col[0]]):
                spp = str(int(row[spp_col[0]]))

            classification = classify_result(row)
            k_str = f"{k:.4f}" if not pd.isna(k) else "N/A"
            p_str = f"{p:.2e}" if not pd.isna(p) else "N/A"

            if classification == 'error':
                status = "? ERROR"
                n_neg_error += 1
                note = "RELAX failed"
            elif classification == 'k_zero_artifact':
                status = "⚠ K≈0 ARTIFACT"
                n_neg_artifact += 1
                note = "Numerical boundary"
            elif classification == 'k_high_artifact':
                status = f"⚠ K>{K_UNRELIABLE_THRESHOLD:.0f} ARTIFACT"
                n_neg_artifact += 1
                note = "Unreliable estimate"
            elif classification == 'pass':
                status = "✓ PASS"
                n_neg_pass += 1
                n_neg_clean += 1
                note = ""
            else:
                status = "✗ FAIL"
                n_neg_clean += 1
                if not pd.isna(k) and k < 0.5:
                    note = "Unexpected relaxation"
                else:
                    note = "K outside range"

            lines.append(f"{gene:<12} {source:<14} {k_str:<12} {p_str:<14} "
                         f"{spp:<6} {status:<22} {note}")

        lines.append("")

        # Summary
        total_pass = n_pos_pass + n_neg_pass
        total_clean = n_pos_clean + n_neg_clean
        total_error = n_pos_error + n_neg_error
        total_artifact = n_pos_artifact + n_neg_artifact
        pct = 100 * total_pass / total_clean if total_clean > 0 else 0

        lines.append("SUMMARY:")
        lines.append(f"  Positive controls: {n_pos_pass}/{n_pos_clean} passed "
                     f"(excl. {n_pos_error} errors, {n_pos_artifact} artifacts)")
        lines.append(f"  Negative controls: {n_neg_pass}/{n_neg_clean} passed "
                     f"(excl. {n_neg_error} errors, {n_neg_artifact} artifacts)")
        lines.append(f"  Overall: {total_pass}/{total_clean} ({pct:.0f}%) — "
                     f"excluding {total_error} errors and {total_artifact} artifacts")
        lines.append("")

    # Cross-mode comparison
    lines.append("=" * 80)
    lines.append("CROSS-MODE COMPARISON")
    lines.append("=" * 80)
    lines.append("")
    lines.append(f"{'Gene':<12} {'Type':<10} {'Gen K':<10} {'Gen p':<12} "
                 f"{'Gen':<8} {'Spc K':<10} {'Spc p':<12} {'Spc':<8}")
    lines.append("-" * 82)

    for gene in GENE_ORDER_POS + GENE_ORDER_NEG:
        ctype = 'positive' if gene in ALL_POSITIVE_GENES else 'negative'
        gen_row = df[(df['mode'] == 'general') & (df['gene'] == gene)]
        spc_row = df[(df['mode'] == 'specific') & (df['gene'] == gene)]

        def fmt_row(row_data):
            if len(row_data) == 0:
                return '—', '—', 'MISS'
            r = row_data.iloc[0]
            k_s = f"{r['K']:.4f}" if not pd.isna(r['K']) else "N/A"
            p_s = f"{r['p_value']:.2e}" if not pd.isna(r['p_value']) else "N/A"
            cl = classify_result(r)
            if cl == 'error':
                st = 'ERROR'
            elif cl == 'pass':
                st = 'PASS'
            elif cl in ('k_zero_artifact', 'k_high_artifact'):
                st = 'ARTFCT'
            else:
                st = 'FAIL'
            return k_s, p_s, st

        gk, gp, gs = fmt_row(gen_row)
        sk, sp, ss = fmt_row(spc_row)

        # Highlight discordant results
        discord = ""
        if gs != ss and gs not in ('MISS', 'ERROR', 'ARTFCT') and ss not in ('MISS', 'ERROR', 'ARTFCT'):
            discord = " ← DISCORDANT"

        lines.append(f"{gene:<12} {ctype:<10} {gk:<10} {gp:<12} {gs:<8} "
                     f"{sk:<10} {sp:<12} {ss:<8}{discord}")

    lines.append("")

    # Write report
    report_text = "\n".join(lines)
    report_file = os.path.join(output_dir, "validation_report.txt")
    with open(report_file, 'w') as f:
        f.write(report_text)

    print(report_text)
    print(f"\nReport saved: {report_file}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    if len(sys.argv) > 1:
        results_file = sys.argv[1]
    else:
        results_file = DEFAULT_RESULTS

    if not os.path.exists(results_file):
        print(f"ERROR: Results file not found: {results_file}")
        print(f"Expected: {results_file}")
        print("Usage: python3 visualize_control_validation_v4.py [results.tsv]")
        sys.exit(1)

    print("=" * 60)
    print("CONTROL VALIDATION VISUALIZATION v4 (for v21 gene set)")
    print("=" * 60)
    print(f"\nInput: {results_file}")

    output_dir = "results/11_figures/validation"
    os.makedirs(output_dir, exist_ok=True)

    df = load_results(results_file)
    print(f"Loaded {len(df)} results across {df['mode'].nunique()} modes")
    print(f"Genes: {sorted(df['gene'].unique())}")

    # Create all figures
    print("\n--- Creating visualizations ---")

    print("\n[1] Paired barplot (general vs specific)...")
    create_paired_barplot(df, output_dir)

    print("\n[2] Forest plot...")
    create_forest_plot(df, output_dir)

    print("\n[3] Per-mode barplots...")
    for mode in ['general', 'specific']:
        if mode in df['mode'].values:
            create_barplot(df, mode, output_dir)

    print("\n[4] Generating report...")
    generate_report(df, output_dir)

    print("\n" + "=" * 60)
    print("All outputs saved to:", output_dir)
    print("=" * 60)


if __name__ == "__main__":
    main()
