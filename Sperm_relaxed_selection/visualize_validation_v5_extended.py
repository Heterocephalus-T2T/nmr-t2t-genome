#!/usr/bin/env python3
"""
visualize_validation_v5_extended.py

Control-gene validation figures for the v4/v22 run. Preserves the v3 layout
(paired barplot, forest plot, per-mode barplots, text report) via the base
module, and ADDS:

  (A) Three mode views for every figure: general, specific, and both
      (choose with --view general|specific|both; default builds all three).
  (B) Manual-grounded K~0 rendering. Per the HyPhy RELAX docs (omega_T =
      omega_R ^ K; K is an exponent with permitted range [0,50]), K~0 with a
      significant p-value is MAXIMAL RELAXATION at the parameter boundary, NOT
      an artifact. We render it in the relaxed-colour family with a distinct
      "boundary" marker, labelled "K~0 (strongly relaxed)", instead of grey
      "artifact". Consistent with the atlas interpretation.
  (C) Three new figures:
        1. classification_vs_K  - K value + significance per gene, showing that
           the relaxed/not classification rests on the p-value (stable) rather
           than the exact K (unstable on flat surfaces).
        2. expected_vs_observed - scatter of expected vs observed K with
           significance overlaid; correct-quadrant genes validate the method.
        3. confusion_matrix     - 2x2 control-type x significance; the cleanest
           statement that the method flags relaxation in positives, not in
           negatives, using the stable significance test.

All figures: PDF + 600 dpi PNG + SVG.

Usage:
    python3 visualize_validation_v5_extended.py [results.tsv] [--view both] [--outdir DIR]

Default input: results/control_gene_validation_v22/validation_results.tsv
"""

import os, sys, argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np

# ------------------------------------------------------------------ config
DEFAULT_RESULTS = "results/09_relax/validation/general/relax_results.tsv"
K_ZERO_THRESHOLD = 0.0001       # at/below this, K is at the relaxation boundary
K_UNRELIABLE_THRESHOLD = 5.0    # above this, K is genuinely unreliable (intensification blow-up)

ALL_POSITIVE_GENES = {
    "CRYAA": {"desc": "Lens crystallin aA"}, "CRYAB": {"desc": "Lens crystallin aB"},
    "CRYBB1": {"desc": "Lens crystallin bB1"}, "RHO": {"desc": "Rhodopsin"},
    "GRK1": {"desc": "Rhodopsin kinase"}, "OPN1SW": {"desc": "Blue opsin"},
    "OR1A1": {"desc": "Olfactory receptor"}, "OR2W1": {"desc": "Olfactory receptor"},
}
ALL_NEGATIVE_GENES = {
    "ATP5F1A": {"desc": "ATP synthase F1a"}, "RPL13": {"desc": "Ribosomal L13"},
    "RPS13": {"desc": "Ribosomal S13"}, "OPN1LW": {"desc": "Red/green opsin"},
    "VCP": {"desc": "Valosin-containing protein"}, "GPI": {"desc": "G6P isomerase"},
    "REEP5": {"desc": "Receptor accessory 5"}, "VPS29": {"desc": "Vacuolar sorting 29"},
    "SNRPD3": {"desc": "snRNP D3"},
}
GENE_ORDER_POS = ["CRYAA", "CRYAB", "CRYBB1", "RHO", "GRK1", "OPN1SW", "OR1A1", "OR2W1"]
GENE_ORDER_NEG = ["ATP5F1A", "RPL13", "RPS13", "OPN1LW", "VCP", "GPI", "REEP5", "VPS29", "SNRPD3"]

COLORS = {
    'pos_pass': '#2ecc71', 'pos_fail': '#e74c3c',
    'neg_pass': '#3498db', 'neg_fail': '#e67e22',
    'strong_relax': '#16a085',   # teal: K~0 strongly-relaxed (NEW, relaxed family)
    'k_high': '#95a5a6',         # grey: genuine K>5 unreliable
    'error': '#d5d8dc', 'missing': '#f0f0f0',
    'positive': '#27ae60', 'negative': '#7f8c8d',
}

def save_all(fig, output_dir, name):
    """Save PDF + 600dpi PNG + SVG."""
    base = os.path.join(output_dir, name)
    fig.savefig(base + '.pdf', bbox_inches='tight')
    fig.savefig(base + '.png', dpi=600, bbox_inches='tight')
    fig.savefig(base + '.svg', bbox_inches='tight')
    print(f"  Saved: {base}.(pdf|png|svg)")
    plt.close(fig)

# ------------------------------------------------------------------ data + classify
def load_results(filepath):
    df = pd.read_csv(filepath, sep='\t')
    # v5 pipeline layout: Gene/K/p_value/LRT/Scheme/Set -> gene/mode
    if 'Gene' in df.columns and 'gene' not in df.columns:
        df = df.rename(columns={'Gene': 'gene', 'Scheme': 'mode'})
    if 'mode' not in df.columns:
        df['mode'] = 'general'
    if 'gene' not in df.columns and 'gene_og' in df.columns:
        df['gene'] = df['gene_og'].str.extract(r'^([A-Z0-9]+)_OG')
    if 'control_type' not in df.columns or df['control_type'].isna().all():
        pos, neg = set(ALL_POSITIVE_GENES), set(ALL_NEGATIVE_GENES)
        df['control_type'] = df['gene'].apply(
            lambda g: 'positive' if g in pos else ('negative' if g in neg else 'unknown'))
    return df

def classify_result(row):
    """
    pass / fail / strong_relax / k_high_artifact / error

    NOTE (manual-grounded): K~0 is NOT lumped with K>5 as a generic 'artifact'.
      - K~0  -> 'strong_relax': maximal relaxation at the parameter boundary.
                For POSITIVE controls with p<0.05 this is a PASS-like success
                (relaxation detected). Rendered distinctly but in relaxed family.
      - K>5  -> 'k_high_artifact': genuine unreliable/intensification blow-up.
    """
    k, p, ctype = row['K'], row['p_value'], row['control_type']
    if pd.isna(k) or pd.isna(p):
        return 'error'
    if k > K_UNRELIABLE_THRESHOLD:
        return 'k_high_artifact'
    if k < K_ZERO_THRESHOLD:
        # boundary: strongly relaxed. Success depends on significance + control type.
        return 'strong_relax'
    if ctype == 'positive':
        return 'pass' if (k < 1 and p < 0.05) else 'fail'
    else:
        return 'pass' if (p >= 0.05 or (0.5 <= k <= 2.0)) else 'fail'

def is_relaxed_significant(row):
    """The STABLE classification: relaxed (K<1 incl. boundary) AND significant."""
    k, p = row['K'], row['p_value']
    if pd.isna(k) or pd.isna(p):
        return False
    return (k < 1.0) and (p < 0.05)

def get_color(row, classification):
    if classification == 'strong_relax':
        # relaxed family; positive control success or negative control's unexpected relaxation
        return COLORS['strong_relax']
    if classification == 'k_high_artifact':
        return COLORS['k_high']
    if classification == 'error':
        return COLORS['error']
    if row['control_type'] == 'positive':
        return COLORS['pos_pass'] if classification == 'pass' else COLORS['pos_fail']
    return COLORS['neg_pass'] if classification == 'pass' else COLORS['neg_fail']

def sig_label(p):
    if pd.isna(p): return ''
    if p < 0.001: return '***'
    if p < 0.01: return '**'
    if p < 0.05: return '*'
    return ''

# ------------------------------------------------------------------ (C1) classification vs K
def fig_classification_vs_K(df, mode, output_dir, suffix):
    """Per gene: K value bar + significance + relaxed/not call. Shows the call
    rests on significance, not exact K."""
    data = df[df['mode'] == mode].copy() if mode != 'both' else df.copy()
    if mode == 'both':
        # average is meaningless; use general for the K bar, annotate both sig
        data = df[df['mode'] == 'general'].copy()
    if len(data) == 0:
        print(f"  (classification_vs_K) no data for {mode}"); return

    order = [g for g in GENE_ORDER_POS + GENE_ORDER_NEG if g in data['gene'].values]
    data = data.set_index('gene').loc[order].reset_index()

    fig, ax = plt.subplots(figsize=(max(13, len(order) * 0.85), 6.5))
    x = np.arange(len(order))

    kdisp = data['K'].clip(upper=K_UNRELIABLE_THRESHOLD).fillna(0)
    colors, edgecs = [], []
    for _, r in data.iterrows():
        relaxed_sig = is_relaxed_significant(r)
        if r['control_type'] == 'positive':
            # positive controls SHOULD be relaxed+sig
            colors.append(COLORS['pos_pass'] if relaxed_sig else COLORS['pos_fail'])
        else:
            # negative controls should NOT be relaxed+sig
            colors.append(COLORS['neg_fail'] if relaxed_sig else COLORS['neg_pass'])
        edgecs.append('black')

    bars = ax.bar(x, kdisp, color=colors, edgecolor=edgecs, linewidth=0.6)

    # mark K~0 boundary bars with a teal cap line + label
    for bar, (_, r) in zip(bars, data.iterrows()):
        k, p = r['K'], r['p_value']
        if pd.isna(k):
            ax.annotate('ERR', xy=(bar.get_x()+bar.get_width()/2, 0.02),
                        ha='center', va='bottom', fontsize=7, color='red', fontweight='bold')
            continue
        # significance stars + relaxed/not tick
        stars = sig_label(p)
        relaxed_sig = is_relaxed_significant(r)
        tick = '\u2713' if (relaxed_sig == (r['control_type'] == 'positive')) else '\u2717'
        if k < K_ZERO_THRESHOLD:
            # draw a small teal marker at the base to indicate boundary/strong relaxation
            ax.plot(bar.get_x()+bar.get_width()/2, 0.03, marker='v',
                    color=COLORS['strong_relax'], markersize=8, zorder=6)
            lbl = f"K~0{stars}"
        elif k > K_UNRELIABLE_THRESHOLD:
            lbl = f">{K_UNRELIABLE_THRESHOLD:.0f}{stars}"
        else:
            lbl = f"{k:.2f}{stars}"
        ax.annotate(lbl, xy=(bar.get_x()+bar.get_width()/2, min(k, K_UNRELIABLE_THRESHOLD)),
                    xytext=(0, 3), textcoords='offset points', ha='center', va='bottom',
                    fontsize=7, rotation=0)

    ax.axhline(y=1.0, color='black', linestyle='--', linewidth=1.5, alpha=0.7, label='K=1')
    ax.axhspan(0, 1, alpha=0.06, color='green')

    n_pos = len([g for g in order if g in ALL_POSITIVE_GENES])
    if 0 < n_pos < len(order):
        ax.axvline(n_pos - 0.5, color='grey', linestyle=':', linewidth=1)

    ax.set_xticks(x); ax.set_xticklabels(order, rotation=45, ha='right', fontsize=10, fontweight='bold')
    ax.set_ylabel('K (relaxation parameter)', fontsize=12)
    mode_lbl = 'General' if mode == 'general' else ('Specific' if mode == 'specific' else 'General (K) + both-mode sig')
    ax.set_title(f'Classification by significance test \u2014 {mode_lbl}\n'
                 '(bar colour = correct/incorrect call; call uses p<0.05 & K<1, not exact K)',
                 fontsize=12)
    ax.set_ylim(0, min(kdisp.max() * 1.3 if kdisp.max() > 0 else 2, 2.6))

    legend = [
        mpatches.Patch(color=COLORS['pos_pass'], label='Positive: correctly relaxed+sig'),
        mpatches.Patch(color=COLORS['pos_fail'], label='Positive: missed (not sig/relaxed)'),
        mpatches.Patch(color=COLORS['neg_pass'], label='Conserved: correctly not-relaxed'),
        mpatches.Patch(color=COLORS['neg_fail'], label='Conserved: unexpectedly relaxed'),
        mlines.Line2D([], [], marker='v', color=COLORS['strong_relax'], linestyle='None',
                      markersize=8, label='K~0 (strongly relaxed, boundary)'),
        mlines.Line2D([], [], color='black', linestyle='--', label='K=1 (neutral)'),
    ]
    ax.legend(handles=legend, loc='upper right', fontsize=7.5, framealpha=0.9)
    ax.text(0.01, 0.98, 'K~0 = maximal relaxation (exponent boundary),\n'
                        'not an artifact. Call rests on the stable\n'
                        'significance test, * p<0.05 ** p<0.01 *** p<0.001',
            transform=ax.transAxes, fontsize=7, va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.yaxis.grid(True, linestyle=':', alpha=0.4); ax.set_axisbelow(True)
    plt.tight_layout()
    save_all(fig, output_dir, f'validation_classification_vs_K_{suffix}')

# ------------------------------------------------------------------ (C2) expected vs observed
def fig_expected_vs_observed(df, mode, output_dir, suffix):
    """Scatter: expected K region vs observed K, significance as fill."""
    modes = ['general', 'specific'] if mode == 'both' else [mode]
    fig, ax = plt.subplots(figsize=(8.5, 7))

    # expected: positive -> target <1 (plot expected at 0.4); negative -> ~1
    exp_pos, exp_neg = 0.4, 1.0
    mode_marker = {'general': 'o', 'specific': 's'}

    for m in modes:
        d = df[df['mode'] == m]
        for _, r in d.iterrows():
            k, p = r['K'], r['p_value']
            if pd.isna(k): continue
            ctype = r['control_type']
            exp = exp_pos if ctype == 'positive' else exp_neg
            kd = min(max(k, 0.0), 3.0)
            sig = (not pd.isna(p)) and p < 0.05
            # colour by correctness
            if ctype == 'positive':
                correct = is_relaxed_significant(r)
                col = COLORS['pos_pass'] if correct else COLORS['pos_fail']
            else:
                correct = not is_relaxed_significant(r)
                col = COLORS['neg_pass'] if correct else COLORS['neg_fail']
            mfc = col if sig else 'white'
            ax.scatter(exp, kd, marker=mode_marker.get(m, 'o'), s=90,
                       facecolor=mfc, edgecolor=col, linewidth=1.6, zorder=5, alpha=0.9)
            ax.annotate(r['gene'], (exp, kd), fontsize=6, xytext=(4, 3),
                        textcoords='offset points', color=col)

    # guide lines
    ax.axhline(1.0, color='black', linestyle='--', linewidth=1.2, alpha=0.6)
    ax.axvline(0.7, color='grey', linestyle=':', linewidth=1, alpha=0.5)
    ax.axhspan(0, 1, alpha=0.05, color='green')
    ax.text(0.4, 2.7, 'Expected: relaxed (K<1)\npositive controls', ha='center', fontsize=8, color=COLORS['positive'])
    ax.text(1.0, 2.7, 'Expected: neutral (K~1)\nconserved controls', ha='center', fontsize=8, color=COLORS['negative'])

    ax.set_xlabel('Expected K (by control type)', fontsize=12)
    ax.set_ylabel('Observed K', fontsize=12)
    ax.set_xlim(0, 1.5); ax.set_ylim(-0.05, 3.0)
    view_lbl = {'general': 'General', 'specific': 'Specific', 'both': 'General + Specific'}[mode]
    ax.set_title(f'Expected vs observed K \u2014 {view_lbl}\n(filled = significant p<0.05; colour = correct call)', fontsize=12)

    legend = [
        mpatches.Patch(color=COLORS['pos_pass'], label='Positive correct'),
        mpatches.Patch(color=COLORS['pos_fail'], label='Positive incorrect'),
        mpatches.Patch(color=COLORS['neg_pass'], label='Conserved correct'),
        mpatches.Patch(color=COLORS['neg_fail'], label='Conserved incorrect'),
    ]
    if mode == 'both':
        legend += [mlines.Line2D([], [], marker='o', color='grey', linestyle='None', label='General'),
                   mlines.Line2D([], [], marker='s', color='grey', linestyle='None', label='Specific')]
    legend += [mlines.Line2D([], [], marker='o', color='grey', linestyle='None', markerfacecolor='grey', label='Significant (filled)'),
               mlines.Line2D([], [], marker='o', color='grey', linestyle='None', markerfacecolor='white', label='Not sig (open)')]
    ax.legend(handles=legend, loc='upper left', fontsize=7.5, framealpha=0.9)
    ax.grid(True, linestyle=':', alpha=0.35); ax.set_axisbelow(True)
    plt.tight_layout()
    save_all(fig, output_dir, f'validation_expected_vs_observed_{suffix}')

# ------------------------------------------------------------------ (C3) confusion matrix
def fig_confusion_matrix(df, mode, output_dir, suffix):
    """2x2: control type (positive/conserved) x call (relaxed+sig / not).
    Uses the stable significance-based classification."""
    modes = ['general', 'specific'] if mode == 'both' else [mode]
    n = len(modes)
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 4.8), squeeze=False)

    for idx, m in enumerate(modes):
        ax = axes[0][idx]
        d = df[df['mode'] == m]
        # counts
        tp = fp = tn = fn = 0
        cells = {'TP': [], 'FN': [], 'FP': [], 'TN': []}
        for _, r in d.iterrows():
            if pd.isna(r['K']) or pd.isna(r['p_value']):
                continue
            relaxed_sig = is_relaxed_significant(r)
            if r['control_type'] == 'positive':
                if relaxed_sig: tp += 1; cells['TP'].append(r['gene'])
                else: fn += 1; cells['FN'].append(r['gene'])
            else:
                if relaxed_sig: fp += 1; cells['FP'].append(r['gene'])
                else: tn += 1; cells['TN'].append(r['gene'])

        mat = np.array([[tp, fn], [fp, tn]])
        # colour: correct cells green-ish, error cells red-ish
        cell_colors = np.array([['#d5f5e3', '#fadbd8'], ['#fadbd8', '#d6eaf8']])
        ax.imshow(np.zeros_like(mat), cmap='Greys', vmin=0, vmax=1)
        for i in range(2):
            for j in range(2):
                ax.add_patch(mpatches.Rectangle((j-0.5, i-0.5), 1, 1,
                             facecolor=cell_colors[i][j], edgecolor='black', lw=1.2))
                key = [['TP', 'FN'], ['FP', 'TN']][i][j]
                genes = cells[key]
                gtxt = '\n'.join(genes) if genes else ''
                ax.text(j, i-0.12, str(mat[i][j]), ha='center', va='center',
                        fontsize=20, fontweight='bold')
                ax.text(j, i+0.28, gtxt, ha='center', va='center', fontsize=5.5, color='#555')
                ax.text(j, i-0.34, key, ha='center', va='center', fontsize=7, color='#777', style='italic')

        ax.set_xlim(-0.5, 1.5); ax.set_ylim(1.5, -0.5)
        ax.set_xticks([0, 1]); ax.set_xticklabels(['Called relaxed\n(K<1 & p<0.05)', 'Not called'], fontsize=8)
        ax.set_yticks([0, 1]); ax.set_yticklabels(['Positive\ncontrol', 'Conserved\ncontrol'], fontsize=8)
        sens = tp / (tp + fn) if (tp + fn) else 0
        spec = tn / (tn + fp) if (tn + fp) else 0
        ax.set_title(f'{m.capitalize()} mode\nsensitivity={sens:.0%}  specificity={spec:.0%}', fontsize=10)

    fig.suptitle('Validation confusion matrix (significance-based call)', fontsize=12, y=1.02)
    plt.tight_layout()
    save_all(fig, output_dir, f'validation_confusion_matrix_{suffix}')

# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('results', nargs='?', default=DEFAULT_RESULTS)
    ap.add_argument('--view', choices=['general', 'specific', 'both', 'all'], default='all',
                    help="which mode view(s) to build; 'all' builds general, specific, and both")
    ap.add_argument('--outdir', default=None)
    args = ap.parse_args()

    if not os.path.exists(args.results):
        sys.exit(f"ERROR: results file not found: {args.results}")

    outdir = args.outdir or 'results/11_figures/validation'
    figdir = os.path.join(outdir, 'figures_v5')
    os.makedirs(figdir, exist_ok=True)

    print("=" * 64)
    print("CONTROL VALIDATION FIGURES v5 (extended)")
    print("=" * 64)
    print(f"Input:  {args.results}")
    print(f"Output: {figdir}")

    df = load_results(args.results)
    print(f"Loaded {len(df)} rows; modes: {sorted(df['mode'].unique())}")

    views = ['general', 'specific', 'both'] if args.view == 'all' else [args.view]
    for v in views:
        suffix = v
        print(f"\n--- Building view: {v} ---")
        fig_classification_vs_K(df, v, figdir, suffix)
        fig_expected_vs_observed(df, v, figdir, suffix)
        fig_confusion_matrix(df, v, figdir, suffix)

    print("\n" + "=" * 64)
    print("NOTE: These are the NEW figures (classification-vs-K, expected-vs-observed,")
    print("confusion matrix), rendered for general/specific/both. For the v3-style")
    print("paired-barplot + forest + per-mode + text report, run the base script")
    print("(visualize_control_validation_v4_v2.py) on the same results file.")
    print("K~0 is rendered as strongly-relaxed (boundary), not artifact, per RELAX docs.")
    print("=" * 64)

if __name__ == '__main__':
    main()
