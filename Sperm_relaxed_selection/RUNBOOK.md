# Reproduction runbook

All commands are run from the project root on a single workstation with 44
threads. HyPhy must be v2.5.84: RELAX estimates are not identical across
versions, and every number in the report was produced with that build.

    hyphy --version          # must report 2.5.84

---

## 1. Data preparation and orthology

    python3 part0_assess_data_quality_v3.py 44
    python3 part1_isoform_selection_v10.py
    python3 part2a_rename_samples_v2.py
    python3 part2b_clean_proteins_v2.py
    bash    part2c_run_orthofinder.sh
    python3 part4_og_symbol_map.py

## 2. Gene lists

    python3 part5a_validation_gene_list.py
    Rscript part5b_extract_go_gene_lists_v4.R
    python3 part5c_ddr_gene_list.py --xlsx Balmus_DDR_genes_Apr2026.xlsx

`part5b` writes the union of the three sperm GO terms plus the 17 targets. The
symbol count depends on the installed `org.Hs.eg.db`; record its version, as the
set has grown substantially between releases:

    Rscript -e 'cat(as.character(packageVersion("org.Hs.eg.db")))'

## 3. Orthogroup selection

    python3 part5d_select_orthogroups.py --gene_list gene_lists/validation.tsv \
        --outdir results/05_sets/validation --min_species 4 --min_single_copy_pct 0

    python3 part5d_select_orthogroups.py --gene_list results/11_ddr_atlas/ddr_gene_list.tsv \
        --outdir results/05_sets/ddr --min_species 15 --min_single_copy_pct 40

    python3 part5d_select_orthogroups.py --gene_list results/10_extended_sperm_genes/go_sperm_genes.tsv \
        --outdir results/05_sets/sperm --min_species 15 --min_single_copy_pct 40 \
        --rescue_list gene_lists/targets17.tsv

`--rescue_list` is required for the sperm set. The 17 targets are curated
candidates and are exempt from the admission filters; without it ZAN (27%
single-copy orthogroup) and SPEM3 (11 species) are dropped and the set is 721
genes rather than 723. The console prints `rescued past filters (targets):
SPEM3, ZAN`.

The validation set deliberately uses permissive filters (4 species, no
single-copy requirement) so that the control panel stresses the pipeline rather
than being pre-cleaned.

## 4. CDS assembly, alignment and trimming

    python3 part6_assemble_cds.py --set validation
    python3 part6_assemble_cds.py --set sperm
    python3 part6_assemble_cds.py --set ddr

    python3 part7_macse.py   --set validation --threads 44 --timeout_hours 3
    python3 part8_clipkit.py --set validation --threads 44

    python3 part7_macse.py   --set sperm --threads 44 --timeout_hours 3
    python3 part7_macse.py   --set ddr   --threads 44 --timeout_hours 3
    python3 part8_clipkit.py --set sperm --threads 44
    python3 part8_clipkit.py --set ddr   --threads 44

A small number of very long genes exceed the default 4 GB MACSE heap and fail.
Rerun those alone with a large heap, then rerun ClipKIT to absorb them:

    nohup python3 part7_macse.py --set sperm --threads 1 --timeout_hours 24 --xmx 80g \
        --only APOB,DNAH1,DNAH11,DNAH5,FSIP2,SYNE1,ZAN > logs/sperm_giants.log 2>&1 &
    nohup python3 part7_macse.py --set ddr --threads 1 --timeout_hours 24 --xmx 80g \
        --only HUWE1,KMT2C > logs/ddr_giants.log 2>&1 &

    python3 part8_clipkit.py --set sperm --threads 44
    python3 part8_clipkit.py --set ddr   --threads 44

Confirm all 17 targets are present before continuing. `targets17` has no
alignment stage of its own; it reads the sperm set's alignments, so this check
must pass first:

    for g in $(tail -n+2 gene_lists/targets17.tsv); do
      test -s results/08_clipkit/sperm/clipkit_trimmed/${g}_NT.clipkit.fa \
        && echo "$g: IN" || echo "$g: MISSING"; done

## 5. Selection analysis (RELAX)

    python3 part9_relax.py --set validation --threads 44
    python3 part9_relax.py --set targets17  --threads 44
    python3 part9_relax.py --set sperm      --threads 44 --relax_timeout 18000
    python3 part9_relax.py --set ddr        --threads 44 --relax_timeout 18000

Defaults are the frozen specification: `--models All`, one starting point,
outlier removal on, selection exemption restricted to the two curated focal
species. Do not pass `--models`, `--starting_points` or `--exempt` for a
reproduction run; they exist for the diagnostic analyses described below.

The raised timeout for the two atlas sets prevents a handful of very long genes
being recorded as timeouts. If any gene does time out, rerunning the same
command resumes: completed genes are skipped and only the unfinished ones are
recomputed.

**Collate each set once.** In batch mode `part9` both computes and writes the
per-scheme results tables. Genes that fail leave an error file but no JSON, so
every subsequent pass re-attempts them at full cost. Once a set has completed,
treat its tables as final and do not rerun it.

## 6. Classification and multiple-testing correction

    python3 part10_summarise.py
    python3 part11_fdr.py --set sperm --shared
    python3 part11_fdr.py --set targets17
    python3 part11_fdr.py --set ddr
    python3 part11_fdr.py --set validation

`part10` applies the boundary rules (K <= 0.01, K >= 5) and writes
`results/10_summary/{set}_summary.tsv`. `part11_fdr` adds Benjamini-Hochberg
q-values computed within each scheme over the testable genes, and `--shared`
additionally restricts the analysis to the genes testable in all three
single-species schemes and runs the paired McNemar tests.

**Preserve the validation summary.** The reported control-panel values come from
one designated draw. Any rerun of `part9 --set validation` followed by `part10`
will overwrite it, so keep a copy:

    cp results/10_summary/validation_summary.tsv \
       results/10_summary/validation_summary_DESIGNATED.tsv

## 7. Replicate stability of the control panel

RELAX is a single-fit estimator and its optimum is not fully reproducible for
genes with a flat likelihood surface. The control panel is therefore reported as
one designated draw plus replicates. To regenerate the stability analysis:

    mkdir -p results/09_relax/validation_replicates
    for i in 1 2 3 4; do
      rm -rf results/09_relax/validation
      python3 part9_relax.py --set validation --threads 44
      cp results/09_relax/validation/general/relax_results.tsv \
         results/09_relax/validation_replicates/val_run$i.tsv
    done

    # restore the designated draw as the on-disk result
    cut -f1-6 results/10_summary/validation_summary_DESIGNATED.tsv \
        > results/09_relax/validation/general/relax_results.tsv

    python3 part11b_validation_stability.py \
        results/09_relax/validation_replicates/val_run*.tsv \
        --official results/10_summary/validation_summary_DESIGNATED.tsv \
        --threshold both

Replicate draws will not reproduce exactly; that variability is the finding, and
the reported consensus is over draws rather than from any single one.

## 8. Figures

    # gene presence/absence across the 23 species (Fig 1)
    python3 05_gene_presence_absence.py

    # per-scheme K-value dot plots for the 17 targets (Figs 2-6)
    python3 visualize_17gene_relax_v7.py

    # K-value heatmap across schemes (Fig 7)
    python3 06_relax_heatmap_bathyergids.py --fdr

    # protein disruption panels for the 17 targets (Figs 8-14)
    python3 01_extract_focal_v5.py \
        --selected_from results/09_relax/sperm/allmolerat_vs_nonmolerat/pruned
    Rscript 02_plot_disruptions_v5.R

    # control-gene validation (Fig 15)
    python3 visualize_control_validation_v4.py
    python3 visualize_validation_v5_extended.py --view general

    # control-gene disruption panels
    python3 01_extract_controls_v1.py \
        --input_dir results/07_macse/validation/macse_aligned \
        --selected_from results/09_relax/validation/general/pruned
    Rscript 02_plot_disruptions_v5.R \
        --input_dir results/07_macse/validation/protein_disruption_figures \
        --output_dir results/07_macse/validation/protein_disruption_figures/figures

    # sociality gradient on the shared gene set (Fig 18)
    python3 fig_gradient_intersection_5scheme.py \
        --merged results/10_summary/sperm_summary.tsv \
        --outdir results/11_figures/atlas_sperm

    # atlas-wide volcano, gradient and overlap panels
    python3 visualize_relax_atlas_v6.py --set sperm --fdr
    python3 visualize_relax_atlas_v6.py --set ddr --fdr

`--selected_from` makes each reference row in the disruption panels show the
ortholog that RELAX actually analysed, rather than whichever copy appears first
in the alignment file.

## 9. Notes on reproducibility

Gene symbols follow current HGNC usage on disk. Two of the 17 targets were
renamed after the study was designed: GARIN3 (formerly FAM71B) and GARIN1A
(formerly FAM71F2). Files carry the current symbols; figures and tables print
the legacy ones.

`targets17` and `sperm` share alignments, so a gene's K value is identical in
both; the sets differ only in which schemes are applied.

Boundary-value genes (K <= 0.01 or K >= 5) are excluded from significance
testing and from the FDR denominator throughout, because a maximum-likelihood
estimate on a parameter boundary invalidates the chi-square approximation of the
likelihood-ratio statistic.
