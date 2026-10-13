# ===========================================================================
# chr15 haplotype segment (segA-segG) k-mer coverage analysis:
# dimensionality reduction (PCA/UMAP/t-SNE), correlation structure between
# segments, manual genotype calling (h1h1/h1h2/h2h2) based on segment
# coverage ratios, and statistical comparison of segment "fold expansion"
# between genotypes.
#
# ===========================================================================

# Get chr15 A-F coverage plots
setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure2/k-mer_analysis/binned_synteny/")

# ---------------------------------------------------------------------
# Load per-sample k-mer counts for each chr15 synteny segment (segA-segG),
# using sample ID as rownames and dropping the now-redundant sample column
# ---------------------------------------------------------------------
seg_cov <- read.table("mHetGla_chr15_binned_synteny_kmer_counts.txt", header = TRUE, as.is = TRUE, sep = " ")
rownames(seg_cov) <- seg_cov$sample
seg_cov <- seg_cov[, colnames(seg_cov) != "sample"]
origin <- seg_cov$origin  # (extracted but not directly used below - origin is re-pulled from seg_cov$origin later)

library(uwot)
library(Matrix)
library(irlba)
library(ggplot2)
library(ggrepel)

# ---------------------------------------------------------------------
# Restrict to the 7 chr15 segments (A-G) and normalize each segment's
# k-mer depth by segG's depth, giving a relative coverage ratio per
# segment per sample
# ---------------------------------------------------------------------
seg_cov_kmers <- seg_cov[, paste0("seg", LETTERS[1:7])]
seg_cov_kmers_depth <- seg_cov_kmers
for (i in names(seg_cov_kmers_depth)) {
  seg_cov_kmers_depth[, i] <- seg_cov_kmers_depth[, i] / seg_cov_kmers_depth[, "segG"]
}

# ---------------------------------------------------------------------
# PCA on normalized segment coverage ratios (all 7 segments, including
# segG which is now a constant column of 1s post-normalization)
# ---------------------------------------------------------------------
seg_cov_pca <- prcomp(seg_cov_kmers_depth[, paste0("seg", LETTERS[1:7])])
pca_cov <- seg_cov_pca$x[, 1:2]
loadings <- seg_cov_pca$rotation
pca_cov <- as.data.frame(pca_cov)

# UMAP embedding of samples using the first 5 principal components
umap_res <- umap(seg_cov_pca$x[, 1:5])

# Reorder seg_cov to match the sample order returned by umap() (rownames
# of umap_res are inherited from seg_cov_pca$x, so this should already
# match, but this re-indexing enforces it explicitly)
seg_cov <- seg_cov[rownames(umap_res), ]

library(Rtsne)
library(pheatmap)

# t-SNE embedding
tsne_out <- Rtsne(seg_cov_kmers_depth, perplexity = 6, scale = FALSE)

seg_cov_kmers_depth1 <- seg_cov_kmers_depth
# seg_cov_kmers_depth$segE <- log1p(seg_cov_kmers_depth$segE)

# Heatmap of normalized segment coverage per sample, clustered by row (sample)
pheatmap(seg_cov_kmers_depth, clustering_method = "ward.D2", cluster_cols = FALSE)

# ---------------------------------------------------------------------
# Correlation structure between segments A-F (segG excluded here since
# it's the normalization reference and would trivially correlate/be
# constant), visualized as a circular correlation plot
# ---------------------------------------------------------------------
cor_cov <- cor(seg_cov_kmers_depth[, paste0("seg", LETTERS[1:6])])
library(corrplot)
pdf(file = "mHetGla_chr15_segment_corrplot.pdf")
print(corrplot(cor_cov, "circle", type = "upper", tl.cex = 1.5,
               addCoef.col = "white"))
dev.off()

# ---------------------------------------------------------------------
# Build a combined per-sample data frame: UMAP coordinates, key segment
# coverage ratios (raw and log-transformed segE), PCA coordinates,
# sample label, and origin metadata

# ---------------------------------------------------------------------
df <- data.frame(UMAP1 = umap_res[, 1],
                 UMAP2 = umap_res[, 2],
                 segB = seg_cov_kmers_depth$segB,
                 segC = seg_cov_kmers_depth$segC,
                 segE = seg_cov_kmers_depth$segE,
                 segD = seg_cov_kmers_depth$segD,
                 logsegE = log1p(seg_cov_kmers_depth$segE),
                 PC1 = pca_cov$PC1,
                 PC2 = pca_cov$PC2,
                 label = rownames(umap_res),
                 origin = seg_cov$origin)

library(ggplot2)

# UMAP scatter colored by sample origin
g1 <- ggplot(df, aes(UMAP1, UMAP2, color = origin, label = label)) +
  geom_point(size = 5) + theme_classic()
ggsave("mHetGla_chr15_synteny_kmer_umap_v0619.png", g1, width = 15, height = 10)

# Scatter of log(segE) vs segB, colored by origin - exploring whether
# these two segments separate samples by genotype/origin
g2 <- ggplot(df, aes(logsegE, segB, color = origin, label = label)) +
  geom_point(size = 3, alpha = 0.5) + theme_classic()
ggsave("mHetGla_chr15_synteny_kmer_logsegE_vs_segB_v0619.png", g2, width = 6, height = 6)
ggsave("mHetGla_chr15_synteny_kmer_logsegE_vs_segB_v0619.pdf", g2, width = 6, height = 4)

# Additional exploratory segment-pair scatterplots
g3 <- ggplot(df, aes(segD, segB, color = origin, label = label)) +
  geom_point(size = 3, alpha = 0.5) + theme_classic()
g4 <- ggplot(df, aes(logsegE, segC, color = origin, label = label)) +
  geom_point(size = 3, alpha = 0.5) + theme_classic()
ggsave("mHetGla_chr15_synteny_kmer_segD_vs_segB_v0619.png", g3, width = 6, height = 6)
ggsave("mHetGla_chr15_synteny_kmer_segD_vs_segB_v0619.pdf", g3, width = 6, height = 6)
ggsave("mHetGla_chr15_synteny_kmer_logsegE_vs_segC_v0619.png", g4, width = 6, height = 6)
ggsave("mHetGla_chr15_synteny_kmer_logsegE_vs_segC_v0619.pdf", g4, width = 6, height = 6)

# ggsave("~/Desktop/NMR_multi_assembly/v090226/chr15_sequences/meryl/plot.hap2.png", p, width=15, height=10)

# ---------------------------------------------------------------------
# PCA scatterplots (PC2 vs PC1), colored variously by origin, logsegE,
# and segB, then combined into a single multi-panel figure via patchwork
# ---------------------------------------------------------------------
pca1 <- ggplot(df, aes(PC2, PC1, color = origin, label = label)) +
  geom_point(size = 3) + theme_bw() +
  ggtitle("sample origin")
pca2 <- ggplot(df, aes(PC2, PC1, color = logsegE, label = label)) +
  geom_point(size = 3) + theme_bw() +
  ggtitle("logsegE")
pca3 <- ggplot(df, aes(PC2, PC1, color = segB, label = label)) +
  geom_point(size = 3) + theme_bw() +
  ggtitle("segB")
pca1 + pca2 + pca3  # requires patchwork to already be loaded/attached for `+` composition to work here

# Repeat of the logsegE vs segB scatter, now with theme_bw() and
# a title, combined side-by-side with the UMAP plot (g1) above
g2 <- ggplot(df, aes(logsegE, segB, color = origin, label = label)) +
  geom_point(size = 3) + theme_bw() + geom_text_repel() + ggtitle("Normalized")
g1 + g2

# ---------------------------------------------------------------------
# Load an alternative/unnormalized coverage table for a subset of
# Kenyan samples, restrict df to just those samples, and overwrite
# their origin labels with values from this new table
# ---------------------------------------------------------------------
unnorm <- read.table("unnormalized_samples_kenya.txt", header = T, as.is = T)
df <- df[unnorm$sample, ]
df$origin <- unnorm$origin

g1 <- ggplot(df, aes(segE, segB, color = origin, label = label)) +
  geom_point(size = 3) + theme_bw() + geom_text_repel() + ggtitle("Normalized")
g2 <- ggplot(unnorm, aes(segE, segB, color = origin, label = sample)) +
  geom_point(size = 3) + theme_bw() + ggtitle("Unnormalized")

library(patchwork)


library(ggplot2)
g1 + g2  # duplicate of the combined plot above

library(ggrepel)
g1 <- ggplot(df, aes(UMAP1, UMAP2, color = segD, label = label)) +
  geom_point(size = 3) + theme_bw() +
  scale_color_viridis_c() +
  geom_text_repel()

unnorm$origin <- factor(unnorm$origin)
g1 <- ggplot(unnorm, aes(segE, segB, color = origin, label = sample)) +
  geom_point(size = 3) + theme_bw() +
  geom_text_repel()

# ---------------------------------------------------------------------
# UMAP colored by each segment's coverage ratio individually, to look
# for spatial structure in the embedding associated with specific
# segments (segD appears twice - g0 and g3 are identical plots)
# ---------------------------------------------------------------------
g0 <- ggplot(df, aes(UMAP1, UMAP2, color = segD)) +
  geom_point(size = 3) + theme_bw() +
  scale_color_viridis_c()
g2 <- ggplot(df, aes(UMAP1, UMAP2, color = segB)) +
  geom_point(size = 3) + theme_bw() +
  scale_color_viridis_c()
g3 <- ggplot(df, aes(UMAP1, UMAP2, color = segD)) +
  geom_point(size = 3) + theme_bw() +
  scale_color_viridis_c()
g4 <- ggplot(df, aes(UMAP1, UMAP2, color = segC)) +
  geom_point(size = 3) + theme_bw() +
  scale_color_viridis_c()
g0 + g1 + g2 + g3 + g4

# Additional correlation visualizations (circle plot, heatmap-style with
# numbers, and a numbered upper-triangle version)
library(corrplot)
corrplot(cor_cov)
g0 <- pheatmap(cor_cov, display_numbers = TRUE, fontsize = 14)
corrplot(cor_cov, type = "upper", tl.cex = 1.5,
         addCoef.col = "black")

library(patchwork)
g1 + g2

g2 <- ggplot(df, aes(UMAP1, UMAP2, color = segD)) +
  geom_point(size = 3) + theme_bw() +
  scale_color_viridis_c()

###
g2 <- ggplot(df, aes(logsegE, segB, color = origin, label = label)) +
  geom_point(size = 3, alpha = 0.5) + theme_classic()
#seg_cov$logsegE <- log1p(seg_cov$segE)


# --- 1. Your data ---

# ---------------------------------------------------------------------
# Rebuild df for the full sample set
# ---------------------------------------------------------------------
df <- data.frame(UMAP1 = umap_res[, 1],
                 UMAP2 = umap_res[, 2],
                 segB = seg_cov_kmers_depth$segB,
                 segC = seg_cov_kmers_depth$segC,
                 segE = seg_cov_kmers_depth$segE,
                 segD = seg_cov_kmers_depth$segD,
                 logsegE = log1p(seg_cov_kmers_depth$segE),
                 PC1 <- pca_cov$PC1,
                 PC2 <- pca_cov$PC2,
                 label = rownames(umap_res),
                 origin = seg_cov$origin)

# ---------------------------------------------------------------------
# Manually call a 3-genotype classification per sample based on
# logsegE and segB coverage ratio thresholds:
# - logsegE < 1.5              -> "h1h1" (homozygous haplotype 1)
# - else if segB < 0.4         -> "h2h2" (homozygous haplotype 2)
# - otherwise                  -> "h1h2" (heterozygous)
# These thresholds appear to have been chosen by eye from the scatter
# plots generated above (e.g. logsegE vs segB), not computed statistically.
# ---------------------------------------------------------------------
df$genotype <- ""
for (i in 1:nrow(df)) {
  if (df[i, "logsegE"] < 1.5) {
    df$genotype[i] <- "h1h1"
  } else if (df[i, "segB"] < 0.4) {
    df$genotype[i] <- "h2h2"
  } else {
    df$genotype[i] <- "h1h2"
  }
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# --- 1. Your data ---
# Assumes df has columns: genotype, origin, segC, segE (numeric)

# --- 2. Segments to analyze, h1h1 as reference for both ---
seg_cols <- c("segC", "segE")

# ---------------------------------------------------------------------
# --- 3. Keep all three genotypes, reshape to long format ---
# Restrict to samples with a called genotype, order genotype levels
# logically (h1h1 -> h1h2 -> h2h2), and pivot segC/segE into long format
# for faceted plotting/comparison
# ---------------------------------------------------------------------
df_long <- df %>%
  filter(genotype %in% c("h1h1", "h1h2", "h2h2")) %>%
  mutate(genotype = factor(genotype, levels = c("h1h1", "h1h2", "h2h2"))) %>%
  select(genotype, origin, all_of(seg_cols)) %>%
  pivot_longer(cols = all_of(seg_cols), names_to = "segment", values_to = "value")

# ---------------------------------------------------------------------
# --- 4. Fold expansion relative to h1h1 mean (same ref for both segments) ---
# For each segment, compute the mean coverage ratio among h1h1 samples
# as the reference baseline, then express every sample's value as a
# fold-change relative to that baseline
# ---------------------------------------------------------------------
df_fold <- df_long %>%
  group_by(segment) %>%
  mutate(
    ref_mean = mean(value[genotype == "h1h1"], na.rm = TRUE),
    fold_expansion = value / ref_mean
  ) %>%
  ungroup() %>%
  mutate(facet_label = paste0(segment, " (ref = h1h1)"))

# --- 5. Summary stats ---
# Mean, SD, N, and SEM of fold expansion per segment x genotype
summary_stats <- df_fold %>%
  group_by(segment, genotype) %>%
  summarise(
    mean_fold = mean(fold_expansion, na.rm = TRUE),
    sd_fold   = sd(fold_expansion, na.rm = TRUE),
    n         = sum(!is.na(fold_expansion)),
    sem_fold  = sd_fold / sqrt(n),
    .groups = "drop"
  )
print(summary_stats)

# ---------------------------------------------------------------------
# --- 6. t-tests: h1h1 vs h2h2 only, per segment ---
# For each segment, run a two-sample t-test comparing fold expansion
# between h1h1 and h2h2 genotypes (h1h2 excluded from this comparison),
# report the mean difference and 95% CI (oriented as h2h2 - h1h1),
# then BH-adjust p-values across segments and assign significance stars
# ---------------------------------------------------------------------
t_test_h1h1_h2h2 <- df_fold %>%
  filter(genotype %in% c("h1h1", "h2h2")) %>%
  mutate(genotype = factor(genotype, levels = c("h1h1", "h2h2"))) %>%
  group_by(segment) %>%
  group_modify(~ {
    tt <- t.test(fold_expansion ~ relevel(.x$genotype, ref = "h1h1"), data = .x)
    tibble(
      mean_diff = diff(tt$estimate) * -1,
      ci_low    = -tt$conf.int[2],
      ci_high   = -tt$conf.int[1],
      t_stat    = tt$statistic,
      df        = tt$parameter,
      p_value   = tt$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    fold_ci_low  = ci_low + 1,
    fold_ci_high = ci_high + 1,
    p_adj        = p.adjust(p_value, method = "BH"),
    sig          = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )
print(t_test_h1h1_h2h2)

# --- 7. Map origin -> color ---
# Same manual origin color scheme used throughout the earlier k-mer
# depth heatmap scripts, for visual consistency across figures
origin_colors <- setNames(
  scales::hue_pal()(length(unique(df_fold$origin))),
  unique(df_fold$origin)
)
origin_colors["kenya-mixed"] <- "#C77CFF"
origin_colors["kenya-lab"]   <- "#00BFC4"
origin_colors["ethiopia-east"] <- "#F8766D"
origin_colors["ethiopia-south"] <- "#7CAE00"



# ---------------------------------------------------------------------
# --- 8. Boxplot ---
# Faceted boxplot of fold expansion (segC and segE) by genotype, with
# individual sample points jittered and colored by origin, plus mean
# +/- SE overlay. Two versions saved: linear scale and log10 scale.
# ---------------------------------------------------------------------
p <- ggplot(df_fold, aes(x = genotype, y = fold_expansion)) +
  geom_boxplot(aes(fill = genotype), width = 0.5, outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(color = origin), width = 0.15, size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5, fill = "white", color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.6) +
  facet_wrap(~ facet_label, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("h1h1" = "white", "h1h2" = "white", "h2h2" = "white")) +
  scale_color_manual(values = origin_colors) +
  labs(
    x = "Genotype",
    y = "Fold expansion (relative to h1h1)",
    title = "Fold expansion: segC and segE",
    fill = "Genotype",
    color = "Origin"
  ) +
  theme_classic(base_size = 8)
print(p)
ggsave("segC_segE_fold_expansion.pdf", p, width = 7, height = 3, dpi = 300, units = "in")

p_log <- p + scale_y_log10()
print(p_log)
ggsave("segC_segE_fold_expansion_log.png", p_log, width = 7, height = 5, dpi = 300, units = "in")
