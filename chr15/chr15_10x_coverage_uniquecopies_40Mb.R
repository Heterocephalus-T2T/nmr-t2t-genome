library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(dendextend)
library(patchwork)
options(scipen = 999)

# ===========================================================================
# k-mer read-depth / uniqueness profiling across chr15 bins, for two
# haplotype assemblies (hap1, hap2), across a panel of NMR samples.
# Produces:
#   1. Hierarchical clustering dendrograms of samples by normalized k-mer
#      depth profile along chr15 (per haplotype)
#   2. Heatmap-style tile plots of normalized depth along chr15, samples
#      ordered by clustering, colored/labeled by sample origin
#   3. The same two outputs but using "unique k-mer fraction covered"
#      instead of normalized depth
# ===========================================================================

setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure2/k-mer_analysis/coverage_plots/")

# Full list of NMR sample IDs to be profiled
samples_full = c("nmr_1","nmr_10","nmr_11","nmr_12","nmr_13","nmr_14",
                 "nmr_1901_10","nmr_1901_11","nmr_1901_12","nmr_1901_13",
                 "nmr_1901_14","nmr_1901_15","nmr_1901_16","nmr_1901_19",
                 "nmr_1901_20","nmr_1901_21","nmr_1901_24","nmr_1901_25",
                 "nmr_1901_3","nmr_1901_5","nmr_1901_7","nmr_1901_9","nmr_2",
                 "nmr_3","nmr_4","nmr_5","nmr_6","nmr_7","nmr_8","nmr_9")

# Load per-sample origin/metadata key (used later to color-code samples
# by their geographic/colony origin, e.g. "kenya-lab", "kenya-mixed")
seg_cov <- read.table("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure2/k-mer_analysis/binned_synteny/mHetGla_chr15_binned_synteny_kmer_counts.txt",
                      header = TRUE, as.is = TRUE, sep = " ")
origin_key <- seg_cov

# ===========================================================================
# SECTION 1: Hap1 vs. Hap2 - NORMALIZED K-MER DEPTH
# ===========================================================================

# ---------------------------------------------------------------------
# Hap1: load per-sample chr15 k-mer depth tables, normalize each
# sample's depth by its own median (to control for coverage differences
# between samples), and stack into one long-format table
# ---------------------------------------------------------------------
setwd("hap1/")
lf <- list.files(pattern = ".chr15.tsv")
samples_full <- gsub(".hap1.chr15.tsv", "", lf)
samples_test = c("nmr_1", "nmr_2", "nmr_3", "nmr_4", "nmr_5")
samples <- samples_full
t <- NULL
index = 0
for (s in samples) {
  for (h in c("hap1")) {
    
    fn <- sprintf("%s.%s.chr15.tsv", s, h)
    data <- read.table(fn, col.names = c("chromosome", "start", "end", "unique_kmer_count", "unique_kmer_covered", "unique_fraction_covered", "kmers_used", "kmer_depth"))
    m <- median(data$kmer_depth)
    data$normalized_depth = data$kmer_depth / m
    data$sample = s
    index = index + 1
    t = rbind(t, data)
  }
}

# ---------------------------------------------------------------------
# Reshape to a sample x genomic-bin matrix of normalized depth, for
# hierarchical clustering of samples by their chr15 depth profile
# ---------------------------------------------------------------------
mat <- t %>%
  mutate(start = as.factor(start)) %>%         # ensure ordered columns
  select(sample, start, normalized_depth) %>%
  tidyr::pivot_wider(names_from = start, values_from = normalized_depth)
rownames_mat <- mat$sample
mat_num <- as.matrix(mat[, -1, drop = FALSE])
rownames(mat_num) <- rownames_mat

# cluster samples using single-linkage on the first 27 Mb of chr15
# (bins in 500kb steps), based on Euclidean distance between depth profiles
d <- dist(mat_num[, as.character(seq(0, 27000000, 500000))], method = "euclidean")
hc <- hclust(d, method = "single")

# Get the order of samples from the dendrogram (used to order the y-axis
# of the heatmap plot so similar samples are grouped together)
sample_order <- rownames(mat_num)[hc$order]
print(sample_order)
dend <- rotate(hc, hc$order)
png("Hap1_kmer_dendogram_plot.png")
plot(dend)
dev.off()

# Attach sample origin metadata to the long-format depth table
t1 <- t %>%
  left_join(origin_key, by = "sample")

# ---------------------------------------------------------------------
# Build the plotting data frame: order samples per clustering, restrict
# to the first 40 Mb of chr15 for display
# ---------------------------------------------------------------------
plot_df <- t1 %>%
  mutate(sample = factor(sample, levels = sample_order))
rf <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
r <- rf(32)
lo = 0
hi = 2
plot_df <- plot_df[plot_df$end <= 40e6, ]

##
##### Try to add sample origin.
##

sample_levels <- sample_order   # length = number of unique samples
# If you're not sure what order ggplot is currently using for plot_df,
# extract it directly so you match the existing plot exactly:
sample_levels <- levels(factor(plot_df$sample))   # only if plot_df$sample is already a factor with the order you want
# OR, if it's not yet a factor and the row order in the data IS your intended order:
# sample_levels <- unique(plot_df$sample)

### 1. Get one origin per sample (should be constant across rows for a given sample)
sample_origin <- plot_df %>%
  distinct(sample, origin)

# 2. Get the y-axis order ggplot will actually use
#    (if sample is already a factor with your fixed order, use its levels;
#     otherwise this falls back to alphabetical, which is ggplot's default)
sample_levels <- levels(factor(plot_df$sample))
# 3. Reorder origin lookup to match that exact order
sample_origin <- sample_origin[match(sample_levels, sample_origin$sample), ]
# 4. Map origin -> color
origin_colors <- setNames(
  scales::hue_pal()(length(unique(sample_origin$origin))),
  unique(sample_origin$origin)
)
# Manually override two origin categories with specific colors (rather
# than relying on the automatic hue palette for these)
origin_colors["kenya-mixed"] <- "#C77CFF"
origin_colors["kenya-lab"]   <- "#00BFC4"
"#7CAE00" # ethiopia-south - 
"#F8766D" # ethiopia-east
# Then re-expand to per-sample labels:
label_colors <- origin_colors[sample_origin$origin]

##
# Heatmap of normalized k-mer depth along chr15 (hap1), samples ordered
# by clustering, y-axis labels colored by sample origin. Full-size and
# "thin" versions saved as both PNG and PDF; legend-less version for
# figure panels, separate version with legend for reference
p <- ggplot(plot_df, aes(start, sample, fill = normalized_depth)) +
  geom_tile() +
  scale_fill_gradientn(colors = r, na.value = "white", breaks = c(lo, hi), limits = c(lo, hi)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = label_colors)) +
  ggtitle(sprintf("Read span profile %s", h)) +
  theme(legend.position = "none")
ggsave("plot.hap1_40Mbp.png", p, width = 15, height = 10)
ggsave("plot.hap1_40Mbp.pdf", p, width = 15, height = 10)
ggsave("plot.hap1_40Mbp_thin.png", p, width = 7.5, height = 5)
ggsave("plot.hap1_40Mbp_thin.pdf", p, width = 7.5, height = 5)

p1 <- ggplot(plot_df, aes(start, sample, fill = normalized_depth)) +
  geom_tile() +
  scale_fill_gradientn(colors = r, na.value = "white", breaks = c(lo, hi), limits = c(lo, hi)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = label_colors)) +
  ggtitle(sprintf("Read span profile %s", h))
ggsave("plot.hap1_40Mbp_legend.pdf", p1, width = 7.5, height = 5)

# Save the hap1 clustering order for reuse later (in the "unique k-mer" section)
sample_order_h1 <- sample_order

##
##### Hap2
##
setwd("/Users/dsokolowski/Library/CloudStorage/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure2/k-mer_analysis/coverage_plots")
setwd("hap2/")
lf <- list.files(pattern = ".chr15.tsv")
samples_full <- gsub(".hap2.chr15.tsv", "", lf)
samples_test = c("nmr_1", "nmr_2", "nmr_3", "nmr_4", "nmr_5")
samples <- samples_full
t <- NULL
index = 0
for (s in samples) {
  for (h in c("hap2")) {
    
    fn <- sprintf("%s.%s.chr15.tsv", s, h)
    data <- read.table(fn, col.names = c("chromosome", "start", "end", "unique_kmer_count", "unique_kmer_covered", "unique_fraction_covered", "kmers_used", "kmer_depth"))
    m <- median(data$kmer_depth)
    data$normalized_depth = data$kmer_depth / m
    data$sample = s
    index = index + 1
    t = rbind(t, data)
  }
}

# convert data frame to matrix (for clustering) - identical procedure to hap1 above
mat <- t %>%
  mutate(start = as.factor(start)) %>%         # ensure ordered columns
  select(sample, start, normalized_depth) %>%
  tidyr::pivot_wider(names_from = start, values_from = normalized_depth)
rownames_mat <- mat$sample
mat_num <- as.matrix(mat[, -1, drop = FALSE])
rownames(mat_num) <- rownames_mat

# cluster
d <- dist(mat_num[, as.character(seq(0, 27000000, 500000))], method = "euclidean")
hc <- hclust(d, method = "single")

# Get the order of samples from the dendrogram (hap2-specific clustering)
sample_order <- rownames(mat_num)[hc$order]
print(sample_order)
sample_order_h2 <- sample_order

######
dend <- rotate(hc, hc$order)
png("hap2_kmer_dendogram_plot.png")
plot(dend)
dev.off()

# ---------------------------------------------------------------------

plot_df <- t %>%
  mutate(sample = factor(sample, levels = sample_order))
rf <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
r <- rf(32)
lo = 0
hi = 2
plot_df <- plot_df[plot_df$end <= 40e6, ]

sample_levels <- sample_order_h2   # length = number of unique samples
sample_levels <- levels(factor(plot_df$sample))
sample_levels <- unique(plot_df$sample)

sample_origin <- plot_df %>%
  distinct(sample, origin)

sample_levels <- levels(factor(plot_df$sample))
sample_origin <- sample_origin[match(sample_levels, sample_origin$sample), ]
origin_colors <- setNames(
  scales::hue_pal()(length(unique(sample_origin$origin))),
  unique(sample_origin$origin)
)
origin_colors["kenya-mixed"] <- "#C77CFF"
origin_colors["kenya-lab"]   <- "#00BFC4"
label_colors <- origin_colors[sample_origin$origin]

##
# Same heatmap as hap1, but for hap2 depth data
p <- ggplot(plot_df, aes(start, sample, fill = normalized_depth)) +
  geom_tile() +
  scale_fill_gradientn(colors = r, na.value = "white", breaks = c(lo, hi), limits = c(lo, hi)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = label_colors)) +
  ggtitle(sprintf("Read span profile %s", h)) +
  theme(legend.position = "none")
ggsave("plot.hap2_40Mbp.png", p, width = 15, height = 10)
ggsave("plot.hap2_40Mbp.pdf", p, width = 15, height = 10)
ggsave("plot.hap2_40Mbp_thin.png", p, width = 7.5, height = 5)
ggsave("plot.hap2_40Mbp_thin.pdf", p, width = 7.5, height = 5)

p1 <- ggplot(plot_df, aes(start, sample, fill = normalized_depth)) +
  geom_tile() +
  scale_fill_gradientn(colors = r, na.value = "white", breaks = c(lo, hi), limits = c(lo, hi)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = label_colors)) +
  ggtitle(sprintf("Read span profile %s", h))
ggsave("plot.hap2_40Mbp_legend.pdf", p1, width = 7.5, height = 5)

# ===========================================================================
# SECTION 2: Hap1 vs. Hap2 - UNIQUE K-MER FRACTION COVERED
# (Same overall structure as Section 1, but using "unique_fraction_covered"
# instead of "normalized_depth" as the metric, and reusing the hap1/hap2
# sample orderings computed in Section 1 rather than re-clustering)
# ===========================================================================

setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure2/k-mer_analysis/coverage_plots/")
samples_full = c("nmr_1","nmr_10","nmr_11","nmr_12","nmr_13","nmr_14",
                 "nmr_1901_10","nmr_1901_11","nmr_1901_12","nmr_1901_13",
                 "nmr_1901_14","nmr_1901_15","nmr_1901_16","nmr_1901_19",
                 "nmr_1901_20","nmr_1901_21","nmr_1901_24","nmr_1901_25",
                 "nmr_1901_3","nmr_1901_5","nmr_1901_7","nmr_1901_9","nmr_2",
                 "nmr_3","nmr_4","nmr_5","nmr_6","nmr_7","nmr_8","nmr_9")

# ---------------------------------------------------------------------
# Hap1: reload chr15 tables (identical loading logic to Section 1) to
# rebuild "t" with unique_fraction_covered available
# ---------------------------------------------------------------------
setwd("hap1/")
lf <- list.files(pattern = ".chr15.tsv")
samples_full <- gsub(".hap1.chr15.tsv", "", lf)
samples_test = c("nmr_1", "nmr_2", "nmr_3", "nmr_4", "nmr_5")
samples <- samples_full
t <- NULL
index = 0
for (s in samples) {
  for (h in c("hap1")) {
    
    fn <- sprintf("%s.%s.chr15.tsv", s, h)
    data <- read.table(fn, col.names = c("chromosome", "start", "end", "unique_kmer_count", "unique_kmer_covered", "unique_fraction_covered", "kmers_used", "kmer_depth"))
    m <- median(data$kmer_depth)
    data$normalized_depth = data$kmer_depth / m
    data$sample = s
    index = index + 1
    t = rbind(t, data)
  }
}

# convert data frame to matrix (for clustering)
mat <- t %>%
  mutate(start = as.factor(start)) %>%         # ensure ordered columns
  select(sample, start, unique_fraction_covered) %>%
  tidyr::pivot_wider(names_from = start, values_from = unique_fraction_covered)
rownames_mat <- mat$sample
mat_num <- as.matrix(mat[, -1, drop = FALSE])
rownames(mat_num) <- rownames_mat

# cluster
d <- dist(mat_num[, as.character(seq(0, 27000000, 500000))], method = "euclidean")
hc <- hclust(d, method = "single")

# Get the order of samples from the dendrogram
# use the hap1 depth-based ordering so both heatmaps show samples in
# the same order
# sample_order <- rownames(mat_num)[hc$order]
sample_order <- sample_order_h1
print(sample_order)
dend <- rotate(hc, hc$order)
png("Hap1_kmer_dendogram_plot_unique.png")
plot(dend)
dev.off()

######### Plot with origin
t1 <- t %>%
  left_join(origin_key, by = "sample")

# plot
plot_df <- t1 %>%
  mutate(sample = factor(sample, levels = sample_order))
rf <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
r <- rf(32)
lo = 0
hi = 1
plot_df <- plot_df[plot_df$end <= 40e6, ]

##
##### Try to add sample origin.
##
sample_levels <- sample_order   # length = number of unique samples
sample_levels <- levels(factor(plot_df$sample))
sample_levels <- unique(plot_df$sample)

### 1. Get one origin per sample (should be constant across rows for a given sample)
sample_origin <- plot_df %>%
  distinct(sample, origin)

sample_levels <- levels(factor(plot_df$sample))
sample_origin <- sample_origin[match(sample_levels, sample_origin$sample), ]
origin_colors <- setNames(
  scales::hue_pal()(length(unique(sample_origin$origin))),
  unique(sample_origin$origin)
)
origin_colors["kenya-mixed"] <- "#C77CFF"
origin_colors["kenya-lab"]   <- "#00BFC4"
label_colors <- origin_colors[sample_origin$origin]

##
# Heatmap of unique k-mer fraction covered along chr15 (hap1)
p <- ggplot(plot_df, aes(start, sample, fill = unique_fraction_covered)) +
  geom_tile() +
  scale_fill_gradientn(colors = r, na.value = "white", breaks = c(lo, hi), limits = c(lo, hi)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = label_colors)) +
  ggtitle(sprintf("k-mer profile %s", h)) +
  theme(legend.position = "none")
ggsave("plot.hap1_40Mbp_unique.png", p, width = 15, height = 10)
ggsave("plot.hap1_40Mbp_unique.pdf", p, width = 15, height = 10)
ggsave("plot.hap1_40Mbp_thin_unique.png", p, width = 7.5, height = 5)
ggsave("plot.hap1_40Mbp_thin_unique.pdf", p, width = 7.5, height = 5)

p1 <- ggplot(plot_df, aes(start, sample, fill = unique_fraction_covered)) +
  geom_tile() +
  scale_fill_gradientn(colors = r, na.value = "white", breaks = c(lo, hi), limits = c(lo, hi)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = label_colors)) +
  ggtitle(sprintf("k-mer profile %s", h))
ggsave("plot.hap1_40Mbp_legend_unique.pdf", p1, width = 7.5, height = 5)

###
sample_order_h1 <- sample_order

##
##### Hap2
##
setwd("/Users/dsokolowski/Library/CloudStorage/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure2/k-mer_analysis/coverage_plots")
setwd("hap2/")
lf <- list.files(pattern = ".chr15.tsv")
samples_full <- gsub(".hap2.chr15.tsv", "", lf)
samples_test = c("nmr_1", "nmr_2", "nmr_3", "nmr_4", "nmr_5")
samples <- samples_full
t <- NULL
index = 0
for (s in samples) {
  for (h in c("hap2")) {
    
    fn <- sprintf("%s.%s.chr15.tsv", s, h)
    data <- read.table(fn, col.names = c("chromosome", "start", "end", "unique_kmer_count", "unique_kmer_covered", "unique_fraction_covered", "kmers_used", "kmer_depth"))
    m <- median(data$kmer_depth)
    data$normalized_depth = data$kmer_depth / m
    data$sample = s
    index = index + 1
    t = rbind(t, data)
  }
}

# convert data frame to matrix (for clustering)
mat <- t %>%
  mutate(start = as.factor(start)) %>%         # ensure ordered columns
  select(sample, start, unique_fraction_covered) %>%
  tidyr::pivot_wider(names_from = start, values_from = unique_fraction_covered)
rownames_mat <- mat$sample
mat_num <- as.matrix(mat[, -1, drop = FALSE])
rownames(mat_num) <- rownames_mat

# cluster
d <- dist(mat_num[, as.character(seq(0, 27000000, 500000))], method = "euclidean")
hc <- hclust(d, method = "single")

# Get the order of samples from the dendrogram
# sample_order <- rownames(mat_num)[hc$order]
sample_order <- sample_order_h2   # reuse hap2 depth-based ordering for consistency
print(sample_order)
dend <- rotate(hc, hc$order)
png("hap2_kmer_dendogram_plot_uni.png")
plot(dend)
dev.off()

######### Plot with origin
t1 <- t %>%
  left_join(origin_key, by = "sample")

# plot
plot_df <- t1 %>%
  mutate(sample = factor(sample, levels = sample_order))
rf <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
r <- rf(32)
lo = 0
hi = 1
plot_df <- plot_df[plot_df$end <= 40e6, ]

##
##### Try to add sample origin.
##
sample_levels <- sample_order   # length = number of unique samples
sample_levels <- levels(factor(plot_df$sample))
sample_levels <- unique(plot_df$sample)

### 1. Get one origin per sample (should be constant across rows for a given sample)
sample_origin <- plot_df %>%
  distinct(sample, origin)

sample_levels <- levels(factor(plot_df$sample))
sample_origin <- sample_origin[match(sample_levels, sample_origin$sample), ]
origin_colors <- setNames(
  scales::hue_pal()(length(unique(sample_origin$origin))),
  unique(sample_origin$origin)
)
origin_colors["kenya-mixed"] <- "#C77CFF"
origin_colors["kenya-lab"]   <- "#00BFC4"
label_colors <- origin_colors[sample_origin$origin]

##
# Heatmap of unique k-mer fraction covered along chr15 (hap2)
p <- ggplot(plot_df, aes(start, sample, fill = unique_fraction_covered)) +
  geom_tile() +
  scale_fill_gradientn(colors = r, na.value = "white", breaks = c(lo, hi), limits = c(lo, hi)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = label_colors)) +
  ggtitle(sprintf("k-mer profile %s", h)) +
  theme(legend.position = "none")
ggsave("plot.hap2_40Mbp_unique.png", p, width = 15, height = 10)
ggsave("plot.hap2_40Mbp_unique.pdf", p, width = 15, height = 10)
ggsave("plot.hap2_40Mbp_thin_unique.png", p, width = 7.5, height = 5)
ggsave("plot.hap2_40Mbp_thin_unique.pdf", p, width = 7.5, height = 5)

p1 <- ggplot(plot_df, aes(start, sample, fill = unique_fraction_covered)) +
  geom_tile() +
  scale_fill_gradientn(colors = r, na.value = "white", breaks = c(lo, hi), limits = c(lo, hi)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(colour = label_colors)) +
  ggtitle(sprintf("k-mer profile %s", h))
ggsave("plot.hap2_40Mbp_legend_unique.pdf", p1, width = 7.5, height = 5)