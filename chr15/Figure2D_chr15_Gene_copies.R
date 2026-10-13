# ===========================================================================
# chr15 "supergene" orthogroup analysis for mHetGla1.1 pri/alt haplotypes
# Compares OrthoFinder results run on two different ID-naming schemes,
# filters orthogroups localized to a chr15 region of interest, cross-
# references pseudogene calls, and produces summary plots/heatmaps plus
# targeted GFF/FASTA extraction for specific genes of interest (Antxrl,
# Adgrv1, Mipol1, etc.)
# ===========================================================================


# ---------------------------------------------------------------------
# Main vs. ID analysis.
# Load OrthoFinder results for version "i" from two separate OrthoFinder
# runs: one using the main naming scheme ("v090326_results"), one using
# an alternate ID-based naming scheme ("v070326_id_results")
# ---------------------------------------------------------------------
i = 2
# ---------------------------------------------------------------------
# Originally there was version 1-3, where 1 was all genes, 2 was genes with ORF and decent score
# and 3 was genes whose best model was identified from stringtie. After some testing version 2 was the clear best
# for identifying interpretable orthogroups 


# --- Load "main" version OrthoFinder results ---
setwd(paste0("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure2/Figure2DE_chr15_genes/v090326_results/version", i))
orthogene_counts <- read.table(paste0("version", i, "_Orthogroups.GeneCount.tsv"), header = T, as.is = T, sep = "\t")
print(nrow(orthogene_counts))
orthogene <- read.table(paste0("version", i, "_Orthogroups.tsv"), header = T, as.is = T, sep = "\t")
rownames(orthogene) <- orthogene$Orthogroup

# --- Load "ID" version OrthoFinder results (used downstream) ---
setwd(paste0("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure2/Figure2DE_chr15_genes/v070326_id_results/version", i))
orthogene_counts <- read.table(paste0("version", i, "_Orthogroups.GeneCount.tsv"), header = T, as.is = T, sep = "\t")
print(nrow(orthogene_counts))
orthogene <- read.table(paste0("version", i, "_Orthogroups.tsv"), header = T, as.is = T, sep = "\t")
rownames(orthogene) <- orthogene$Orthogroup
rownames(orthogene_counts) <- orthogene_counts$Orthogroup

## ---------------------------------------------------------------------
## Get chr15 orthogenes
##
## Pri mikado gene IDs go up to: mikado.chr15G978
## Alt mikado gene IDs go up to: mikado.chr15G918
## Keep orthogroups where >50% of the orthogroup's genes (for either pri
## or alt) fall within that chr15 gene-ID range - i.e. orthogroups
## enriched for genes physically located in the chr15 region of interest
## ---------------------------------------------------------------------
library(data.table)
library(stringr)
library(Biostrings)
dt <- fread(paste0("version", i, "_Orthogroups.tsv"))
#rownames(dt) <- dt$Orthogroup

# Helper: given a character vector of comma-separated mikado gene ID
# lists (one string per orthogroup), determine for each orthogroup
# whether more than 50% of its genes have a chr15G<N> ID with N <= max_gene
check_column <- function(x, max_gene) {
  res <- logical(length(x))
  for (i in seq_along(x)) {
    if (is.na(x[i]) || x[i] == "") {
      res[i] <- FALSE
      next
    }
    genes <- str_split(x[i], ",")[[1]]
    genes <- str_trim(genes)
    
    # Extract the numeric gene index from IDs like "...chr15G123"
    nums <- str_extract(genes, "(?<=chr15G)\\d+")
    nums <- str_extract(nums, "\\d+")
    nums <- as.numeric(nums)
    
    in_range <- sum(nums >= 1 & nums <= max_gene, na.rm = TRUE)
    res[i] <- (in_range / length(genes)) > 0.5
  }
  res
}

# Flag orthogroups as chr15-region candidates based on alt and pri gene ID ranges
dt[, alt_condition := check_column(mHetGla.alt.mikadoID.version2, 918)]
dt[, pri_condition := check_column(mHetGla.pri.mikadoID.version2, 978)]

# Keep orthogroups flagged by either haplotype's condition
filtered <- dt[alt_condition | pri_condition]
fwrite(filtered, "chr15_filtered_orthogroups_geneID.txt", sep = "\t")

# ---------------------------------------------------------------------
# Name the ones where we can.
# Translate mikado gene IDs within each filtered orthogroup into
# human-readable gene symbols ("unique_gene"), using per-haplotype
# symbol/metrics lookup tables
# ---------------------------------------------------------------------
setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure2/Figure2DE_chr15_genes")

pri_symbol_metrics <- read.table("mHetGla1.1.pri_gene_symbols_metrics.tsv", header = T, as.is = T, sep = "\t")
rownames(pri_symbol_metrics) <- pri_symbol_metrics$mikado_id

alt_symbol_metrics <- read.table("mHetGla.alt_gene_symbols_metrics.tsv", header = T, as.is = T, sep = "\t")
rownames(alt_symbol_metrics) <- alt_symbol_metrics$mikado_id

# Split pri mikado IDs per orthogroup, strip whitespace, then map each
# ID to its gene symbol via the pri lookup table
filtered_pri <- strsplit(filtered$mHetGla.pri.mikadoID.version2, ",")
filtered_pri1 <- lapply(filtered_pri, function(x) gsub(" ", "", x))
filtered_pri_name <- lapply(filtered_pri1, function(x) {
  y <- c()
  for (i in 1:length(x)) {
    y[i] <- pri_symbol_metrics[x[i], "unique_gene"]
  }
  return(y)
})

# Same process for alt mikado IDs
filtered_alt <- strsplit(filtered$mHetGla.alt.mikadoID.version2, ",")
filtered_alt1 <- lapply(filtered_alt, function(x) gsub(" ", "", x))
filtered_alt_name <- lapply(filtered_alt1, function(x) {
  y <- c()
  for (i in 1:length(x)) {
    y[i] <- alt_symbol_metrics[x[i], "unique_gene"]
  }
  return(y)
})

# Collapse the per-gene symbol lists back into comma-separated strings,
# one per orthogroup, and add as new columns
filtered$mHetGla.alt.version2 <- unlist(lapply(filtered_alt_name, function(x) paste0(x, collapse = ",")))
filtered$mHetGla.pri.version2 <- unlist(lapply(filtered_pri_name, function(x) paste0(x, collapse = ",")))

# Merge in orthogroup gene counts, keep a simplified column subset, and
# save the combined chr15 supergene orthogroup table
filtered_counts <- orthogene_counts[filtered$Orthogroup, ]
filtered_simple <- filtered[, c(1, 2, 3, 4, 9, 10)]
filtered_orthogroup <- merge(filtered_simple, filtered_counts, by = "Orthogroup")

write.table(filtered_orthogroup, file = "chr15_supergene_pri_alt.txt", quote = F, row.names = T, col.names = T, sep = "\t")
saveRDS(filtered_orthogroup, "chr15_supergene_pri_alt.rds")


# ---------------------------------------------------------------------
# Load the manually curated version of the filtered chr15 orthogroup table
# ---------------------------------------------------------------------
named_filtered_orthogroups <- read.delim("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure2/Figure2DE_chr15_genes/chr15_orthgroups_v13326.txt")

# Sanity checks: these sums should be non-trivial only if the
# "noPseudogenes.counts" columns already exist in this file (they are
# only computed later in this script) - if not, R will throw an error
# indicating these columns are undefined, so this only works because
# named_filtered_orthogroups already includes those columns from a prior run
sum(named_filtered_orthogroups$mHetGla.pri.mikadoID.version2 - named_filtered_orthogroups$mHetGla.pri.noPseudogenes.counts)
sum(named_filtered_orthogroups$mHetGla.alt.mikadoID.version2 - named_filtered_orthogroups$mHetGla.alt.noPseudogenes.counts)

# ---------------------------------------------------------------------
# Get pseudogenes from pGre (pseudogene re-annotation pipeline output)
# for the pri haplotype, and compute, per orthogroup, how many of its
# pri genes are NOT pseudogenes
# ---------------------------------------------------------------------
mHetGla1.pri_pseudogenes <- read.table("mHetGla1.pri_annotation_metrics_with_pseudogenes_v120426.tsv.gz", header = T, as.is = T, sep = "\t")

pri_noPseudo <- c()
pri_noPseudo_count <- c()
for (i in 1:nrow(named_filtered_orthogroups)) {
  genes <- strsplit(named_filtered_orthogroups$mHetGla.pri.version2[i], ",")[[1]]
  mHetGla1.pri_pseudogenes_orthogene <- mHetGla1.pri_pseudogenes[mHetGla1.pri_pseudogenes$Name %in% genes, ]
  mHetGla1.pri_pseudogenes_orthogene <- mHetGla1.pri_pseudogenes_orthogene[!duplicated(mHetGla1.pri_pseudogenes_orthogene$Name), ]
  gene_nopseudo <- mHetGla1.pri_pseudogenes_orthogene$Name[mHetGla1.pri_pseudogenes_orthogene$gene_type != "pseudogene"]
  pri_noPseudo_count[i] <- length(gene_nopseudo)
  pri_noPseudo[i] <- paste0(gene_nopseudo, collapse = ",")
}

# Simplify orthogroup labels: take only the first symbol before ";"
# (some orthogroup "Name" fields contain multiple semicolon-separated labels)
label <- named_filtered_orthogroups$Name
label <- unlist(lapply(strsplit(label, ";"), function(x) return(x[1])))

# Build a data frame comparing total pri gene count per orthogroup vs.
# the count excluding pseudogenes, for a scatter/bar comparison
df2 <- as.data.frame(cbind(named_filtered_orthogroups$mHetGla.pri.mikadoID.version2,
                           pri_noPseudo_count))
df2$V1 <- as.numeric(df2$V1)
df2$pri_noPseudo_count <- as.numeric(df2$pri_noPseudo_count)
colnames(df2) <- c("Orthogroup_Count", "Orthogroup_noPseudo")
#rownames(df2) <- label

# Ratio of non-pseudogene count to total count (1 = no pseudogenes in group)
ratio <- df2$Orthogroup_noPseudo / df2$Orthogroup_Count
ratio[is.na(ratio)] <- 1
names(ratio) <- label
table(ratio == 1)

# Bar plot of orthogroups where the ratio is below 1 (i.e., contain
# at least one pseudogene) - shows how "pseudogene-diluted" each group is
ratio1 <- ratio[ratio < 1]
barplot(ratio1, horiz = TRUE, las = 2)

#df2$Orthogene[is.na(df2$Orthogene)] <- ""

# Scatter plot: total orthogroup gene count vs. non-pseudogene count,
# labeled by orthogroup name (pri haplotype)
library(ggplot2); library(ggrepel)
g2 <- ggplot(df2, aes(Orthogroup_Count, Orthogroup_noPseudo, label = Orthogene)) +
  geom_point() + theme_bw() +
  geom_text_repel(max.overlaps = 100) + ggtitle("mHetGla1.1 Pri")
#rownames(df2) <- named_filtered_orthogroup

# Base R equivalent scatter plots (with/without label argument, which
# base plot() ignores - label is not a valid plot() parameter, so this
# is effectively identical to the second call)
plot(named_filtered_orthogroups$mHetGla.pri.mikadoID.version2, pri_noPseudo_count,
     xlab = "chr15 Orthogene count", ylab = "chr15 Orthogene count noPseudogene",
     main = "mHetGla1.1pri - chr15", label = label)
plot(named_filtered_orthogroups$mHetGla.pri.mikadoID.version2, pri_noPseudo_count,
     xlab = "chr15 Orthogene count", ylab = "chr15 Orthogene count noPseudogene",
     main = "mHetGla1.1pri - chr15")

# ---------------------------------------------------------------------
# Same analysis but for alt haplotype
# ---------------------------------------------------------------------
mHetGla1.alt_pseudogenes <- read.table("mHetGla1.alt_annotation_metrics_with_pseudogenes_v050526.tsv", header = T, as.is = T, sep = "\t")

alt_noPseudo <- c()
alt_noPseudo_count <- c()
for (i in 1:nrow(named_filtered_orthogroups)) {
  genes <- strsplit(named_filtered_orthogroups$mHetGla.alt.version2[i], ",")[[1]]
  mHetGla1.alt_pseudogenes_orthogene <- mHetGla1.alt_pseudogenes[mHetGla1.alt_pseudogenes$Name %in% genes, ]
  mHetGla1.alt_pseudogenes_orthogene <- mHetGla1.alt_pseudogenes_orthogene[!duplicated(mHetGla1.alt_pseudogenes_orthogene$Name), ]
  gene_nopseudo <- mHetGla1.alt_pseudogenes_orthogene$Name[mHetGla1.alt_pseudogenes_orthogene$gene_type != "pseudogene"]
  alt_noPseudo_count[i] <- length(gene_nopseudo)
  alt_noPseudo[i] <- paste0(gene_nopseudo, collapse = ",")
}

label <- named_filtered_orthogroups$Name
label <- unlist(lapply(strsplit(label, ";"), function(x) return(x[1])))
# Here, unlike the pri version above, labels are de-duplicated by
# appending a row index suffix (needed since rownames(df2) <- label
# is used below, and rownames must be unique)
label <- paste0(label, "-", 1:length(label))

df2 <- as.data.frame(cbind(named_filtered_orthogroups$mHetGla.alt.mikadoID.version2,
                           alt_noPseudo_count, label))
df2$V1 <- as.numeric(df2$V1)
df2$alt_noPseudo_count <- as.numeric(df2$alt_noPseudo_count)
colnames(df2) <- c("Orthogroup_Count", "Orthogroup_noPseudo", "Orthogene")
rownames(df2) <- label

ratio <- df2$Orthogroup_noPseudo / df2$Orthogroup_Count
ratio[is.na(ratio)] <- 1
names(ratio) <- label
table(ratio == 1)
ratio1 <- ratio[ratio < 1]
barplot(ratio1, horiz = TRUE, las = 2)

df2$Orthogene[is.na(df2$Orthogene)] <- ""

library(ggplot2); library(ggrepel)
g2 <- ggplot(df2, aes(Orthogroup_Count, Orthogroup_noPseudo, label = Orthogene)) +
  geom_point() + theme_bw() +
  geom_text_repel(max.overlaps = 100) + ggtitle("mHetGla1.1 alt")

plot(named_filtered_orthogroups$mHetGla.alt.mikadoID.version2, alt_noPseudo_count,
     xlab = "chr15 Orthogene count", ylab = "chr15 Orthogene count noPseudogene",
     main = "mHetGla1.1alt - chr15", label = label)
plot(named_filtered_orthogroups$mHetGla.alt.mikadoID.version2, alt_noPseudo_count,
     xlab = "chr15 Orthogene count", ylab = "chr15 Orthogene count noPseudogene",
     main = "mHetGla1.1alt - chr15")

##
###### Add in counts for no pseudogenes
###

# Attach the computed non-pseudogene gene lists/counts (both haplotypes)
# back onto the main orthogroup table
named_filtered_orthogroups$mHetGla.pri.noPseudogenes <- pri_noPseudo
named_filtered_orthogroups$mHetGla.pri.noPseudogenes.counts <- pri_noPseudo_count
named_filtered_orthogroups$mHetGla.alt.noPseudogenes <- alt_noPseudo
named_filtered_orthogroups$mHetGla.alt.noPseudogenes.counts <- alt_noPseudo_count

# ---------------------------------------------------------------------
# Aggregate non-pseudogene counts by orthogroup Name (in case multiple
# rows share the same Name), for pri and alt separately, then merge
# into a single pri-vs-alt copy-number comparison table
# ---------------------------------------------------------------------
pri_aggregated <- aggregate(named_filtered_orthogroups$mHetGla.pri.noPseudogenes.counts ~ named_filtered_orthogroups$Name, FUN = "sum")
alt_aggregated <- aggregate(named_filtered_orthogroups$mHetGla.alt.noPseudogenes.counts ~ named_filtered_orthogroups$Name, FUN = "sum")
colnames(alt_aggregated) <- c("name", "alt")
colnames(pri_aggregated) <- c("name", "pri")

merged <- merge(alt_aggregated, pri_aggregated, by = "name")
rownames(merged) <- merged$name
merged <- merged[, c("pri", "alt")]
merged <- merged[order(rowSums(merged), decreasing = T), ]

# Split into genes with real symbols ("noID") vs. genes still labeled
# by chromosome/mikado-style placeholder IDs ("ID")
merged_noID <- merged[!grepl("chr", rownames(merged)), ]
# Restrict to genes with an actual pri/alt copy-number discrepancy
merged_no_even <- merged_noID[merged_noID$pri != merged_noID$alt, ]
merged_no_even <- merged_no_even[!grepl("ENSHGLG",rownames(merged_no_even)),]
# ---------------------------------------------------------------------
# Heatmaps: pri vs alt copy number (non-pseudogene) for named genes
# with uneven copy number between haplotypes
# ---------------------------------------------------------------------
library(pheatmap)
my_colors <- colorRampPalette(c("white", "red"))(100)

pdf(file = "chr15_uneven_symbol_copynumber_nopseudo_v20260619.pdf")
print(pheatmap(t(log2(merged_no_even + 1)), display_numbers = t(merged_no_even), fontsize = 12, cluster_rows = FALSE,
               cluster_cols = FALSE, fontsize_number = 14, cellwidth = 20, , cellheight = 20, color = my_colors))
dev.off()

png(file = "chr15_uneven_symbol_copynumber_nopseudo_v20260506.png", width = 8, height = 11, units = "in", res = 300)
print(pheatmap(log2(merged_no_even + 1), display_numbers = merged_no_even, fontsize = 12, cluster_rows = FALSE,
               cluster_cols = FALSE, fontsize_number = 14, cellwidth = 50, color = my_colors))
dev.off()

# Full (all named genes, not just uneven ones) heatmap, dropping rows
# where both pri and alt copy number are zero
merged_noID <- merged_noID[rowSums(merged_noID) > 0, ]
png(file = "chr15_symbol_copynumber_nopseudo_v20260506.png", width = 8, height = 11, units = "in", res = 300)
print(pheatmap(merged_noID[rowSums(merged_noID) > 0, ], display_numbers = merged_noID, fontsize = 6, cluster_rows = FALSE,
               cluster_cols = FALSE, fontsize_number = 6, cellwidth = 10, color = my_colors))
dev.off()

# Same heatmap, but for genes still identified only by chromosome/ID
# (i.e., unnamed genes), restricted to non-zero rows
mergedID <- merged[grepl("chr", rownames(merged)), ]
mergedID <- mergedID[rowSums(mergedID) > 0, ]
png(file = "chr15_chr_copynumber_nopseudo_v20260506.png", width = 8, height = 11, units = "in", res = 300)
print(pheatmap(mergedID, display_numbers = mergedID, fontsize = 6, cluster_rows = FALSE,
               cluster_cols = FALSE, fontsize_number = 6, cellwidth = 10, color = my_colors))
dev.off()

# Uneven copy-number heatmap including chr-ID-only genes (full set,
# not restricted to named genes only)
merged_noeven_withID <- merged[merged$pri != merged$alt, ]
merged_nochr <- merged[grepl("chr", rownames(merged)), ]
pdf(file = "chr15_uneven_copy_number_nopseudo.pdf", width = 22)
print(pheatmap(log2(merged_noeven_withID + 1), display_numbers = merged_noeven_withID))
dev.off()

# ---------------------------------------------------------------------
# GO enrichment (via g:Profiler) on the named-gene sets: all named
# genes vs. those with uneven pri/alt copy number
# ---------------------------------------------------------------------
noID <- gprofiler2::gost(rownames(merged_noID), organism = "mmusculus", evcodes = T, significant = F)
noEven <- gprofiler2::gost(rownames(merged_no_even), organism = "mmusculus", evcodes = T)

##
##### Antrxl (Antxrl) analysis
##
# setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Supplementary Figures/chr15_genes_and_antrx/")

alt_aa <- readAAStringSet("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Supplementary Figures/chr15_genes_and_antrx/mHetGla.alt.version2.faa")
pri_aa <- readAAStringSet("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Supplementary Figures/chr15_genes_and_antrx/mHetGla.pri.version2.faa")

# Pull out the specific non-pseudogene mikado IDs for the Antxrl-related
# orthogroup entries (Gm5460 = alt symbol, Antxrl = pri symbol), extract
# their protein sequences, and write out per-haplotype FASTA files
alt_antrxl <- named_filtered_orthogroups[named_filtered_orthogroups$Name == "Gm5460", "mHetGla.alt.noPseudogenes"][[1]]
pri_antrxl <- named_filtered_orthogroups[named_filtered_orthogroups$Name %in% "Antxrl", "mHetGla.pri.noPseudogenes"][[1]]

alt_Antxrl3_aa <- alt_aa[strsplit(alt_antrxl, ",")[[1]]]
pri_Antxrl_aa <- pri_aa[strsplit(pri_antrxl, ",")[[1]]]

writeXStringSet(alt_Antxrl3_aa, file = "alt_Antxrl3.faa")
writeXStringSet(pri_Antxrl_aa, file = "pri_Antxrl.faa")

## ---------------------------------------------------------------------
## Load in InterProScan results for pri Antxrl proteins, restrict to
## Pfam domain hits, and plot a presence/absence heatmap of domains
## per protein
## ---------------------------------------------------------------------
pri_Antxrl_interpo <- read.table("AntrxlMipol1Adgrv1/pri_Antxrl.faa.tsv", header = FALSE, as.is = T, sep = "\t")
pri_Antxrl_interpo_pfam <- pri_Antxrl_interpo[pri_Antxrl_interpo$V4 == "Pfam", ]
domains <- pri_Antxrl_interpo_pfam[, c("V1", "V6")]
tab_domains <- table(domains[, 1], domains[, 2])
dim(tab_domains)

colnames(tab_domains) <- c("Antrx_ECM", "vWFa")
my_colors <- colorRampPalette(c("white", "black"))(2)
rownames(tab_domains) <- gsub("mikado.", "", rownames(tab_domains))

pdf(file = "pri_Antrxl_domains.pdf")
print(pheatmap(tab_domains, cluster_rows = F, cluster_cols = F, cellwidth = 20, color = my_colors, fontsize = 6))
dev.off()

# Same domain heatmap, but for the alt haplotype Antxrl3 proteins
alt_Antxrl_interpo <- read.table("AntrxlMipol1Adgrv1/alt_Antxrl3.faa.tsv", header = FALSE, as.is = T, sep = "\t")
alt_Antxrl_interpo_pfam <- alt_Antxrl_interpo[alt_Antxrl_interpo$V4 == "Pfam", ]
domains <- alt_Antxrl_interpo_pfam[, c("V1", "V6")]
tab_domains <- table(domains[, 1], domains[, 2])
dim(tab_domains)
colnames(tab_domains) <- c("Antrx_ECM", "vWFa")
my_colors <- colorRampPalette(c("white", "black"))(2)
rownames(tab_domains) <- gsub("mikado.", "", rownames(tab_domains))

pdf(file = "alt_Antxrl3_domains.pdf")
pheatmap(tab_domains, cluster_rows = F, cluster_cols = F, cellwidth = 20, color = my_colors, fontsize = 6)
dev.off()

# ---------------------------------------------------------------------
# Extract protein sequences for all chr15-region genes (by mikado ID,
# not gene symbol) from the earlier chr15 supergene orthogroup table,
# and write per-haplotype FASTA files
# ---------------------------------------------------------------------
mikado_pri <- unlist(strsplit(filtered_orthogroup$mHetGla.pri.version2, ","))
pri <- mikado_pri[grep("chr15G", mikado_pri)]
mikado_alt <- unlist(strsplit(filtered_orthogroup$mHetGla.alt.version2, ","))
alt <- mikado_alt[grep("chr15G", mikado_alt)]

library(Biostrings)
alt_aa <- readAAStringSet("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Supplementary Figures/chr15_genes_and_antrx/mHetGla.alt.version2.faa")
pri_aa <- readAAStringSet("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Supplementary Figures/chr15_genes_and_antrx/mHetGla.pri.version2.faa")

alt_aa_IDonly <- alt_aa[alt]
pri_aa_IDonly <- pri_aa[pri]
writeXStringSet(pri_aa_IDonly, file = "pri_chr15_gene_ID_only.faa")
writeXStringSet(alt_aa_IDonly, file = "alt_chr15_gene_ID_only.faa")

##
##### Antrxl gene filtering
##
writeXStringSet(pri_aa_IDonly, file = "pri_chr15_gene_symbol_only.faa")
writeXStringSet(alt_aa_IDonly, file = "alt_chr15_gene_symbol_only.faa")

###

# ---------------------------------------------------------------------
# Identify orthogroups with a strongly uneven gene count between pri
# and alt (>3x difference), restricted to orthogroups with a
# meaningful minimum count (>10 in either haplotype) - candidate
# haplotype-specific expansions/collapses
# ---------------------------------------------------------------------
double_group <- named_filtered_orthogroups[
  named_filtered_orthogroups$mHetGla.pri.mikadoID.version2 > 3 * named_filtered_orthogroups$mHetGla.alt.mikadoID.version2 |
    named_filtered_orthogroups$mHetGla.alt.mikadoID.version2 > 3 * named_filtered_orthogroups$mHetGla.pri.mikadoID.version2, ]
double_group <- double_group[double_group$mHetGla.alt.mikadoID.version2 > 10 | double_group$mHetGla.pri.mikadoID.version2 > 10, ]

# Separate one-directional subsets: pri-expanded and alt-expanded orthogroups
double_group_pri <- named_filtered_orthogroups[named_filtered_orthogroups$mHetGla.pri.mikadoID.version2 > 3 * named_filtered_orthogroups$mHetGla.alt.mikadoID.version2, ]
double_group_alt <- named_filtered_orthogroups[named_filtered_orthogroups$mHetGla.alt.mikadoID.version2 > 3 * named_filtered_orthogroups$mHetGla.pri.mikadoID.version2, ]
double_group_pri <- double_group_pri[double_group_pri$mHetGla.pri.mikadoID.version2 > 2, ]

# ---------------------------------------------------------------------
# Pull the actual gene annotation records (from the full genome GFFs)
# for genes belonging to these uneven-copy-number orthogroups,
# restricted to chr15 genes below 26 Mb
# ---------------------------------------------------------------------
library(rtracklayer)
mHetGla1.pri.anno <- readGFF("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/mHetGla1.pri.MT.v020326.gff.gz")

mHetGla1.pri.anno$Name %in% double_group$mHetGla.pri.version2

pri_chr15_genes_OGs_enriched <- unique(unlist(strsplit(double_group$mHetGla.pri.version2, ",")))
pri_chr15_genes_OGs_enriched <- pri_chr15_genes_OGs_enriched[!is.na(pri_chr15_genes_OGs_enriched)]

mHetGla1.pri.anno.OG <- mHetGla1.pri.anno[mHetGla1.pri.anno$Name %in% pri_chr15_genes_OGs_enriched, ]
mHetGla1.pri.anno.OG <- mHetGla1.pri.anno.OG[mHetGla1.pri.anno.OG$seqid == "chr15", ]
mHetGla1.pri.anno.OG <- mHetGla1.pri.anno.OG[mHetGla1.pri.anno.OG$end < 2.6e7, ]
mHetGla1.pri.anno.OG <- mHetGla1.pri.anno.OG[mHetGla1.pri.anno.OG$type == "gene", ]

# Same extraction for the alt haplotype annotation
mHetGla1.alt.anno <- readGFF("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/mHetGla1.1.alt/mHetGla1.alt.edited.geneSymbol.gff.gz")
alt_chr15_genes_OGs_enriched <- unique(unlist(strsplit(double_group$mHetGla.alt.version2, ",")))
alt_chr15_genes_OGs_enriched <- alt_chr15_genes_OGs_enriched[!is.na(alt_chr15_genes_OGs_enriched)]

mHetGla1.alt.anno.OG <- mHetGla1.alt.anno[mHetGla1.alt.anno$Name %in% alt_chr15_genes_OGs_enriched, ]
mHetGla1.alt.anno.OG <- mHetGla1.alt.anno.OG[mHetGla1.alt.anno.OG$seqid == "chr15", ]
mHetGla1.alt.anno.OG <- mHetGla1.alt.anno.OG[mHetGla1.alt.anno.OG$end < 2.6e7, ]
mHetGla1.alt.anno.OG <- mHetGla1.alt.anno.OG[mHetGla1.alt.anno.OG$type == "gene", ]

# Export the full uneven-copy-number gene sets (pri and alt) as GFF3
export(mHetGla1.pri.anno.OG, con = "mHetGla1.pri.anno.OG.gt10.gff", format = "GFF3")
export(mHetGla1.alt.anno.OG, con = "mHetGla1.alt.anno.OG.gt10.gff", format = "GFF3")

# ---------------------------------------------------------------------
# Manually correct one orthogroup's label (row 1) and pull out specific
# gene records of interest by hardcoded row index into double_group:
# row 1 = Antrxl3 (alt), row 2 = Antrxl/Antxrl (pri),
# row 3 = Adgrv1 (alt), row 4 = Mipol1 (pri)
#
# ---------------------------------------------------------------------
double_group$Name[1] <- "Antrxl3"

Antrxlpri <- mHetGla1.pri.anno.OG[mHetGla1.pri.anno.OG$Name %in% strsplit(double_group[2, "mHetGla.pri.version2"], ",")[[1]], ]
Antrxlalt <- mHetGla1.alt.anno.OG[mHetGla1.alt.anno.OG$Name %in% strsplit(double_group[2, "mHetGla.alt.version2"], ",")[[1]], ]

Mipol1pri <- mHetGla1.pri.anno.OG[mHetGla1.pri.anno.OG$Name %in% strsplit(double_group[4, "mHetGla.pri.version2"], ",")[[1]], ]
Mipol1alt <- mHetGla1.alt.anno.OG[mHetGla1.alt.anno.OG$Name %in% strsplit(double_group[4, "mHetGla.alt.version2"], ",")[[1]], ]

Antrxl3pri <- mHetGla1.pri.anno.OG[mHetGla1.pri.anno.OG$Name %in% strsplit(double_group[1, "mHetGla.pri.version2"], ",")[[1]], ]
Antrxl3alt <- mHetGla1.alt.anno.OG[mHetGla1.alt.anno.OG$Name %in% strsplit(double_group[1, "mHetGla.alt.version2"], ",")[[1]], ]

Adgrv1pri <- mHetGla1.pri.anno.OG[mHetGla1.pri.anno.OG$Name %in% strsplit(double_group[3, "mHetGla.pri.version2"], ",")[[1]], ]
Adgrv1alt <- mHetGla1.alt.anno.OG[mHetGla1.alt.anno.OG$Name %in% strsplit(double_group[3, "mHetGla.alt.version2"], ",")[[1]], ]

# Write minimal BED files (seqid, start, end) for each gene set of interest
write.table(Antrxl3alt[c(1, 4, 5)], file = "AntrxlMipol1Adgrv1/Antrxl3alt.bed", quote = F,
            row.names = F, col.names = F, sep = "\t")
write.table(Antrxl3pri[c(1, 4, 5)], file = "AntrxlMipol1Adgrv1/Antrxl3pri.bed", quote = F,
            row.names = F, col.names = F, sep = "\t")

write.table(Antrxlpri[c(1, 4, 5)], file = "AntrxlMipol1Adgrv1/Antrxlpri.bed", quote = F,
            row.names = F, col.names = F, sep = "\t")
write.table(Antrxlalt[c(1, 4, 5)], file = "AntrxlMipol1Adgrv1/Antrxlalt.bed", quote = F,
            row.names = F, col.names = F, sep = "\t")

write.table(Mipol1pri[c(1, 4, 5)], file = "AntrxlMipol1Adgrv1/Mipol1pri.bed", quote = F,
            row.names = F, col.names = F, sep = "\t")

write.table(Adgrv1pri[c(1, 4, 5)], file = "AntrxlMipol1Adgrv1/Adgrv1pri.bed", quote = F,
            row.names = F, col.names = F, sep = "\t")
write.table(Adgrv1alt[c(1, 4, 5)], file = "AntrxlMipol1Adgrv1/Adgrv1alt.bed", quote = F,
            row.names = F, col.names = F, sep = "\t")

# Export full GFF3 records for each gene set of interest.

export(Antrxlpri, con = "Antrxlpri.gff", format = "GFF3")
export(Antrxl3alt, con = "Antrxl3alt.gff", format = "GFF3") 
export(Mipol1pri, con = "Mipol1pri.gff", format = "GFF3")
export(Adgrv1alt, con = "Adgrv1alt.gff", format = "GFF3")

export(Antrxlpri, con = "Antrxlpri.gff", format = "GFF3")
export(Antrxl3alt, con = "Antrxl3alt.gff", format = "GFF3")  
export(Mipol1pri, con = "Mipol1pri.gff", format = "GFF3")
export(Adgrv1alt, con = "Adgrv1alt.gff", format = "GFF3")