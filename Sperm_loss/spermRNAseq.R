# ==============================================================================
# NMR Sperm RNA-seq Analysis: Sub vs. exsub Differential Expression,
# Cross-Species (Mouse vs. NMR) Comparison, and Relaxed-Selection Overlap
#
# Overview:
#   1. Load NMR + mouse sperm RNA-seq count data, normalize (RPKM/CPM)
#   2. QC: sample correlation, PCA
#   3. DE analysis: subordinate (sub) vs. exsub (bre) NMR sperm
#   4. Annotate DE results with gene names/status from custom annotation metrics
#   5. Pathway enrichment (gprofiler2) on subordinate- vs. exsub-biased genes
#   6. Spermatogenesis-specific heatmaps, with mature-sperm marker highlighting
#   7. Sensitivity analysis: repeat DE controlling for sperm/epithelial
#      contamination (Prm2/Krt5 admixture proxy)
#   8. Cross-species comparison: identify 1:1 mouse-NMR orthologs (OrthoFinder),
#      compare expression, DE (with and without exsub samples), and
#      cross-validate edgeR vs. manual t-test fold changes
#   9. Integrate with dN/dS "RELAX" selection analysis to test whether
#      DE genes overlap with relaxed/intensified selection categories
#      (ActivePathways enrichment)
#  10. Non-coding RNA burden comparison between mouse and NMR sperm
# ==============================================================================


# ------------------------------------------------------------------------
# SECTION 0: One-off annotation prep (commented out / not run)
# Adds a clean gene_id column to exon-level GFF entries for featureCounts
# by stripping trailing isoform suffixes from the transcript ID.
# ------------------------------------------------------------------------
## library(rtracklayer)
## gff <- readGFF("~/Desktop/NMR_multi_assembly/v090226/annotations/MT_annotations/mHetGla1.pri.MT.v020326.gff")
## gff$gene_id <- sub("^([^.]*\\.[^.]*)\\..*", "\\1", gff$ID)
## export(gff, con = "mHetGla1.1.pri.exon_geneID.gff", format = "gff3")



# ------------------------------------------------------------------------
# SECTION 1: Setup and libraries
# ------------------------------------------------------------------------
library(RUVSeq)
library(edgeR)
library(pheatmap)
library(gprofiler2)
library(patchwork)
library(ggplot2)
setwd("~/Desktop/NMR_multi_assembly/sperm_rnaseq")


# ------------------------------------------------------------------------
# SECTION 2: Load raw count data (mouse + NMR sperm)
# ------------------------------------------------------------------------
mouseSperm <- read.table("Mousesperm.txt", header = TRUE, as.is = TRUE)
mouse_counts <- mouseSperm[, grep("sorted.bam", colnames(mouseSperm))]

nmrSperm <- read.table("NMRsperm.txt", header = TRUE, as.is = TRUE)
nmr_counts <- nmrSperm[, grep("sorted.bam", colnames(nmrSperm))]
rownames(nmr_counts) <- nmrSperm$Geneid


# ------------------------------------------------------------------------
# SECTION 3: NMR RPKM normalization + low-expression filter
# Keep only genes with median RPKM > 1 across the 5 NMR samples
# ------------------------------------------------------------------------
d <- DGEList(counts = nmr_counts)
d$genes <- data.frame(Length = nmrSperm$Length)
rpkm_norm <- rpkm(d)
rpkm_norm <- rpkm_norm[rowMedians(rpkm_norm) > 1, ]
colnames(rpkm_norm) <- gsub("_sperm_rightepy.sorted.bam", "", colnames(rpkm_norm))


# ------------------------------------------------------------------------
# SECTION 4: Sample QC - correlation heatmap + PCA
# 5 NMR samples: 2 subordinate ("sub"), 3 exsub ("bre")
# ------------------------------------------------------------------------
library(ggplot2)
library(ggfortify)
pca_res <- prcomp(t(rpkm_norm), scale. = TRUE)
library(corrplot)

colnames(rpkm_norm) <- c("BaD7Sub", "CD34Sub", "963ABre", "A32FBre", "8521Bre")

pdf(file = "NMR_sperm_rpkm_correlation.pdf")
print(corrplot(cor(rpkm_norm), addCoef.col = "white", number.cex = 2))
dev.off()

scores <- as.data.frame(pca_res$x)
scores$condition <- c("sub", "sub", "bre", "bre", "bre")
rownames(scores) <- c("BaD7Sub", "CD34Sub", "963ABre", "A32FBre", "8521Bre")
scores$label <- rownames(scores)

# PCA plot with sample labels
p1 <- autoplot(pca_res, data = scores, colour = "condition", size = 4,
               label = TRUE, label.size = 4, label.vjust = 2) +
  theme_bw() +
  theme(text = element_text(size = 14), axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), legend.title = element_text(size = 13),
        legend.text = element_text(size = 11))

# PCA plot without labels (cleaner for figures), custom sub/bre colors
p2 <- autoplot(pca_res, data = scores, colour = "condition", size = 4,
               label = FALSE, label.size = 4, label.vjust = 2) +
  theme_bw() +
  theme(text = element_text(size = 10), axis.title = element_text(size = 10),
        axis.text = element_text(size = 10), legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values = c("bre" = "purple", "sub" = "orange"))

pdf(file = "NMR_sperm_rpkm_pca.pdf"); print(p1); dev.off()
ggsave(filename = "NMR_sperm_rpkm_pca.png", plot = p2, width = 4, height = 3)


# ------------------------------------------------------------------------
# SECTION 5: Differential expression - Subordinate vs. exsub (NMR)
# Standard edgeR GLM workflow: upperquartile norm, common+tagwise dispersion,
# LRT test on the "condition" coefficient (sub = reference level)
# ------------------------------------------------------------------------
nmr_counts_rpkm_gt1 <- nmr_counts[rownames(rpkm_norm), ]
scores$condition <- factor(scores$condition, levels = c("sub", "bre"))
design <- model.matrix(~condition, data = scores)

y <- DGEList(counts = nmr_counts_rpkm_gt1, group = scores$condition)
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
top <- topTags(lrt, n = nrow(nmr_counts_rpkm_gt1))$table


# ------------------------------------------------------------------------
# SECTION 6: Annotate DE table with gene names/type/completeness
# Joins DE results to a custom NMR annotation metrics file (includes
# pseudogenes); picks the highest-scoring entry per gene when duplicated
# ------------------------------------------------------------------------
metrics <- read.table(
  "~/Desktop/NMR_multi_assembly/v090226/annotations/mHetGla1.pri_annotation_metrics_with_pseudogenes_v120426.tsv.gz",
  header = TRUE, as.is = TRUE, sep = "\t"
)
met1 <- metrics[order(metrics$score, decreasing = TRUE), ]
met_gene <- met1[!duplicated(met1$parent), ]
rownames(met_gene) <- met_gene$parent
met_gene_in <- met_gene[rownames(top), ]

top$Name <- met_gene_in$Name
top$Status <- met_gene_in$gene_type
table(top$FDR < 0.05 & abs(top$logFC > 2))
top$is_complete <- met_gene_in$is_complete
top$min_score <- met_gene_in$minScore

library(dplyr)
df <- top

# Classify each gene by significance direction (bre-biased / sub-biased / NS)
df$significance <- "NS"
df$significance[df$FDR < 0.05 & df$logFC > 1] <- "bre"
df$significance[df$FDR < 0.05 & df$logFC < -1] <- "sub"


# ------------------------------------------------------------------------
# SECTION 7: Volcano plots - full gene set
# ------------------------------------------------------------------------
pv1 <- ggplot(df, aes(x = logFC, y = -log10(FDR), color = significance, shape = Status)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("bre" = "purple", "sub" = "orange", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_classic() +
  labs(x = "log2 Fold Change", y = expression(-log[10]("Adjusted p-value")), color = NULL)

pdf(file = "NMR_sperm_sub_vs_bre_volcano.pdf"); print(pv1); dev.off()
ggsave(filename = "NMR_sperm_sub_vs_bre_volcano.png", plot = pv1, width = 4, height = 4)

pv2 <- pv1 + theme(legend.position = "none")
ggsave(filename = "NMR_sperm_sub_vs_bre_volcano_nolegend.png", plot = pv2, width = 4, height = 4)

saveRDS(df, file = "NMR_sperm_sub_vs_bre.rds")


# ------------------------------------------------------------------------
# SECTION 8: Flag genes with real gene symbols vs. unnamed/mikado/copy IDs
# "Symbol" == TRUE means the gene has a curated name (not a mikado/ENSHGL
# placeholder ID and not a segmental-duplication "-copy" entry)
# ------------------------------------------------------------------------
df$Symbol <- !(grepl("mikado", df$Name) | grepl("ENSHGL", df$Name) | grepl("-copy", df$Name))
df_name <- df[!(grepl("mikado", df$Name) | grepl("ENSHGL", df$Name)), ]
df_noname <- df[(grepl("mikado", df$Name) | grepl("ENSHGL", df$Name)), ]

table(df_name$Status, df_name$significance)
table(df$Symbol, df$significance)
table(df$Status, df$significance)

# Stringent gene set: complete gene models, type == "gene" (not pseudogene/ncRNA),
# and has a real symbol - used for downstream labeling/pathway analysis
df_stringent <- df[df$is_complete == "True" & df$Status == "gene" & df$Symbol == TRUE, ]

# Top 10 most sub-biased + top 10 most bre-biased genes, for labeling on plots
topgenes <- c(df_stringent$Name[df_stringent$logFC < -1][1:10],
              df_stringent$Name[df_stringent$logFC > 1][1:10])
df_stringent$topgenes <- ""
df_stringent$topgenes[df_stringent$Name %in% topgenes] <- df_stringent$Name[df_stringent$Name %in% topgenes]


# ------------------------------------------------------------------------
# SECTION 9: Volcano plots - stringent gene set, with top-gene labels
# ------------------------------------------------------------------------
library(ggrepel)

pv2 <- ggplot(df_stringent, aes(x = logFC, y = -log10(FDR), color = significance, shape = Status)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("bre" = "purple", "sub" = "orange", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_classic() +
  labs(x = "log2 Fold Change", y = expression(-log[10]("Adjusted p-value")), color = NULL) +
  geom_text_repel(data = df_stringent, aes(label = topgenes), size = 3, max.overlaps = Inf) +
  theme(legend.position = "none")

pdf(file = "NMR_sub_vs_bre_stringent.pdf"); print(pv2); dev.off()
ggsave(filename = "NMR_sub_vs_bre_stringent.png", plot = pv2, width = 4, height = 4)

# Same plot, unlabeled, theme_bw variant
pv3 <- ggplot(df_stringent, aes(x = logFC, y = -log10(FDR), color = significance, shape = Status)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("bre" = "purple", "sub" = "orange", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_bw() +
  labs(x = "log2 Fold Change", y = expression(-log[10]("Adjusted p-value")), color = NULL)

ggsave(filename = "NMR_sub_vs_bre_stringent.png", plot = pv3, height = 4, width = 4)
print(pv3); dev.off()

saveRDS(df_stringent, file = "NMR_sperm_sub_vs_bre_stringentgene.rds")


# ------------------------------------------------------------------------
# SECTION 10: Pathway enrichment (gprofiler2) on sub- vs. bre-biased genes
# Background restricted to the stringent gene set itself (not genome-wide)
# ------------------------------------------------------------------------
library(gprofiler2)
library(patchwork)

SUB_DE <- gost(df_stringent$Name[df_stringent$significance == "sub"],
               organism = "mmusculus", custom_bg = df_stringent$Name, evcodes = TRUE)
BRE_DE <- gost(df_stringent$Name[df_stringent$significance == "bre"],
               organism = "mmusculus", custom_bg = df_stringent$Name, evcodes = TRUE)

saveRDS(SUB_DE, file = "NMR_sub_vs_bre_sub_pathways.rds")
saveRDS(BRE_DE, file = "NMR_sub_vs_bre_bre_pathways.rds")

# Top-5 enriched terms, sub- and bre-biased side by side
gp1 <- ggplot(SUB_DE$result %>% head(5),
              aes(x = reorder(term_name, -p_value), y = -log10(p_value), color = source)) +
  geom_col() + coord_flip() + theme_bw() + theme(legend.position = "none") +
  labs(x = "Pathway", y = "-log10(FDR)") +
  scale_color_manual(values = c("GO:BP" = "gold"))

gp2 <- ggplot(BRE_DE$result %>% head(5),
              aes(x = reorder(term_name, -p_value), y = -log10(p_value), color = source)) +
  geom_col() + coord_flip() + theme_bw() + theme(legend.position = "none") +
  labs(x = "", y = "-log10(FDR)") +
  scale_color_manual(values = c("GO:CC" = "purple"))

pdf(file = "NMR_sperm_pathway_enrichment.pdf", width = 16); print(gp1 + gp2); dev.off()
ggsave("NMR_sperm_pathway_enrichment.png", gp1 + gp2, height = 3.5, width = 7)


# ------------------------------------------------------------------------
# SECTION 11: Spermatogenesis-gene heatmap
# Pulls the specific "spermatogenesis" GO term's gene set from the sub-biased
# enrichment result and plots per-sample RPKM
# ------------------------------------------------------------------------
rpkm_norm_stringent <- rpkm_norm[rownames(df_stringent), ]
rownames(rpkm_norm_stringent) <- df_stringent$Name

DEsplit <- strsplit(SUB_DE$result[SUB_DE$result$term_name == "spermatogenesis", "intersection"], ",")[[1]]
sperms <- rpkm_norm_stringent[DEsplit, ]
colnames(sperms) <- c("BaD7Sub", "CD34Sub", "963Abre", "A32Fbre", "8521bre")

pheatmap(t(sperms), scale = "column", fontsize_row = 10, fontsize_col = 10, cellwidth = 10, cellheight = 10)

# Add per-gene median/mean RPKM as a row annotation
annotation_col <- data.frame(RPKM_med = rowMedians(sperms), RPKM_mea = rowMeans(sperms))
p2 <- pheatmap(t(sperms), annotation_col = annotation_col, scale = "column",
               fontsize_row = 8, fontsize_col = 6, cellwidth = 6, cellheight = 6)
pdf(file = "spermatogenesis_pheatmap.pdf", width = 8); print(p2); dev.off()
ggsave(file = "spermatogenesis_pheatmap.pdf", p2, width = 12)

spermatogenesis_genes <- pheatmap(sperms, annotation_row = annotation_col, scale = "row",
                                  fontsize_row = 8, fontsize_col = 8, cellwidth = 10, cellheight = 10)
write.table(rownames(sperms), file = "Sub_spermatogenesis.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Flag genes that also appear in an external mature-sperm marker list
# (from a published snRNA-seq study, PMID 30404016) with an asterisk
mature_sperm <- "Catsperz, Ropn1l, Efcab9, Tssk6, Tnp2, Acrbp, Cabyr, Odf1, Odf2, Spag16, Spag6l, Dnali1, Drc1, Rsph6a, Ttll5, Mns1, Pacrg, Zpbp2, Spata9, Hook1, Ccdc146, Cfap221, Lrrc46, Tesk2"
mature <- read.table(
  "~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Supplementary Figures/sperm_sub_vs_bre/mature_sperm_markers_snrnaseq",
  header = TRUE, as.is = TRUE, sep = "\t"
)
rownames(sperms)[rownames(sperms) %in% mature$PMID30404016] <-
  paste0(rownames(sperms)[rownames(sperms) %in% mature$PMID30404016], "*")

pdf(file = "spermatogenesis_genes_nmr.pdf", height = 16); print(spermatogenesis_genes); dev.off()

spermatogenesis_genes_rotate <- pheatmap(t(sperms), annotation_col = annotation_col, scale = "column",
                                         fontsize_row = 10, fontsize_col = 10, cellwidth = 10, cellheight = 10)
pdf(file = "spermatogenesis_genes_nmr_rotate.pdf", width = 16); print(spermatogenesis_genes_rotate); dev.off()

# Quick look at epithelial contamination marker genes (not saved to file)
pheatmap(rpkm_norm_stringent[rownames(rpkm_norm_stringent) %in%
                               c("Krt5", "Krt14", "Krt8", "Krt18", "Car12", "Cftr"), ])
df_stringent[df_stringent$Name %in% c("Krt5", "Krt14", "Krt8", "Krt18", "Car12", "Cftr"), ]


# ------------------------------------------------------------------------
# SECTION 12: Sperm-purity / epithelial-contamination proxy
# Uses Prm2 (protamine 2, mature sperm marker) vs. Krt5 (epithelial marker)
# as a proxy for sample purity, per the standard markers of good sperm
# preps (protamine ratio, transition protein depletion, absence of
# epithelial markers)
# ------------------------------------------------------------------------
genes <- rpkm_norm[c("mikado.chr23G692", "mikado.chr16G686"), ]
rownames(genes) <- c("Prm2", "Krt5")
colnames(genes) <- c("BaD7_sub", "CD34_sub", "963A_bre", "A32F_bre", "8521_bre")

df_stringent[df_stringent$Name %in% c("Crisp1", "Crisp2", "Crisp3", "Cd9", "Cd63"), ]
df[grep("Prm1", df$Name), ]
df[grep("Prm", df$Name), ]

markers <- pheatmap(log2(genes[c("Prm2", "Krt5"), ] + 1))
pdf(file = "Sperm_vs_epithelial_marker.pdf"); print(markers); dev.off()

# Z-scored contamination proxies, to be added as DE model covariates
PRM2_z <- scale(log2(genes["Prm2", ] + 1))
KRT5_z <- scale(log2(genes["Krt5", ] + 1))
rpkm_norm["mikado.chr23G691", ] # Prm1 (probably ncRNA copy)
rpkm_norm["mikado.chr23G102", ] # Prm1 segmental-duplication copy

scores$PRM2_z <- PRM2_z
scores$KRT5_z <- KRT5_z


# ------------------------------------------------------------------------
# SECTION 13: DE re-run controlling for sperm/epithelial contamination
# Same edgeR workflow as Section 5, but design now includes PRM2_z and
# KRT5_z as covariates alongside condition
# ------------------------------------------------------------------------
design <- model.matrix(~ condition + PRM2_z + KRT5_z, data = scores)
y <- DGEList(counts = nmr_counts_rpkm_gt1, group = scores$condition)
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
top <- topTags(lrt, n = nrow(nmr_counts_rpkm_gt1))$table

# Re-annotate as in Section 6
metrics <- read.table(
  "~/Desktop/NMR_multi_assembly/v090226/annotations/mHetGla1.pri_annotation_metrics_with_pseudogenes_v120426.tsv.gz",
  header = TRUE, as.is = TRUE, sep = "\t"
)
met1 <- metrics[order(metrics$score, decreasing = TRUE), ]
met_gene <- met1[!duplicated(met1$parent), ]
rownames(met_gene) <- met_gene$parent
met_gene_in <- met_gene[rownames(top), ]

top$Name <- met_gene_in$Name
top$Status <- met_gene_in$gene_type
table(top$FDR < 0.05 & abs(top$logFC > 2))
top$is_complete <- met_gene_in$is_complete
top$min_score <- met_gene_in$minScore

df <- top
df$significance <- "NS"
df$significance[df$FDR < 0.05 & df$logFC > 1] <- "bre"
df$significance[df$FDR < 0.05 & df$logFC < -1] <- "sub"
df$logFC[df$logFC < -30] <- -30  # cap extreme fold changes for plotting

admixture_volcano <- ggplot(df, aes(x = logFC, y = -log10(FDR), color = significance, shape = Status)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("bre" = "purple", "sub" = "orange", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_bw() +
  labs(x = "log2 Fold Change", y = expression(-log[10]("Adjusted p-value")), color = NULL)

pdf(file = "sperm_vs_bre_admixture_control.pdf"); print(admixture_volcano); dev.off()

# Repeat stringent-set filtering + labeling, as in Sections 8-9
df$Symbol <- !(grepl("mikado", df$Name) | grepl("ENSHGL", df$Name) | grepl("-copy", df$Name))
df_name <- df[!(grepl("mikado", df$Name) | grepl("ENSHGL", df$Name)), ]
df_noname <- df[(grepl("mikado", df$Name) | grepl("ENSHGL", df$Name)), ]
table(df_name$Status, df_name$significance)
table(df$Symbol, df$significance)
table(df$Status, df$significance)

df_stringent <- df[df$is_complete == "True" & df$Status == "gene" & df$Symbol == TRUE, ]
topgenes <- c(df_stringent$Name[df_stringent$logFC < -1][1:10],
              df_stringent$Name[df_stringent$logFC > 1][1:10])
df_stringent$topgenes <- ""
df_stringent$topgenes[df_stringent$Name %in% topgenes] <- df_stringent$Name[df_stringent$Name %in% topgenes]

library(ggrepel)
admixture_control_stringent <- ggplot(df_stringent, aes(x = logFC, y = -log10(FDR), color = significance, shape = Status)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("bre" = "purple", "sub" = "orange", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_classic() + theme(legend.position = "none") +
  labs(x = "log2 Fold Change", y = expression(-log[10]("Adjusted p-value")), color = NULL) +
  geom_text_repel(data = df_stringent, aes(label = topgenes), size = 3, max.overlaps = Inf)

library(gprofiler2)
pdf(file = "NMR_sperm_vs_bre_admixture_stringent.pdf"); print(admixture_control_stringent); dev.off()
ggsave("admixture_control_stringent.png", admixture_control_stringent, width = 4, height = 4)

SUB_DE <- gost(df_stringent$Name[df_stringent$significance == "sub"],
               organism = "mmusculus", custom_bg = df_stringent$Name, evcodes = TRUE)
BRE_DE <- gost(df_stringent$Name[df_stringent$significance == "bre"],
               organism = "mmusculus", custom_bg = df_stringent$Name, evcodes = TRUE)
BRE_DEGs <- df_stringent[df_stringent$significance == "bre", ]

admixture_sub <- ggplot(SUB_DE$result %>% head(5),
                        aes(x = reorder(term_name, -p_value), y = -log10(p_value), color = source)) +
  geom_col() + coord_flip() + theme_classic() +
  labs(x = "Pathway", y = "-log10(p-value)") +
  scale_color_manual(values = c("GO:BP" = "gold", "GO:CC" = "gold", "HP" = "gold", "REAC" = "gold"))

pdf(file = "NMR_sperm_vs_bre_admixture_sub_pathways.pdf", width = 8); print(admixture_sub); dev.off()


# ------------------------------------------------------------------------
# SECTION 14: Mouse sperm - internal QC and outlier detection
# 3 wild-type mouse replicates; checks for a within-species outlier sample
# ------------------------------------------------------------------------
mouseSperm_counts <- mouseSperm[, grep(".bam", colnames(mouseSperm))]
colnames(mouseSperm_counts) <- c("WT1", "WT2", "W3")
rownames(mouseSperm_counts) <- mouseSperm$Geneid

d <- DGEList(counts = mouseSperm_counts)
d$genes <- data.frame(Length = mouseSperm$Length)
rpkm_norm <- rpkm(d)
rpkm_norm <- rpkm_norm[rowMedians(rpkm_norm) > 1, ]
mouseSperm_counts_gt1 <- mouseSperm_counts[rownames(rpkm_norm), ]

mouse_corplot <- corrplot::corrplot(cor(rpkm_norm), addCoef.col = "white", number.cex = 3)
pdf(file = "mouse_withinspecies_corrplot.pdf")
corrplot::corrplot(cor(rpkm_norm), addCoef.col = "white", number.cex = 3)
dev.off()

# Arbitrary a/b/a grouping used purely to run a contrast that isolates
# whichever mouse sample is the outlier (WT2 vs WT1/WT3)
scores <- as.data.frame(c("a", "b", "a"))
colnames(scores) <- "condition"
rownames(scores) <- colnames(mouseSperm_counts)
design <- model.matrix(~condition, data = scores)

y <- DGEList(counts = mouseSperm_counts_gt1, group = scores$condition)
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
top <- topTags(lrt, n = nrow(nmr_counts_rpkm_gt1))$table

topgenes <- c(rownames(top)[top$logFC < -1][1:10])
top$significance <- "NS"
top$significance[top$FDR < 0.05 & top$logFC > 1] <- "WT3"
top$significance[top$FDR < 0.05 & top$logFC < -1] <- "WT1/WT5"
top$topgenes <- ""
top$topgenes[rownames(top) %in% topgenes] <- topgenes

mouse_outlier_MA <- ggplot(top, aes(x = logCPM, y = logFC, color = significance)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("WT3" = "purple", "WT1/WT5" = "orange", "NS" = "grey70")) +
  theme_classic() +
  labs(x = "logCPM", y = "logFC", color = NULL) +
  geom_text_repel(data = top, aes(label = topgenes), size = 3, max.overlaps = Inf)

ggsave("Mouse_outlier_MA.pdf", plot = mouse_outlier_MA, width = 8, height = 8)
ggsave("Mouse_outlier_MA_legend.png", plot = mouse_outlier_MA, width = 8, height = 8)
mouse_outlier_MA1 <- mouse_outlier_MA + theme(legend.position = "none")
ggsave("Mouse_outlier_MA_nolegend.png", plot = mouse_outlier_MA1, width = 8, height = 8)
pdf(file = "Mouse_outlier_MA_plot.pdf"); print(mouse_outlier_MA); dev.off()

# Pathway enrichment for genes driving the outlier signal
mouse_DE <- gost(rownames(top)[top$significance == "WT1/WT5"],
                 organism = "mmusculus", custom_bg = rownames(top), evcodes = TRUE)
mouse_DE$result <- mouse_DE$result[mouse_DE$result$source != "KEGG", ]

mouse_outlier <- ggplot(mouse_DE$result %>% head(5),
                        aes(x = reorder(term_name, -p_value), y = -log10(p_value), color = source)) +
  geom_col() + coord_flip() + theme_classic() + theme(legend.position = "none") +
  labs(x = "Pathway", y = "-log10(p-value)") +
  scale_color_manual(values = c("GO:BP" = "gold", "GO:CC" = "purple", "GO:MF" = "turquoise",
                                "HP" = "gold", "REAC" = "magenta"))

ggsave("mouse_outlier_GO.pdf", mouse_outlier, width = 4, height = 4)
ggsave("mouse_outlier_GO.png", mouse_outlier, width = 4, height = 4)
pdf(file = "mouse_outlier_gprofiler.pdf"); print(mouse_outlier); dev.off()


# ==============================================================================
# SECTION 15: Cross-species (Mouse vs. NMR) comparison
# ==============================================================================

# --- 15a. Load OrthoFinder results, identify strict 1:1 mouse-to-NMR orthologs ---
i <- 2
setwd(paste0("~/Desktop/NMR_multi_assembly/fig4_orthofinder/v070326_id_results/version2"))

orthogene_counts <- read.table(paste0("version", i, "_Orthogroups.GeneCount.tsv"),
                               header = TRUE, as.is = TRUE, sep = "\t")
print(nrow(orthogene_counts))

orthogene <- read.table(paste0("version", i, "_Orthogroups.tsv"),
                        header = TRUE, as.is = TRUE, sep = "\t")
rownames(orthogene) <- orthogene$Orthogroup
orthogene.df <- orthogene[, c("Mmus_proteins", "mHetGla.pri.mikadoID.version2")]

rownames(orthogene_counts) <- orthogene_counts[, 1]
orthogene_counts_species <- orthogene_counts[, c("Mmus_proteins", "mHetGla.pri.mikadoID.version2")]
colnames(orthogene.df) <- colnames(orthogene_counts_species) <- c("mouse", "nmr")

# Orthogroups where mouse has exactly 1 gene and NMR has >=1 (candidate 1:1 or 1:many)
mouse_eq1 <- orthogene_counts_species[orthogene_counts_species$mouse == 1 & orthogene_counts_species$nmr > 0, ]
orthogene_mmus1 <- orthogene.df[rownames(mouse_eq1), ]

# --- 15b. Collapse NMR multi-copy hits down to a single best gene per orthogroup ---
# Rationale: multiple NMR IDs per orthogroup are often pseudogene copies,
# segmental duplicates, or isoforms of the same underlying locus rather than
# true independent paralogs - this resolves them to one representative gene
metrics <- read.table(
  "~/Desktop/NMR_multi_assembly/v090226/annotations/mHetGla1.pri_annotation_metrics_with_pseudogenes_v120426.tsv.gz",
  header = TRUE, as.is = TRUE, sep = "\t"
)
setwd("/Users/dsokolowski/Desktop/NMR_multi_assembly/sperm_rnaseq")
dim(orthogene_mmus1)

other <- strsplit(orthogene_mmus1$nmr, ", ")
nmrID_out <- nmrID_out_lens <- c()
lens <- lengths(other)

for (i in 1:length(other)) {
  gene <- other[[i]]
  if (lens[i] == 1) {
    # Simple case: only one NMR gene in this orthogroup
    nmrID_out[i] <- gene
    nmrID_out_lens[i] <- 1
  } else {
    # Multiple NMR genes: look up their annotation metrics
    metrics_gene <- metrics[metrics$parent %in% gene, ]
    metrics_gene <- metrics_gene[!duplicated(metrics_gene$parent), ]
    if (!all(metrics_gene$gene_type == "gene")) {
      # Mixed gene/pseudogene/ncRNA set: pick the single highest-scoring entry
      gene <- metrics_gene$parent[which.max(metrics_gene$score)][1]
      nmrID_out[i] <- gene
      nmrID_out_lens[i] <- length(metrics_gene1)
    } else {
      # All entries are real genes: keep them only if there's exactly one,
      # otherwise collapse into a comma-joined multi-gene string
      metrics_gene1 <- metrics_gene[metrics_gene$gene_type == "gene", "parent"]
      if (length(metrics_gene1) > 1) {
        gene <- paste0(metrics_gene1, collapse = ",")
        nmrID_out[i] <- gene
        nmrID_out_lens[i] <- length(metrics_gene1)
      } else {
        gene <- metrics_gene1
        nmrID_out[i] <- gene
        nmrID_out_lens[i] <- 1
      }
    }
  }
}

orthogene_mmus1$nmr_short <- nmrID_out
orthogene_mmus1$nmr_short_lens <- nmrID_out_lens

# --- 15c. Restrict to strict 1:1 orthologs that are actually expressed ---
mouse_counts <- mouseSperm_counts
mmus_one_to_one <- orthogene_mmus1[orthogene_mmus1$nmr_short_lens == 1, ]
mmus_one_to_one_expressed <- mmus_one_to_one[
  mmus_one_to_one$mouse %in% rownames(mouseSperm_counts_gt1) |
    mmus_one_to_one$nmr_short %in% rownames(nmr_counts_rpkm_gt1),
]
mmus_one_to_one_expressed <- mmus_one_to_one_expressed[
  mmus_one_to_one_expressed$mouse %in% rownames(mouse_counts) &
    mmus_one_to_one_expressed$nmr_short %in% rownames(nmr_counts),
]

# --- 15d. Assemble combined count/RPKM matrices across species ---
NMR_counts_1to1 <- nmr_counts[mmus_one_to_one_expressed$nmr_short, ]
Mouse_counts_1to1 <- mouse_counts[mmus_one_to_one_expressed$mouse, ]

d <- DGEList(counts = mouseSperm_counts)
d$genes <- data.frame(Length = mouseSperm$Length)
mouse_rpkm <- rpkm(d)
mouse_rpkm_one2one <- mouse_rpkm[rownames(Mouse_counts_1to1), ]

d <- DGEList(counts = nmr_counts)
d$genes <- data.frame(Length = nmrSperm$Length)
nmr_rpkm <- rpkm(d)
NMR_rpkm_one2one <- nmr_rpkm[rownames(NMR_counts_1to1), ]

NMR_mmus_raw_counts <- as.data.frame(cbind(Mouse_counts_1to1, NMR_counts_1to1))
NMR_mmus_rpkm <- as.data.frame(cbind(mouse_rpkm_one2one, NMR_rpkm_one2one))
colnames(NMR_mmus_raw_counts) <- colnames(NMR_mmus_rpkm) <-
  c("WT1", "WT2", "WT3", "subB", "subC", "bre9", "breA", "bre8")

write.table(NMR_mmus_raw_counts, file = "mmus_mHetGla_one_to_one_orthologs.tsv",
            quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

# --- 15e. Cross-species correlation (sanity check on ortholog mapping quality) ---
pdf(file = "cross_species_spearman.pdf")
print(corrplot(cor(NMR_mmus_rpkm, method = "spearman"), addCoef.col = "orange", number.cex = 2))
dev.off()
pdf(file = "cross_species_spearman_nolabel.pdf")
print(corrplot(cor(NMR_mmus_rpkm, method = "spearman")))
dev.off()
png("cross_species_spearman_nolabel.png", res = 300, width = 4, height = 4)
corrplot(cor(NMR_mmus_rpkm, method = "spearman"))
dev.off()


# ------------------------------------------------------------------------
# SECTION 16: Cross-species DE - Mouse vs. NMR (all NMR samples, incl. exsubs)
# ------------------------------------------------------------------------
species <- as.data.frame(c(rep("mmus", 3), rep("nmr", 5)))
colnames(species) <- "species"
species$species <- factor(species$species, levels = c("nmr", "mmus"))
design <- model.matrix(~species, data = species)

y <- DGEList(counts = NMR_mmus_raw_counts, group = species$species)
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
top <- topTags(lrt, n = nrow(nmr_counts_rpkm_gt1))$table

df <- top
df$significance <- "NS"
df$significance[df$FDR < 0.05 & df$logFC > 2] <- "Mouse"
df$significance[df$FDR < 0.05 & df$logFC < -2] <- "NMR"
df$logFC[df$logFC < -30] <- -30

mouse_nmr_all_volcano <- ggplot(df, aes(x = logFC, y = -log10(FDR), color = significance)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Mouse" = "purple", "NMR" = "orange", "NS" = "grey70")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_classic() + theme(legend.position = "none") +
  labs(x = "log2 Fold Change", y = expression(-log[10]("Adjusted p-value")), color = NULL)

pdf(file = "Mouse_vs_NMR_allDE.pdf"); print(mouse_nmr_all_volcano); dev.off()
ggsave(filename = "Mouse_vs_NMR_allDEv0609.png", plot = mouse_nmr_all_volcano, height = 4, width = 4)

# Pathway enrichment for species-biased orthologs
library(gprofiler2)
df$Name <- rownames(df)
nmr_DE <- gost(df$Name[df$significance == "NMR"], organism = "mmusculus", custom_bg = df$Name, evcodes = TRUE)
mouse_DE <- gost(df$Name[df$significance == "Mouse"], organism = "mmusculus", custom_bg = df$Name, evcodes = TRUE)

nmrpg <- ggplot(nmr_DE$result %>% head(5),
                aes(x = reorder(term_name, -p_value), y = -log10(p_value), color = source)) +
  geom_col() + coord_flip() + theme_classic() + theme(legend.position = "none") +
  labs(x = "Pathway", y = "-log10(FDR)") +
  scale_color_manual(values = c("GO:BP" = "gold", "GO:CC" = "gold", "GO:MF" = "gold",
                                "HP" = "gold", "CORUM" = "gold"))

mouse_DE$result <- mouse_DE$result[mouse_DE$result$source == "GO:BP", ]
mmuspg <- ggplot(mouse_DE$result %>% head(5),
                 aes(x = reorder(term_name, -p_value), y = -log10(p_value), color = source)) +
  geom_col() + coord_flip() + theme_classic() + theme(legend.position = "none") +
  labs(x = "", y = "-log10(FDR)") +
  scale_color_manual(values = c("GO:BP" = "purple", "GO:CC" = "purple", "GO:MF" = "purple",
                                "HP" = "purple", "REAC" = "purple"))

pdf(file = "nmr_vs_mouse_allsample_gprofiler_v0619.pdf", width = 16); print(nmrpg + mmuspg); dev.off()
ggsave("nmr_vs_mouse_allsample_gprofiler.png", plot = nmrpg + mmuspg, width = 16, height = 4)


# ------------------------------------------------------------------------
# SECTION 17: Cross-validate edgeR results with a manual per-gene t-test
# Sanity check: does a simple log2(mean ratio) + t-test agree with edgeR's
# GLM-based fold change and FDR?
# ------------------------------------------------------------------------
tttest_out <- apply(NMR_mmus_rpkm, 1, function(x) {
  p <- t.test(log2(x[1:3] + 1), log2(x[4:8] + 1))
  if (mean(x[4:8]) == 0) { nmr <- 1e-3 } else { nmr <- mean(x[4:8]) }
  log2FC <- log2(mean(x[1:3]) / nmr)
  return(c(p$p.value, log2FC))
})
t_out <- as.data.frame(t(tttest_out))
colnames(t_out) <- c("P", "log2FCManual")
t_out$FDRManual <- p.adjust(t_out$P, method = "fdr")

mouse_vs_nmr_fin <- merge(df, t_out, by = "row.names")
cor(-log10(mouse_vs_nmr_fin$FDR), -log10(mouse_vs_nmr_fin$FDRManual), method = "spearman")

# Classify each gene by which method(s) called it significant
mouse_vs_nmr_fin$DEBoth <- "NS"
mouse_vs_nmr_fin$DEBoth[mouse_vs_nmr_fin$FDR < 0.05 & abs(mouse_vs_nmr_fin$logCPM) > 2 &
                          mouse_vs_nmr_fin$FDRManual < 0.05 & abs(mouse_vs_nmr_fin$log2FCManual) > 2] <- "both"
mouse_vs_nmr_fin$DEBoth[mouse_vs_nmr_fin$FDR > 0.05 & abs(mouse_vs_nmr_fin$logCPM) > 2 &
                          mouse_vs_nmr_fin$FDRManual < 0.05 & abs(mouse_vs_nmr_fin$log2FCManual) > 2] <- "ttest"
mouse_vs_nmr_fin$DEBoth[mouse_vs_nmr_fin$FDR < 0.05 & abs(mouse_vs_nmr_fin$logCPM) > 2 &
                          mouse_vs_nmr_fin$FDRManual > 0.05 & abs(mouse_vs_nmr_fin$log2FCManual) > 2] <- "edgeR"

edgeR_vs_T <- ggplot(mouse_vs_nmr_fin, aes(x = log2FCManual, y = logFC, color = DEBoth)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("NS" = "grey", "both" = "purple", "ttest" = "lightgreen", "edgeR" = "lightblue")) +
  theme_classic()

ggsave("FoldChanges_edgeR_vs_ttest.pdf", edgeR_vs_T)
ggsave("FoldChanges_edgeR_vs_ttest.png", edgeR_vs_T, width = 4, height = 4)

mouse_vs_nmr_fin_allsamples <- mouse_vs_nmr_fin
saveRDS(mouse_vs_nmr_fin_allsamples, file = "mouse_vs_nmr_DE_allsamples.rds")
saveRDS(NMR_mmus_raw_counts, file = "NMR_mmus_raw_counts.rds")
saveRDS(NMR_mmus_rpkm, file = "NMR_mmus_raw_counts.rds")


# ------------------------------------------------------------------------
# SECTION 18: Repeat cross-species DE with exsubs DROPPED
# Tests whether the mouse-vs-NMR signal is robust when restricted to
# only subordinate NMR samples (removing potential exsub-status confound)
# ------------------------------------------------------------------------
species <- as.data.frame(c(rep("mmus", 3), rep("nmr", 2)))
colnames(species) <- "species"
species$species <- factor(species$species, levels = c("nmr", "mmus"))
design <- model.matrix(~species, data = species)

y <- DGEList(counts = NMR_mmus_raw_counts[, 1:5], group = species$species)  # WT1-3 + subB/subC only
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
top <- topTags(lrt, n = nrow(nmr_counts_rpkm_gt1))$table

df <- top
df$significance <- "NS"
df$significance[df$FDR < 0.05 & df$logFC > 2] <- "Mouse"
df$significance[df$FDR < 0.05 & df$logFC < -2] <- "NMR"
df$logFC[df$logFC < -30] <- -30

volcano_nobre <- ggplot(df, aes(x = logFC, y = -log10(FDR), color = significance)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Mouse" = "purple", "NMR" = "orange", "NS" = "grey70")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_classic() +
  labs(x = "log2 Fold Change", y = expression(-log[10]("Adjusted p-value")), color = NULL)

library(gprofiler2)
df$Name <- rownames(df)
pdf(file = "NMR_vs_mouse_volcano_noBRE.pdf"); print(volcano_nobre); dev.off()
saveRDS(df, file = "NMR_vs_mouse_DE_noBRE.rds")

nmr_DE <- gost(df$Name[df$significance == "NMR"], organism = "mmusculus", custom_bg = df$Name, evcodes = TRUE)
mouse_DE <- gost(df$Name[df$significance == "Mouse"], organism = "mmusculus", custom_bg = df$Name, evcodes = TRUE)

gp1 <- ggplot(nmr_DE$result %>% head(20),
              aes(x = reorder(term_name, -p_value), y = -log10(p_value), color = source)) +
  geom_col() + coord_flip() + theme_bw() +
  labs(x = "Pathway", y = "-log10(p-value)") +
  scale_color_manual(values = c("GO:BP" = "red", "GO:CC" = "purple", "GO:MF" = "turquoise",
                                "HP" = "gold", "REAC" = "magenta"))

gp2 <- ggplot(mouse_DE$result %>% head(20),
              aes(x = reorder(term_name, -p_value), y = -log10(p_value), color = source)) +
  geom_col() + coord_flip() + theme_bw() +
  labs(x = "Pathway", y = "-log10(p-value)") +
  scale_color_manual(values = c("GO:BP" = "red", "GO:CC" = "purple", "GO:MF" = "turquoise",
                                "HP" = "gold", "REAC" = "magenta"))

pdf(file = "NMR_vs_mouse_sperm_pathway_enrichment_noBRE.pdf", width = 16); print(gp1 + gp2); dev.off()

# Manual t-test cross-check for the no-exsub comparison
tttest_out <- apply(NMR_mmus_rpkm, 1, function(x) {
  p <- t.test(log2(x[1:3] + 1), log2(x[4:5] + 1))
  if (mean(x[4:5]) == 0) { nmr <- 1e-3 } else { nmr <- mean(x[4:5]) }
  log2FC <- log2(mean(x[1:3]) / nmr)
  return(c(p$p.value, log2FC))
})
t_out <- as.data.frame(t(tttest_out))
colnames(t_out) <- c("P", "log2FCManual")
t_out$FDRManual <- p.adjust(t_out$P, method = "fdr")

mouse_vs_nmr_fin <- merge(df, t_out, by = "row.names")
cor(-log10(mouse_vs_nmr_fin$FDR), -log10(mouse_vs_nmr_fin$FDRManual), method = "spearman")

mouse_vs_nmr_fin$DEBoth <- "NS"
mouse_vs_nmr_fin$DEBoth[mouse_vs_nmr_fin$FDR < 0.05 & abs(mouse_vs_nmr_fin$logCPM) > 2 &
                          mouse_vs_nmr_fin$FDRManual < 0.05 & abs(mouse_vs_nmr_fin$log2FCManual) > 2] <- "both"
mouse_vs_nmr_fin$DEBoth[mouse_vs_nmr_fin$FDR > 0.05 & abs(mouse_vs_nmr_fin$logCPM) > 2 &
                          mouse_vs_nmr_fin$FDRManual < 0.05 & abs(mouse_vs_nmr_fin$log2FCManual) > 2] <- "ttest"
mouse_vs_nmr_fin$DEBoth[mouse_vs_nmr_fin$FDR < 0.05 & abs(mouse_vs_nmr_fin$logCPM) > 2 &
                          mouse_vs_nmr_fin$FDRManual > 0.05 & abs(mouse_vs_nmr_fin$log2FCManual) > 2] <- "edgeR"

noBRE_cor <- ggplot(mouse_vs_nmr_fin, aes(x = log2FCManual, y = logFC, color = DEBoth)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("NS" = "grey", "both" = "purple", "ttest" = "lightgreen", "edgeR" = "lightblue")) +
  theme_classic()
ggsave("FoldChanges_edgeR_vs_ttest_noBRE.pdf", noBRE_cor)
ggsave("FoldChanges_edgeR_vs_ttest_noBRE.png", noBRE_cor, width = 4, height = 4)

mouse_vs_nmr_fin_sub_only <- mouse_vs_nmr_fin
saveRDS(mouse_vs_nmr_fin_sub_only, file = "mouse_vs_nmr_DE_allsamples_noDE.rds")
dim(mouse_vs_nmr_fin_sub_only)
rownames(mouse_vs_nmr_fin_sub_only) == rownames(mouse_vs_nmr_fin_allsamples)


# ------------------------------------------------------------------------
# SECTION 19: Concordance check - with vs. without exsub samples
# Tests whether DE calls are consistent regardless of whether exsub
# samples are included, per gene (robustness check for Section 18 vs. 16-17)
# ------------------------------------------------------------------------
with_and_without_BRE <- as.data.frame(cbind(
  mouse_vs_nmr_fin_sub_only$Row.names, mouse_vs_nmr_fin_sub_only$logFC,
  mouse_vs_nmr_fin_sub_only$DEBoth,
  mouse_vs_nmr_fin_allsamples$logFC, mouse_vs_nmr_fin_allsamples$DEBoth
))
colnames(with_and_without_BRE) <- c("Gene", "logFCsub", "SigSub", "logFCall", "SigAll")
with_and_without_BRE$concordant <- with_and_without_BRE$SigSub == with_and_without_BRE$SigAll
with_and_without_BRE$logFCsub <- as.numeric(with_and_without_BRE$logFCsub)
with_and_without_BRE$logFCall <- as.numeric(with_and_without_BRE$logFCall)

noBRE_cor <- ggplot(with_and_without_BRE, aes(x = logFCall, y = logFCsub, color = SigAll, shape = concordant)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("NS" = "grey", "both" = "purple", "ttest" = "lightgreen", "edgeR" = "lightblue")) +
  theme_bw() + ggtitle("edgeR") + theme(legend.position = "none")

# Same concordance check using the manual t-test fold changes instead of edgeR's
with_and_without_BRE <- as.data.frame(cbind(
  mouse_vs_nmr_fin_sub_only$Row.names, mouse_vs_nmr_fin_sub_only$log2FCManual,
  mouse_vs_nmr_fin_sub_only$DEBoth,
  mouse_vs_nmr_fin_allsamples$log2FCManual, mouse_vs_nmr_fin_allsamples$DEBoth
))
colnames(with_and_without_BRE) <- c("Gene", "logFCsub", "SigSub", "logFCall", "SigAll")
with_and_without_BRE$concordant <- with_and_without_BRE$SigSub == with_and_without_BRE$SigAll
with_and_without_BRE$logFCsub <- as.numeric(with_and_without_BRE$logFCsub)
with_and_without_BRE$logFCall <- as.numeric(with_and_without_BRE$logFCall)

noBRE_cor_manual <- ggplot(with_and_without_BRE, aes(x = logFCall, y = logFCsub, color = SigAll, shape = concordant)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("NS" = "grey", "both" = "purple", "ttest" = "lightgreen", "edgeR" = "lightblue")) +
  theme_bw() + ggtitle("Manual")

pdf(file = "FC_distributions_include_or_exclude_BRE.pdf")
print(noBRE_cor + noBRE_cor_manual)
dev.off()

# Confusion-matrix-style heatmap: DE category with vs. without exsubs
pal <- colorRampPalette(c("white", "red"))
p1 <- pheatmap(table(mouse_vs_nmr_fin_sub_only$DEBoth, mouse_vs_nmr_fin_allsamples$DEBoth),
               color = pal(10), cluster_rows = FALSE, cluster_cols = FALSE,
               display_numbers = TRUE, number_format = "%.0f", fontsize_number = 12)
ggsave("Concordance_heatmap.pdf", p1, width = 4, height = 4)
ggsave("Concordance_heatmap.png", p1, width = 4, height = 4)
table(mouse_vs_nmr_fin_sub_only$DEBoth, mouse_vs_nmr_fin_allsamples$DEBoth)

mouse_vs_nmr_fin_allsamples <- readRDS("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure3/Figure3_KLM_mouse_vs_NMR/mouse_vs_nmr_DE_allsamples.rds")

# --- 20b. ActivePathways enrichment: are Mouse-biased DE genes enriched for
#          any RELAX category, relative to a background of all RELAX-tested genes? ---
apPlts <- list()
for (i in c("Mouse")) {
  MouseDE <- mouse_vs_nmr_fin_allsamples[mouse_vs_nmr_fin_allsamples$significance == i, ]
  head(as.data.frame(MouseDE$FDR, row.names = MouseDE$Row.names, col.names = "FDR"))
  
  mat <- as.matrix(as.data.frame(MouseDE$FDR, row.names = MouseDE$Row.names, col.names = "FDR"))
  colnames(mat) <- "FDR"
  
  rownames(mouse_vs_nmr_fin_allsamples) <- mouse_vs_nmr_fin_allsamples$Row.names
  bg <- rownames(mouse_vs_nmr_fin_allsamples)
  bg2 <- union(makeBackground(relaxed_gmt), bg)
  
  
  AP <- ActivePathways(scores = mat, gmt = "/Users/dsokolowski/Desktop/NMR_multi_assembly/sperm_proteomics/relaxed_sperm_genes_Sept10_26.gmt",
                       background = unique(c(makeBackground(relaxed_gmt), bg2)),
                       significant = 1)
  
  # Parse term IDs back into species/state components for plotting
  AP$species <- unlist(lapply(strsplit(AP$term_id, "_"), function(x) x[1]))
  AP$state <- do.call("rbind", str_split(AP$term_id, "_(?=[A-Z])", n = 2))[, 2]
  AP$state_ol <- paste0(lengths(AP$overlap), "/", AP$term_size)
  AP$species <- factor(AP$species, levels = c("Total", "allmolerat", "nmr"))
  AP <- AP[AP$species %in% c("Total", "allmolerat", "nmr"), ]
  
  apPlts[[i]] <- AP_plt <- ggplot(AP, aes(x = term_id, y = -log10(adjusted_p_val))) +
    geom_col(aes(fill = species), color = NA) +
    geom_col(aes(color = state), fill = NA, linewidth = 1) +
    geom_text(aes(y = 2, label = state_ol), vjust = 0.5, size = 6) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_manual(values = c("RELAXED" = "blue", "INTENSIFIED" = "red", "NOT_SIGNIFICANT" = "grey")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
    coord_flip()
  
  pdf(file = paste0(i, "_DE_relaxed_olap_v0619.pdf")); print(AP_plt); dev.off()
}

ggsave("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure3/Figure3_KLM_mouse_vs_NMR/Mouse_DE_relaxed_olap_fdr_20260910.pdf", AP_plt)
ggsave("relaxed_mouse_DEGs.png", apPlts$Mouse)
apPlts$Mouse + apPlts$NMR


# ------------------------------------------------------------------------
# SECTION 21: Genes that dropped out of the 1:1 ortholog / RELAX-tested set
# Diagnostic: are genes that failed QC/ortholog-matching systematically
# different (e.g., lower expression) from those that made it through?
# ------------------------------------------------------------------------
relax_analysis_out <- relax_analysis[!(rownames(relax_analysis) %in% rownames(NMR_mmus_rpkm)), ]
table(rownames(mouse_rpkm) %in% relax_analysis_out$Gene_mmus)

mouse_rpkm_relaxed_not_one2one <- mouse_rpkm[rownames(mouse_rpkm) %in% relax_analysis_out$Gene_mmus, ]
medians <- rowMedians(log2(mouse_rpkm + 1))
no_one_to_one <- rownames(mouse_rpkm) %in% relax_analysis_out$Gene_mmus

inter <- intersect(rownames(mouse_rpkm), relax_analysis_out$Gene_mmus)
mouse_rpkm_relaxed_not_one2one_test <- mouse_rpkm_relaxed_not_one2one[inter, ]
relax_analysis_out_test <- relax_analysis_out[inter, ]

df_relax <- as.data.frame(table(relax_analysis$Category_nmr_only))
df_relax_noNMR <- as.data.frame(table(relax_analysis_out_test$Category_nmr_only))
colnames(df_relax) <- colnames(df_relax_noNMR) <- c("Category", "Frequency")

noNMR <- ggplot(df_relax_noNMR, aes(x = Category, y = Frequency)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_cartesian(ylim = c(2, NA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
  coord_flip() + theme_bw() + ggtitle("No NMR RNASeq")

allRelax <- ggplot(df_relax, aes(x = Category, y = Frequency)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_cartesian(ylim = c(2, NA)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) +
  coord_flip() + theme_bw() + ggtitle("All Spermatogenesis genes")

noNMR + allRelax
barplot(table(relax_analysis_out_test$Category_nmr_only), las = 2, horiz = TRUE)
pheatmap(table(relax_analysis_out_test$Category_nmr_only, relax_analysis_out_test$Category_bathyergidae),
         display_numbers = TRUE, number_format = "%.0f", fontsize_number = 12)
boxplot(medians ~ no_one_to_one)

# Identify NMR expression levels for genes that dropped out of the 1:1 set
nmr_rpkm_relaxed_not_one2one <- nmr_rpkm[rownames(nmr_rpkm) %in% relax_analysis_out$Gene_mmus, ]

metrics <- read.table(
  "~/Desktop/NMR_multi_assembly/v090226/annotations/mHetGla1.pri_annotation_metrics_with_pseudogenes_v120426.tsv.gz",
  header = TRUE, as.is = TRUE, sep = "\t"
)
metrics_relax_analysis_out <- c()
for (i in relax_analysis_out$Gene_mmus) {
  metrics_relax_analysis_out <- rbind(metrics_relax_analysis_out,
                                      metrics[metrics$Name == relax_analysis_out$Gene_mmus, ])
}

# Manually curated list of genes known from the literature to be "lost"/pseudogenized in NMR
lostgene <- c("Fsip2", "Zan", "Fam71b", "Rimbp3", "Ttll8", "Fsip2", "Cabs1", "Tmprss12",
              "Atp1a4", "Dpy19l2", "Ndc80", "Adig", "Irgc")
df[df$Name %in% lostgene, ]


# ==============================================================================
# SECTION 22: Non-coding RNA expression burden - Mouse vs. NMR sperm
# Tests whether the fraction of expression coming from non-protein-coding
# genes (lncRNA/ncRNA_gene) differs between species and between NMR
# subordinate/exsub status
# ==============================================================================
library(rtracklayer)

nmrSperm <- read.table("NMRsperm.txt", header = TRUE, as.is = TRUE)
nmr_counts <- nmrSperm[, grep("sorted.bam", colnames(nmrSperm))]
rownames(nmr_counts) <- nmrSperm$Geneid

mouse <- readGFF("C57_T2T.gff")
mouse <- mouse[mouse$type == "gene", ]
mouse_gene <- mouse[, c("gene", "gene_biotype")]

nmr_counts <- nmrSperm[, grep("sorted.bam", colnames(nmrSperm))]
rownames(nmr_counts) <- nmrSperm$Geneid

dge_mouse <- DGEList(counts = mouseSperm_counts)
mouse_cpm <- cpm(dge_mouse)
dge_nmr <- DGEList(counts = nmr_counts)
nmr_cpm <- cpm(dge_nmr)

head(metrics)
metrics_gene <- metrics[!duplicated(metrics$parent), ]
rownames(metrics_gene) <- metrics_gene$parent
metrics_gene <- metrics_gene[rownames(nmr_cpm), ]

mouse_gene <- mouse_gene[!duplicated(mouse_gene$gene), ]
rownames(mouse_gene) <- mouse_gene$gene
mouse_gene <- mouse_gene[rownames(mouse_cpm), ]

# Fraction of total CPM coming from non-protein-coding genes, per sample
mmus_nonPCG <- 1 - (colSums(mouse_cpm[mouse_gene$gene_biotype == "protein_coding", ]) /
                      colSums(mouse_cpm[mouse_gene$gene_biotype %in% c("protein_coding", "lncRNA"), ]))
nmr_nonPCG <- 1 - (colSums(nmr_cpm[metrics_gene$gene_type == "gene", ]) /
                     colSums(nmr_cpm[metrics_gene$gene_type %in% c("gene", "ncRNA_gene"), ]))

ncRNA_burden <- as.data.frame(c(mmus_nonPCG, nmr_nonPCG))
ncRNA_burden$state <- factor(c("mouse", "mouse", "mouse", "sub", "sub", "bre", "bre", "bre"),
                             levels = c("sub", "bre", "mouse"))
colnames(ncRNA_burden) <- c("ncRNA", "state")

# One-way ANOVA + Tukey HSD post-hoc: does non-coding burden differ by group?
TukeyHSD(aov(ncRNA_burden$ncRNA ~ ncRNA_burden$state))
boxplot(ncRNA_burden$ncRNA ~ ncRNA_burden$state, xlab = "", ylab = "ncRNA expression fraction", ylim = c(0, .1))