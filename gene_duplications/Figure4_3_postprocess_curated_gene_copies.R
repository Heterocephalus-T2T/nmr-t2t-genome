library(Biostrings)
library(rtracklayer)
library(ggplot2)

# ===========================================================================
#
# Downstream gene-copy-number / segmental-duplication (SD) analysis for
# NMR (mHetGla), building on OrthoFinder orthogroup output. Includes:
#   - extracting orthogroup member genes/proteins
#   - filtering a candidate intra-species SD gene pair list (removing
#     genes with truncated CDS, chr15p region, and known multi-copy
#     gene family artifacts like ORs, ZNFs, histones, VNRs)
#   - cross-referencing against manually curated gene-copy tables for
#     NMR, AMR (Fukomys), and BatSui
#   - visualizing final curated copy-number sets as heatmaps, with MGI
#     mouse phenotype annotation
#
# NOTE - this script has all of the intermediary steps required to convert
# processed segmental duplications into gene copy events. To regenerate the 
# figures in fig4 of Sokolowski et al., 2026, running 
# Figure4_4_B_C_Pms1_structure_and_expression.R is sufficient.
#
# ===========================================================================

# ---------------------------------------------------------------------
# Helper: given an orthogroup ID (row name in `ortho`), extract and
# clean the list of mHetGla pri version2 gene IDs belonging to it
# ---------------------------------------------------------------------
split_ogroup <- function(x) {
  genes <- ortho[x, "mHetGla.pri.version2"]
  y <- unlist(strsplit(genes, ","))
  y <- gsub(" ", "", y)
  return(y)
}


# Quick look at feature type counts in the MT/nuclear annotation GFF
gff <- readGFF("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/mHetGla1.pri.MT.v020326.gff.gz")


# Load the "version2" (intact, one-per-gene) NMR pri protein set and
# gene symbol metrics table
pri <- readAAStringSet("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure4/Figure4A_Duplications/gene_copy_orthofinder/referece_proteins/mHetGla.pri.version2.faa")
geneSyms <- read.table("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/mHetGla1.1.pri_gene_symbols_metrics.tsv",
                       header = TRUE, as.is = TRUE, sep = "\t")



# ---------------------------------------------------------------------
# Restrict the full GFF to genes present in the "version2" protein set
# and export as a simplified GFF3, then summarize StringTie-supported
# mRNA transcripts by tissue of origin (same "tissue from alias" logic
# used in earlier scripts in this conversation)
# ---------------------------------------------------------------------
gff_simple <- gff[gff$Name %in% names(pri), ]
export(gff_simple, "mHetGla1.pri.MT.v020326.version2.simple.gff", "GFF3")

gff_stringtie <- gff_simple[gff_simple$type == "mRNA" & grepl("stringtie", gff_simple$alias), ]
gff_stringtie$tissue <- unlist(lapply(strsplit(gff_stringtie$alias, "\\."), function(x) {
  return(x[1])
}))
gff_stringtie$tissue <- gsub("stringtie_", "", gff_stringtie$tissue)
data <- as.data.frame(table(gff_stringtie$tissue))


# ---------------------------------------------------------------------
# Load a table of intra-species segmental duplication (SD) candidate
# gene pairs (identified per orthogroup) and the corresponding
# OrthoFinder "version2" orthogroup table
# ---------------------------------------------------------------------

setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure4/Figure4A_Duplications")
intraSDs <- read.table("gene_copy_analysis/orthogroup_intraSDs.tsv", header = T, as.is = T, sep = "\t")
ortho <- read.table("gene_copy_orthofinder/Orthofinder/version2_Orthogroups.tsv", header = T, as.is = T, sep = "\t")
rownames(ortho) <- ortho$Orthogroup

# ---------------------------------------------------------------------
# Pri - make it smaller by only keeping genes in OG where the CDS is
# >40% the length of the max (comment says 40%, code below actually
# uses a >0.5 (50%) cutoff - worth reconciling which threshold was
# intended)
#
# For each orthogroup, normalize each member gene's protein length by
# the longest member's length, and keep only genes retaining more than
# half that maximum length (a heuristic to drop severely truncated/
# fragmentary gene models before downstream SD analysis)
# ---------------------------------------------------------------------

widths <- lengths(pri)
names(widths) <- names(pri)
pri_small <- c()
for (i in 1:nrow(ortho)) {
  OG <- ortho$Orthogroup[i]
  genes <- split_ogroup(OG)
  
  width1 <- widths[genes] / max(widths[genes])
  width2 <- width1[width1 > 0.5]
  pri_small <- c(pri_small, names(width2))
}
pri_small_aa <- pri[pri_small]

# ---------------------------------------------------------------------
# Part 1
## Clean Orthogroup based on existence in version2 set
#
# For each candidate SD gene pair, keep it only if BOTH member genes
# survived the length-filtering step above (pri_small_aa). Also record,
# per row, how many total NMR copies of that orthogroup remain post-
# filtering (used for downstream prioritization/sorting).
# ---------------------------------------------------------------------
uni <- sort(unique(intraSDs$orthogroup))
cleanedOG <- rep(TRUE, nrow(intraSDs))
for (i in 1:nrow(intraSDs)) {
  OG <- intraSDs$orthogroup[i]
  genes <- split_ogroup(OG)
  if (!((intraSDs$gene1[i] %in% names(pri_small_aa) & intraSDs$gene2[i] %in% names(pri_small_aa)))) {
    cleanedOG[i] <- FALSE
  }
  
  genes1 <- intersect(genes, names(pri_small_aa))
  intraSDs$nmr_copies_in_og_total[i] <- length(genes1)
}
intraSDsV2 <- intraSDs[cleanedOG, ]

# Filter for chr15p
# Remove SD pairs anchored on chr15 where either partner falls below
# 26 Mb (i.e., within the chr15p "supergene"/segmental-duplication
# region already characterized elsewhere in this analysis - excluded
# here to avoid re-flagging genes already accounted for in that region)
intraSDsV3 <- intraSDsV2[!(intraSDsV2$chrom == "chr15" & (intraSDsV2$end1 < 2.6e7 | intraSDsV2$end2 < 2.6e7)), ]

# Filter for gene families
# Remove SD pairs belonging to well-known large multi-copy gene
# families that are expected to show apparent "duplication" for
# reasons unrelated to genuine lineage-specific expansion (vomeronasal/
# olfactory receptors, zinc fingers, keratins, histones)
intraSDsV4 <- intraSDsV3[!(grepl("VN1R", intraSDsV3$human_symbols) | grepl("OR", substr(intraSDsV3$human_symbols, 1, 2)) | grepl("ZNF", substr(intraSDsV3$human_symbols, 1, 3)) | grepl("KRT", substr(intraSDsV3$human_symbols, 1, 3)) | grepl("VMN1", substr(intraSDsV3$human_symbols, 1, 4)) | grepl("H2A", substr(intraSDsV3$human_symbols, 1, 4)) | grepl("H3-", substr(intraSDsV3$human_symbols, 1, 4)) | grepl("H4C", substr(intraSDsV3$human_symbols, 1, 4))), ]

write.table(intraSDsV4, file = "gene_copy_analysis/clean_orthogroup_intraSDs.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")



# ---------------------------------------------------------------------
# Load the manually curated version of the cleaned intra-SD gene-pair
# table with "Designation" labels. Deduplicate to one row per orthogroup,
# then plot the interpretability x Designation breakdown as a heatmap.
# ---------------------------------------------------------------------
gened <- read.delim("DS - clean_orthogroup_intraSDs.txt", header = T)
gene_nodup <- gened[!duplicated(gened$orthogroup), ]
OGs <- as.matrix(table(gene_nodup$interpretability, gene_nodup$Designation))

# ---------------------------------------------------------------------
# Extract the curated "real" (Designation == "R") duplicated genes,
# sort by total NMR copy count, deduplicate by gene symbol, and save
# both the full curated "R" table and a flat list of gene symbols
# (some of which are semicolon-delimited multi-gene entries) for
# downstream tandem-duplication analysis / GO enrichment
# ---------------------------------------------------------------------
R_genes <- unique(toupper(gened$human_symbols[gened$Designation == "R"]))
gene_nodup <- gene_nodup[gene_nodup$Designation == "R", ]
gene_nodup <- gene_nodup[order(gene_nodup$nmr_copies_in_og_total, decreasing = T), ]
gene_nodup$human_symbols <- toupper(gene_nodup$human_symbols)
gene_nodup <- gene_nodup[!duplicated(gene_nodup$human_symbols), ]

write.table(gene_nodup, file = "gene_copy_analysis/clean_orthogroup_intraSDs_Rgene.txt",
            quote = F, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(unlist(strsplit(R_genes, ";")), file = "gene_copy_analysis/Tandem_Dup.txt", quote = F, row.names = F, col.names = F, sep = "\t")

##
##### One-to-one analysis
##
# ---------------------------------------------------------------------
# Cross-check single-copy ("one-to-one" orthologous) gene support
# against the length-filtered pri_small_aa set, then load a manually
# curated version of this table for downstream integration with the
# "R" (real duplication) gene set above
# ---------------------------------------------------------------------
single_copy <- read.table("gene_copy_analysis/Orthogene_counts_singe_copy_version2.gene_support.tsv",
                          header = TRUE, as.is = TRUE, sep = "\t")
new_copy <- new_copy_count <- c()
for (i in 1:nrow(single_copy)) {
  gene <- strsplit(single_copy$focal_copy_ids_kept[i], ",")[[1]]
  gene1 <- gene[gene %in% names(pri_small_aa)]
  new_copy[i] <- paste0(gene, collapse = ",")
  new_copy_count[i] <- length(gene)
}

curated_one2one <- read.delim("gene_copy_analysis/ds-Orthogene_counts_singe_copy_version2.gene_support.txt",
                              header = TRUE, as.is = TRUE, sep = "\t")
table(curated_one2one$Designation)
barplot(table(curated_one2one$Designation), las = 2, horiz = T)

# Standardize Designation labels (uppercase, strip spaces) and extract
# the "R" (real duplication) subset of the one-to-one curated table
curated_one2one$Designation <- toupper(gsub(" ", "", curated_one2one$Designation))
curated_gene <- curated_one2one$gene.1[curated_one2one$Designation == "R"] 
curated_one2one_R <- curated_one2one[curated_one2one$Designation == "R", ]
curated_one2one_R <- curated_one2one_R[!is.na(curated_one2one_R$gene), ] 

# ---------------------------------------------------------------------
# Combine curated one-to-one "R" genes with the orthogroup-level "R"
# genes (R_genes) and manually added extras into one master duplicated-
# gene list, run enrichment, and save
# ---------------------------------------------------------------------
inter1 <- intersect(toupper(curated_one2one_R$gene), unlist(strsplit(R_genes, ";")))
manually_added <- read.table("~/Desktop/NMR_multi_assembly/fig4_orthofinder/extra_gene_copies_v290326.txt",
                             header = T, as.is = T, sep = "\t")
interchromosomal <- unlist(lapply(strsplit(R_genes, ";"), function(x) return(x[1])))

all <- toupper(unique(c(toupper(curated_one2one_R$gene), R_genes, manually_added$Gene)))
# Manual correction: an apparent misannotated/miscalled gene symbol swap
all[all == "TUBA1A"] <- "TUBA3B"

gost_dups <- gost(all, evcodes = T)
write.table(all, file = "gene_copy_analysis/NMR_gene_copies_v230326.txt",
            quote = F, row.names = F, col.names = F, sep = "\t")

# ---------------------------------------------------------------------
# Reload the manually added extras table
# ---------------------------------------------------------------------
manually_added <- read.table("gene_copy_analysis/extra_gene_copies_v290326.txt",
                             header = T, as.is = T, sep = "\t")


# ---------------------------------------------------------------------
# Cross-species gene copy number comparison: load AMR (Fukomys, listed
# as "bathyergid") copy number data and a curated NMR copy number set,
# clean AMR gene names (strip "-COPY" suffixes), and identify candidate
# genes present in both
# ---------------------------------------------------------------------
AMR_copy <- read.table("gene_copy_analysis/bathyergid_gene_copies_v6_v7_comparison.tsv", header = T, as.is = T, sep = "\t")
curated_copy <- read.delim("gene_copy_analysis/curated_NMR_gene_copy_number_set_v290326.txt", header = T, as.is = T, sep = "\t")
AMR_copy$GENE_NAME <- unlist(lapply(strsplit(toupper(AMR_copy$Gene_Name), "-COPY"), function(x) return(x[1])))
AMR_copy_candidate <- AMR_copy[AMR_copy$GENE_NAME %in% curated_copy$Gene, ]



# ---------------------------------------------------------------------
# African mole-rat (Fukomys/FukAns) BISER segmental duplication (SD)
# calls: quantify total intrachromosomal vs interchromosomal SD bases
# ---------------------------------------------------------------------
bedpe <- read.table("gene_copy_analysis/BiserFukAns.bedpe", header = F, as.is = T, sep = "\t")

intrachr <- bedpe[bedpe$V1 == bedpe$V4, ]
sum(intrachr$V3 - intrachr$V2) + sum(intrachr$V6 - intrachr$V5)

interchr <- bedpe[bedpe$V1 != bedpe$V4, ]
(sum(interchr$V3 - interchr$V2) + sum(interchr$V6 - interchr$V5)) / 1e6

# ---------------------------------------------------------------------
# Clean BatSui GFF for extra layer
# Restrict the BatSui gene-symbol-annotated GFF to genes present in the
# BatSui "version2" (intact) protein set, and export a simplified GFF3
# ---------------------------------------------------------------------
library(rtracklayer)
sym <- readGFF("gene_copy_analysis/BatSui.edited.geneSymbol.gff.gz")
batsuiV2 <- readAAStringSet("gene_copy_orthofinder/referece_proteins/BatSui.version2.faa")
symID <- sym[sym$Name %in% names(batsuiV2), ]
export(object = symID, con = "BatSui.version2.simple.gff", format = "GFF3")

##
##### Visualize final set
##
# ---------------------------------------------------------------------
# Load the final curated NMR gene copy number set (across species: NMR/
# HetGla, AMR/FukAns, BatSui) and plot as a heatmap of copy counts
# ---------------------------------------------------------------------
curated <- read.delim("gene_copy_analysis/curated_NMR_gene_copy_number_set_v010426.txt",
                      header = T, as.is = T, sep = "\t")
mats <- curated[, grep("Copies", names(curated))]
mats  # (autoprint - quick inspection)
rownames(mats) <- curated$Gene
mats  # (autoprint again, post-rownames assignment)

pdf("all_genes.pdf", height = 11)
print(pheatmap(mats, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_row = 6))
dev.off()

# ---------------------------------------------------------------------
# Flag genes located on pseudoautosomal/chrY/autosomal-ambiguous
# regions in ANY of the three species' notes/copy fields, to exclude
# them from the "different copy number" comparison below (since sex-
# chromosome-linked or pseudoautosomal genes can show spurious apparent
# copy-number differences due to how each assembly represents sex
# chromosomes, rather than genuine lineage-specific duplication)
# ---------------------------------------------------------------------
gene <- rep(TRUE, nrow(mats))
for (i in 1:nrow(mats)) {
  if (grepl("chrY", curated$HetGla.Note[i])) {
    gene[i] <- FALSE
  }
  if (grepl("auto", curated$HetGla.Note[i])) {
    gene[i] <- FALSE
  }
  ###
  if (grepl("chrY", curated$FukAns.Copies[i])) {
    gene[i] <- FALSE
  }
  if (grepl("auto", curated$FukAns.Copies[i])) {
    gene[i] <- FALSE
  }
  ###
  if (grepl("chrY", curated$BatSui.Copies[i])) {
    gene[i] <- FALSE
  }
  if (grepl("auto", curated$BatSui.Copies[i])) {
    gene[i] <- FALSE
  }
}
curated_noPseudo <- curated[gene, ]
write.table(curated_noPseudo, file = "gene_copy_analysis/curated_NMR_gene_copy_number_set_noPseudoautosomal_v010426.txt",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

mats_noPseudo <- mats[gene, ]
pdf("no_pseudoautosomal.pdf", height = 11)
print(pheatmap(mats_noPseudo, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_row = 6))
dev.off()

# ---------------------------------------------------------------------
# Another curated gene-copy-events spreadsheet, with binary trait
# annotations (Yes/blank across several biological categories). Build
# a copy-number heatmap with a matching binary annotation sidebar.
# ---------------------------------------------------------------------
geneCopy <- read.delim("gene_copy_analysis/v040426 NNMR gene copy events - Sheet1.tsv",
                       header = T, as.is = T, sep = "\t")
geneCopy <- geneCopy[geneCopy$Gene != "", ]
geneCopy <- geneCopy[!(colnames(geneCopy) %in% c("Behaviour", "Vascular"))]

# Keep only genes flagged "Yes" in at least one trait column
gene_copy_drop <- geneCopy[apply(geneCopy, 1, function(x) {
  return(any(x == "Yes"))
}), ]

annotation_row <- gene_copy_drop[, colnames(geneCopy)[4:10]]
rownames(annotation_row) <- gene_copy_drop$Gene

# Convert Yes/blank trait flags to 1/0 for pheatmap annotation coloring
for (i in colnames(annotation_row)) {
  annotation_row[, i][annotation_row[, i] == "Yes"] <- 1
  annotation_row[, i][annotation_row[, i] != 1] <- 0
}

gene_copy_plot <- gene_copy_drop[, grep("Copies", colnames(gene_copy_drop))]
rownames(gene_copy_plot) <- gene_copy_drop$Gene
white_to_red <- colorRampPalette(c("white", "red"))(10)

# Order both the copy-number matrix and its annotation rows identically
# (alphabetically by gene) so pheatmap rows line up correctly
gene_copy_plot <- gene_copy_plot[order(rownames(gene_copy_plot)), ]
annotation_row <- annotation_row[order(rownames(annotation_row)), ]

ann_colors <- list()
for (i in colnames(annotation_row)) {
  ann_colors[[i]] <- c("0" = "white", "1" = "black")
}

pheatmap(gene_copy_plot, color = white_to_red, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = T, number_format = "%.0f",
         cellwidth = 24, annotation_row = annotation_row, annotation_colors = ann_colors)

# ---------------------------------------------------------------------
# PMS1-specific gene copy matrix visualization (binary presence/absence
# heatmap, black/white)
# ---------------------------------------------------------------------
pms1 <- read.csv("gene_copy_analysis/pms1_matrix.csv", header = T, as.is = T)
rownames(pms1) <- pms1$gene.symbol

pms2 <- pms1 <- pms1[, -1]
white_to_red <- colorRampPalette(c("white", "black"))(2)

pdf(file = "gene_copy_analysis/Figure4C_Pms1_transcript_expression.pdf")
pheatmap(t(pms2), color = white_to_red, cluster_rows = FALSE, cluster_cols = FALSE, number_format = "%.0f",
         cellwidth = 12, cellheight = 12)
dev.off()
