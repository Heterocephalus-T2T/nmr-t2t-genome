library("gprofiler2")
# setwd("~/Desktop/NMR_multi_assembly/summarize_NMR_assemblies/")
# ===========================================================================
# Compare TOGA-inferred gene loss/intactness calls across NMR and other
# species (BatSui, CDMR, etc.) to identify genes lost specifically in
# NMR (naked mole-rat), test for GO/Reactome enrichment among lost gene
# sets (sperm/spermatogenesis-related pathways), then cross-reference
# candidate "lost" sperm genes against gene annotation quality metrics
# to flag possible segmental duplication candidates and testis-specific
# genes.
# ===========================================================================

# ---------------------------------------------------------------------
# Load TOGA gene-loss calls for BatSui, keep only gene-level rows, and
# strip the "gene-" prefix from gene IDs for consistency with other files
# ---------------------------------------------------------------------
BatSui <- read.table("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure3/Figure3C_GeneLoss/toga_loss/BatSui.loss_summ_data.tsv")
BatSui <- BatSui[BatSui$V1 == "GENE", ]
BatSui$V2 <- gsub("gene-", "", BatSui$V2)

# ---------------------------------------------------------------------
# Load TOGA gene-loss summary files for all other species/assemblies in
# the current directory (each file corresponds to one species/lineage),
# apply the same gene-level filtering and ID cleanup, and store each in
# a named list keyed by species (derived from the filename)
# ---------------------------------------------------------------------
lf <- list.files(pattern = "mHetGla1.pri.loss_summ_data.tsv",path = "~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure3/Figure3C_GeneLoss/toga_loss/")
toga_losses <- list()
for (i in lf) {
  togas <- read.table("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure3/Figure3C_GeneLoss/toga_loss/mHetGla1.pri.loss_summ_data.tsv")
  togas <- togas[togas$V1 == "GENE", ]
  togas$V2 <- gsub("gene-", "", togas$V2)
  toga_losses[[i]] <- togas
}
names(toga_losses) <- gsub("_pri_loss_sum_data.tsv", "", names(toga_losses))

# Add BatSui into the same list so it's processed alongside the other species
toga_losses$BatSui <- BatSui

# ---------------------------------------------------------------------
# Build the full union of gene IDs seen across all species' TOGA results
# ---------------------------------------------------------------------
unique_genes <- c()
for (i in names(toga_losses)) {
  unique_genes <- unique(c(unique_genes, toga_losses[[i]]$V2))
}

# ---------------------------------------------------------------------
# Build a gene x species matrix of TOGA loss/intactness calls
# (e.g., "L" = lost, "I" = intact, "UL" = uncertain loss, etc.),
# defaulting to "M" (missing/not tested) where a species has no call
# for a given gene
# ---------------------------------------------------------------------
mat <- matrix("M", nrow = length(unique_genes), ncol = length(toga_losses))
rownames(mat) <- unique_genes
colnames(mat) <- names(toga_losses)

for (i in names(toga_losses)) {
  toga1 <- toga_losses[[i]]
  rownames(toga1) <- toga1$V2
  mat[rownames(toga1), i] <- toga1$V3
}

# Save the combined gene-call matrix for reuse without re-running the
# TOGA file parsing above
write.table(mat, file = "NMR_BatSui_TOGA_gene_assignments_040126.tsv", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

# Reload from disk (as a fresh starting point for the analyses below)
mat <- read.table("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/SummarizeAssemblies/NMR_BatSui_TOGA_gene_assignments_040126.tsv", header = T, as.is = T, sep = "\t")

# NMR and BatSui full loss

# ---------------------------------------------------------------------
# GO/Reactome enrichment on genes lost ("L") across ALL species/columns
# in the matrix
# ---------------------------------------------------------------------
## All Loss
all_loss <- apply(mat, 1, function(x) all(x == "L"))
all_loss_genes <- names(all_loss)[all_loss == TRUE]
all_loss_gprof <- gprofiler2::gost(all_loss_genes, organism = "mmusculus", custom_bg = unique_genes, correction_method = "fdr", evcodes = TRUE)
all_loss_gprof_out <- all_loss_gprof$result[grepl("GO:", all_loss_gprof$result$source) | grepl("REAC", all_loss_gprof$result$source), ]
all_loss_gprof_out <- all_loss_gprof_out[order(all_loss_gprof_out$p_value), ]

# ---------------------------------------------------------------------
# Same enrichment, but relaxing "loss" to include "L" or "UL" (uncertain
# loss) calls across all species
# ---------------------------------------------------------------------
# All loss or UL
all_loss_or_UL <- apply(mat, 1, function(x) all(x == "L" | x == "UL"))
all_loss_or_UL_genes <- names(all_loss_or_UL)[all_loss_or_UL == TRUE]
all_loss_or_UL_gprof <- gprofiler2::gost(all_loss_or_UL_genes, organism = "mmusculus", custom_bg = unique_genes, correction_method = "fdr", evcodes = TRUE)
all_loss_or_UL_gprof_out <- all_loss_or_UL_gprof$result[grepl("GO:", all_loss_or_UL_gprof$result$source) | grepl("REAC", all_loss_or_UL_gprof$result$source), ]
all_loss_or_UL_gprof_out <- all_loss_or_UL_gprof_out[order(all_loss_or_UL_gprof_out$p_value), ]

# ---------------------------------------------------------------------
# Genes lost ("L") in columns 1-4 and intact ("I") in column 5.
# Enrichment run on this set.
# ---------------------------------------------------------------------
# NMR loss
NMR_loss <- apply(mat, 1, function(x) (all(x[1:4] == "L") & x[5] == "I"))
NMR_loss_genes <- names(NMR_loss)[NMR_loss == TRUE]
NMR_loss_genes_gprof <- gprofiler2::gost(NMR_loss_genes, organism = "mmusculus", custom_bg = unique_genes, correction_method = "fdr", evcodes = TRUE)
NMR_loss_genes_gprof_out <- NMR_loss_genes_gprof$result[grepl("GO:", NMR_loss_genes_gprof$result$source) | grepl("REAC", NMR_loss_genes_gprof$result$source), ]
NMR_loss_genes_gprof_out <- NMR_loss_genes_gprof_out[order(NMR_loss_genes_gprof_out$p_value), ]

# Manual inspection of enrichment results, rows 16-30 (by p-value rank)
View(NMR_loss_genes_gprof_out[16:30, c("term_name", "p_value", "intersection_size", "intersection"), ])

# Same "NMR-specific loss" logic, but relaxed to include UL calls in the
# comparison species
NMR_loss_or_UL <- apply(mat, 1, function(x) (all(x[1:4] == "L" | x[1:4] == "UL") & (x[5] == "I")))
NMR_loss_or_UL_genes <- names(NMR_loss_or_UL)[NMR_loss_or_UL == TRUE]
NMR_loss_or_UL_gprof <- gprofiler2::gost(NMR_loss_or_UL_genes, organism = "mmusculus", custom_bg = unique_genes, correction_method = "fdr", evcodes = TRUE)
NMR_loss_or_UL_gprof_out <- NMR_loss_or_UL_gprof$result[grepl("GO:", NMR_loss_or_UL_gprof$result$source) | grepl("REAC", NMR_loss_or_UL_gprof$result$source), ]
NMR_loss_or_UL_gprof_out <- NMR_loss_or_UL_gprof_out[order(NMR_loss_or_UL_gprof_out$p_value), ]

# Manual inspection of enrichment results, rows 7-22 (by p-value rank)
View(NMR_loss_or_UL_gprof_out[7:(7 + 15), c("term_name", "p_value", "intersection_size", "intersection"), ])

# ---------------------------------------------------------------------
# Save all four enrichment result sets together for reuse without
# re-running gprofiler queries
# ---------------------------------------------------------------------
gprofiler_outs <- list(all_loss_gprof_out = all_loss_gprof_out,
                       all_loss_or_UL_gprof_out = all_loss_or_UL_gprof_out,
                       NMR_loss_genes_gprof_out = NMR_loss_genes_gprof_out,
                       NMR_loss_or_UL_gprof_out = NMR_loss_or_UL_gprof_out)
saveRDS(gprofiler_outs, file = "Figure3_Sperm_Loss_Gprofiler_outs.rds")

# Reload (starting point for downstream candidate gene analysis below)
gprofiler_outs <- readRDS("Figure3_Sperm_Loss_Gprofiler_outs.rds")
NMR_loss_genes_gprof_out <- gprofiler_outs$NMR_loss_genes_gprof_out
NMR_loss_or_UL_gprof_out <- gprofiler_outs$NMR_loss_or_UL_gprof_out

# ---------------------------------------------------------------------
# Full losses: pull out genes belonging to two specific, tightly
# sperm-related GO terms from the strict NMR-loss enrichment result
# ---------------------------------------------------------------------
# Full losses
full_loss <- NMR_loss_genes_gprof_out[NMR_loss_genes_gprof_out$term_name %in% c("sperm flagellum", "motile cilium"), ]
full_loss_genes <- unique(unlist(strsplit(full_loss$intersection, ",")))

# Broader candidate set: genes belonging to any of several sperm/
# fertilization-related GO terms, from the relaxed (loss-or-UL)
# enrichment result
full_and_UL_set <- unique(unlist(strsplit(NMR_loss_or_UL_gprof_out[NMR_loss_or_UL_gprof_out$term_name %in% c("sperm flagellum", "motile cilium", "male gamete generation",
                                                                                                             "sperm midpiece", "sexual reproduction", "spermatogenesis",
                                                                                                             "fusion of sperm to egg plasma membrane involved in single fertilization"), "intersection"], ",")))

# ---------------------------------------------------------------------
# Load NMR gene annotation quality metrics/symbols table, to cross-
# reference the candidate lost sperm genes against annotation
# confidence scores and hierarchical gene naming
# ---------------------------------------------------------------------
symbols_score <- read.table("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure2/Figure2DE_chr15_genes/mHetGla1.1.pri_gene_symbols_metrics.tsv", header = TRUE,
                            as.is = TRUE, sep = "\t")

# Quick lookup subset for Adam-family genes
symbols_score_adam <- symbols_score[grep("Adam", symbols_score$unique_gene), ]

# Restrict to genes with a defined hierarchical gene symbol, sorted by
# annotation score (descending = most confident annotations first)
symbols_score_notNA <- symbols_score[!is.na(symbols_score$hierarchical_gene), ]
symbols_score_notNA <- symbols_score_notNA[order(symbols_score_notNA$score, decreasing = TRUE), ]
# symbols_score_notNA <- symbols_score_notNA[!duplicated(symbols_score_notNA$unique_gene),]  # (dedup by unique_gene considered but not applied)

# Drop genes whose "hierarchical_gene" symbol is still an Ensembl
# placeholder ID (ENSHGLG...) rather than a real gene symbol
drop_embl <- symbols_score_notNA[!grepl("ENSHGLG", symbols_score_notNA$hierarchical_gene), ]

# Restrict to genes overlapping the "full loss or UL" sperm-related
# candidate gene set identified above (case-insensitive match)
drop_embl_full_loss <- drop_embl[toupper(drop_embl$hierarchical_gene) %in% toupper(full_and_UL_set), ]

# ---------------------------------------------------------------------
# Flag genes that show testis-specific RNA-seq (StringTie) evidence, as
# these are candidates to exclude from the "gene truly lost" set (since
# testis expression evidence would argue against complete loss)
# ---------------------------------------------------------------------
# Drop GENE in there is testis expression
drop_genes <- c()
for (i in unique(drop_embl_full_loss$hierarchical_gene)) {
  gene <- drop_embl_full_loss[drop_embl_full_loss$hierarchical_gene == i, ]
  print(i)
  if (any(grepl("stringtie_testis", gene$alias))) {
    drop_genes <- c(drop_genes, i)
  }
  # drop_embl_notestis1 <- drop_embl_full_loss[-(grep("stringtie_testis",drop_embl_full_loss$alias)),]
}
keep_genes <- unique(drop_embl_full_loss$hierarchical_gene)[(!unique(drop_embl_full_loss$hierarchical_gene) %in% drop_genes)]

# ---------------------------------------------------------------------
# For the retained candidate genes, identify ones with multiple gene-
# model "copies" in the annotation (possible segmental duplication
# candidates worth manual review)
# ---------------------------------------------------------------------
embl_copies <- c()
for (i in keep_genes) {
  eachcopy <- drop_embl[(toupper(i) == toupper(drop_embl$hierarchical_gene)), ]
  if (nrow(eachcopy) > 1) {
    embl_copies <- rbind(embl_copies, eachcopy)
  }
}

# candidate
dups <- unique(embl_copies$hierarchical_gene)
write.csv(unique(embl_copies$hierarchical_gene), file = "Sperm_SegDup_candidates.csv",
          quote = F, row.names = F, col.names = F)

# All good - though some should be dropped for ISOseq
write.csv(keep_genes, file = "Sperm_SegDup_candidates_full.csv",
          quote = F, row.names = F, col.names = F)
