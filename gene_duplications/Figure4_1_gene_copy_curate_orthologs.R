setwd("~/Desktop/NMR_multi_assembly/v090226/annotations/")
library(pbapply)

# ===========================================================================
# Build three progressively stricter "one representative transcript per
# gene" tables for each assembly (pri, alt, BatSui):
#   Version 1: score > 0, prefer multi-exonic transcript if multiple pass
#   Version 2 ("intact"): score > 16 AND is_complete == "True", same
#              multi-exonic tie-breaking
#   Version 3 ("intact.stringtie"): Version 2 genes further restricted to
#              those with RNA-seq (StringTie) transcriptional support
#
# The same three-step procedure is repeated for pri, alt, and BatSui
# assemblies below, followed by protein FASTA curation for OrthoFinder
# and an assembly-quality (N50) comparison plot.
# ===========================================================================

setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/")
# ---------------------------------------------------------------------
# Version 1: all genes, best transcript (mHetGla pri)
# ---------------------------------------------------------------------
tab <- read.table("mHetGla1.1.pri_gene_symbols_metrics.tsv", header = TRUE, as.is = TRUE, sep = "\t")

# For genes without a resolved hierarchical (symbol-based) gene name,
# fall back to using the raw mikado transcript ID as the gene identifier
tab$hierarchical_gene[is.na(tab$hierarchical_gene)] <- tab$mikado_id[is.na(tab$hierarchical_gene)]

# Version 1: keep multi-exonic copies. and score > 0 genes
tab_fix_copy <- c()
uni <- unique(tab$hierarchical_gene)

### Version 2: intact genes

# For each gene symbol, keep transcripts with score > 0. If exactly one
# passes, return it. If multiple pass, restrict to multi-exonic
# transcripts (as a simple heuristic to avoid partial/fragmentary
# single-exon gene models) and return those if any remain.
single_copy_table <- function(i) {
  # print(i)
  tab1 <- tab[tab$hierarchical_gene == i & tab$score > 0, ]
  if (nrow(tab1) == 1) {
    return(tab1)
  } else {
    tab1 <- tab1[tab1$highest_cds_exon_number > 1, ]
    if (nrow(tab1) > 0) {
      return(tab1)
    }
  }
}
single_copy_genes <- as.data.frame(do.call("rbind", pblapply(uni, single_copy_table)))
single_copy_genes <- single_copy_genes[!is.na(single_copy_genes$tid), ]
write.table(single_copy_genes, file = "mHetGla.pri.version1.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# Reload from disk
single_copy_genes <- read.table("mHetGla.pri.version1.tsv", header = T, as.is = T, sep = "\t")

# ---------------------------------------------------------------------
# Version 2: score > 16 and complete gene. if multiple copies, keep multi-exonic
# ---------------------------------------------------------------------
tab_fix_copy <- c()
uni <- unique(tab$hierarchical_gene)

### Version 2: intact genes
single_copy_table <- function(i) {
  # print(i)
  tab1 <- tab[tab$hierarchical_gene == i & tab$is_complete == "True" & tab$score > 16, ]
  if (nrow(tab1) == 1) {
    return(tab1)
  } else {
    tab1 <- tab1[tab1$highest_cds_exon_number > 1, ]
    if (nrow(tab1) > 0) {
      return(tab1)
    }
  }
}
single_copy_genes <- as.data.frame(do.call("rbind", pblapply(uni, single_copy_table)))
write.table(single_copy_genes, file = "mHetGla.pri.intact.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# ---------------------------------------------------------------------
# Version 3: with transcription (pri)
# Restrict Version 2 genes to those with at least one StringTie-derived
# transcript among their gene models, as evidence of active transcription
# ---------------------------------------------------------------------
mHetGla1.pri <- read.table("mHetGla1.pri.metrics.tsv", header = TRUE, sep = "\t",
                           as.is = TRUE)
genes1 <- unique(single_copy_genes$mikado_id)

stringtie_fun <- function(i) {
  string_trans <- mHetGla1.pri[mHetGla1.pri$parent == i, ]
  return("stringtie" %in% string_trans$original_source)
}
genes_keep <- unlist(pblapply(genes1, stringtie_fun))
head(mHetGla1.pri)

genes_version3 <- genes1[genes_keep == TRUE]
single_copy_genes_version3 <- single_copy_genes[single_copy_genes$mikado_id %in% genes_version3, ]
head(genes1)
write.table(single_copy_genes_version3, file = "mHetGla.pri.intact.stringtie.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

##
##### mHetGla.Alt
##
# ---------------------------------------------------------------------
# Repeat the exact same 3-version procedure for the alt haplotype
# assembly
# ---------------------------------------------------------------------


# Version 1: all genes, best transcript
tab <- read.table("mHetGla.alt_gene_symbols_metrics.tsv", header = TRUE, as.is = TRUE, sep = "\t")
tab$hierarchical_gene[is.na(tab$hierarchical_gene)] <- tab$mikado_id[is.na(tab$hierarchical_gene)]

# Version 1: keep multi-exonic copies. and score > 0 genes
tab_fix_copy <- c()
uni <- unique(tab$hierarchical_gene)

### Version 2: intact genes
single_copy_table <- function(i) {
  # altnt(i)
  tab1 <- tab[tab$hierarchical_gene == i & tab$score > 0, ]
  if (nrow(tab1) == 1) {
    return(tab1)
  } else {
    tab1 <- tab1[tab1$highest_cds_exon_number > 1, ]
    if (nrow(tab1) > 0) {
      return(tab1)
    }
  }
}
single_copy_genes <- as.data.frame(do.call("rbind", pblapply(uni, single_copy_table)))
write.table(single_copy_genes, file = "mHetGla.alt.version1.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# Version 2: score > 16 and complete gene. if multiple copies, keep multi-exonic
tab_fix_copy <- c()
uni <- unique(tab$hierarchical_gene)

### Version 2: intact genes
single_copy_table <- function(i) {
  # altnt(i)
  tab1 <- tab[tab$hierarchical_gene == i & tab$is_complete == "True" & tab$score > 16, ]
  if (nrow(tab1) == 1) {
    return(tab1)
  } else {
    tab1 <- tab1[tab1$highest_cds_exon_number > 1, ]
    if (nrow(tab1) > 0) {
      return(tab1)
    }
  }
}
single_copy_genes <- as.data.frame(do.call("rbind", pblapply(uni, single_copy_table)))
write.table(single_copy_genes, file = "mHetGla.alt.intact.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

single_copy_genes <- read.table("mHetGla.alt.intact.tsv", header = T,
                                as.is = T, sep = "\t")

# ---------------------------------------------------------------------
# Version 3: with transcription (alt)
# ---------------------------------------------------------------------
mHetGla1.alt <- read.table("mHetGla1.alt.metrics.tsv", header = TRUE, sep = "\t",
                           as.is = TRUE)
genes1 <- unique(single_copy_genes$mikado_id)

stringtie_fun <- function(i) {
  string_trans <- mHetGla1.alt[mHetGla1.alt$parent == i, ]
  return("stringtie" %in% string_trans$original_source)
}
genes_keep <- unlist(pblapply(genes1, stringtie_fun))
head(mHetGla1.alt)

genes_version3 <- genes1[genes_keep == TRUE]
single_copy_genes_version3 <- single_copy_genes[single_copy_genes$mikado_id %in% genes_version3, ]
head(genes1)
write.table(single_copy_genes_version3, file = "mHetGla.alt.intact.stringtie.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

##
######
# ---------------------------------------------------------------------
# Repeat Versions 1-2 (no Version 3/stringtie step here) for the BatSui
# assembly
# ---------------------------------------------------------------------
library(pbapply)

tab <- read.table("mBatSui1.1.pri/BatSui_gene_symbols_metrics.tsv", header = TRUE, as.is = TRUE, sep = "\t")
tab$hierarchical_gene[is.na(tab$hierarchical_gene)] <- tab$mikado_id[is.na(tab$hierarchical_gene)]

# Version 1: score > 0. If multiple copies, keep multi-exonic
# Version 2: score > 16 and complete gene. if multiple copies, keep multi-exonic
tab_fix_copy <- c()
uni <- unique(tab$hierarchical_gene)

### Version 2: intact genes
single_copy_table <- function(i) {
  # print(i)
  tab1 <- tab[tab$hierarchical_gene == i & tab$score > 0, ]
  if (nrow(tab1) == 1) {
    return(tab1)
  } else {
    tab1 <- tab1[tab1$highest_cds_exon_number > 1, ]
    if (nrow(tab1) > 0) {
      return(tab1)
    }
  }
}
single_copy_genes <- as.data.frame(do.call("rbind", pblapply(uni, single_copy_table)))
write.table(single_copy_genes, file = "BatSui.version1.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# Version 2: score > 16 and complete gene. if multiple copies, keep multi-exonic
tab_fix_copy <- c()
uni <- unique(tab$hierarchical_gene)

### Version 2: intact genes
single_copy_table <- function(i) {
  # print(i)
  tab1 <- tab[tab$hierarchical_gene == i & tab$is_complete == "True" & tab$score > 16, ]
  if (nrow(tab1) == 1) {
    return(tab1)
  } else {
    tab1 <- tab1[tab1$highest_cds_exon_number > 1, ]
    if (nrow(tab1) > 0) {
      return(tab1)
    }
  }
}
single_copy_genes <- as.data.frame(do.call("rbind", pblapply(uni, single_copy_table)))
write.table(single_copy_genes, file = "BatSui.intact.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

##
##### Part 2, curating the "table" and faa files for orthofinder
##
library(Biostrings)
setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure4/Figure4A_Duplications/gene_copy_orthofinder/referece_proteins/")

##
##### Mouse C57
##
# ---------------------------------------------------------------------
# For the C57 mouse reference proteome: pick the longest protein
# isoform per gene symbol (using AA sequence width as a proxy for
# CDS/protein length), restrict to genes where that longest transcript
# is actually present in the protein FASTA, and write out a cleaned
# one-protein-per-gene FASTA plus a reference lookup table
# ---------------------------------------------------------------------
C57_T2T_protein <- readAAStringSet("C57_T2T.nostop.protein.faa")
c57_widths <- C57_T2T_protein@ranges@width
names(c57_widths) <- C57_T2T_protein@ranges@NAMES

c57_names <- read.table("C57_T2T.table.txt", header = T, as.is = T, sep = "\t")
uni <- unique(c57_names$geneSymbol)

# For each gene symbol, find the transcript ID with the longest protein
longest_c57_transcript <- function(i) {
  transcripts <- c57_names$transcriptID[c57_names$geneSymbol == i]
  y <- names(which.max(c57_widths[transcripts]))
  return(c(i, y))
}
c57_with_longest <- do.call("rbind", pblapply(uni, longest_c57_transcript))
c57.df <- as.data.frame(c57_with_longest, stringsAsFactors = FALSE)
colnames(c57.df) <- c("geneSymbol", "transcriptID")
rownames(c57.df) <- c57.df$transcriptID

# Restrict to genes whose chosen longest transcript is actually present
# in the protein FASTA, clean up "rna-" ID prefixes, and rename
# sequences by gene symbol before writing the final protein FASTA
inter <- intersect(rownames(c57.df), names(C57_T2T_protein))
c57.inter <- c57.df[inter, ]
c57.inter$transcriptID <- gsub("rna-", "", c57.inter$transcriptID)
C57_T2T_protein.inter <- C57_T2T_protein[inter]
names(C57_T2T_protein.inter) <- gsub("rna-", "", names(C57_T2T_protein.inter))

write.table(c57.inter, file = "Mmus_reference.txt",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

names(C57_T2T_protein.inter) <- c57.inter$geneSymbol
Biostrings::writeXStringSet(C57_T2T_protein.inter,
                            file = "Mmus_proteins.faa")

##
##### Hsap T2T
##
# ---------------------------------------------------------------------
# Identical longest-isoform-per-gene procedure, applied to the human
# T2T-CHM13v2 reference proteome
# ---------------------------------------------------------------------
T2T_CHM14v2_protein <- readAAStringSet("T2T_CHM14v2.nostop.protein.faa")
c57_widths <- T2T_CHM14v2_protein@ranges@width
names(c57_widths) <- T2T_CHM14v2_protein@ranges@NAMES

c57_names <- read.table("T2T_CHM14v2.table.txt", header = T, as.is = T, sep = "\t")
uni <- unique(c57_names$geneSymbol)

longest_c57_transcript <- function(i) {
  transcripts <- c57_names$transcriptID[c57_names$geneSymbol == i]
  y <- names(which.max(c57_widths[transcripts]))
  return(c(i, y))
}
c57_with_longest <- do.call("rbind", pblapply(uni, longest_c57_transcript))
c57.df <- as.data.frame(c57_with_longest, stringsAsFactors = FALSE)
colnames(c57.df) <- c("geneSymbol", "transcriptID")
rownames(c57.df) <- c57.df$transcriptID

inter <- intersect(rownames(c57.df), names(T2T_CHM14v2_protein))
c57.inter <- c57.df[inter, ]
c57.inter$transcriptID <- gsub("rna-", "", c57.inter$transcriptID)
T2T_CHM14v2_protein.inter <- T2T_CHM14v2_protein[inter]
names(T2T_CHM14v2_protein.inter) <- gsub("rna-", "", names(T2T_CHM14v2_protein.inter))

write.table(c57.inter, file = "Hsap_reference.txt",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

names(T2T_CHM14v2_protein.inter) <- c57.inter$geneSymbol
Biostrings::writeXStringSet(T2T_CHM14v2_protein.inter,
                            file = "Hsap_proteins.faa")

##
##### NMR
##
library(Biostrings)

# Version 1 - still get rid of single-exon copies, just don't have the cutoff in version 2



for (i in c("pri","alt")) {
  # ---------------------------------------------------------------------
  # Version 1: intersect the cleaned protein FASTA with the Version 1
  # gene table, rename sequences first by mikado ID (saved as one
  # FASTA), then by unique gene symbol (saved as a second FASTA)
  # ---------------------------------------------------------------------
  mHetGla <- readAAStringSet(paste0("mHetGla1.1.", i, ".cleaned.protein.faa"))
  symbols <- read.table(paste0("mHetGla.", i, ".version1.tsv"),
                        header = TRUE, as.is = TRUE, sep = "\t")
  symbols <- symbols[!is.na(symbols$tid), ]
  rownames(symbols) <- symbols$tid
  inter_version1 <- intersect(names(mHetGla), rownames(symbols))
  mHetGla.version1 <- mHetGla[inter_version1]
  symbols.version1 <- symbols[inter_version1, ]
  
  names(mHetGla.version1) <- symbols.version1$mikado_id
  Biostrings::writeXStringSet(mHetGla.version1,
                              file = paste0("mHetGla.", i, ".mikadoID.version1.faa"))
  
  names(mHetGla.version1) <- symbols.version1$unique_gene
  Biostrings::writeXStringSet(mHetGla.version1,
                              file = paste0("mHetGla.", i, ".version1.faa"))
  
  # ---------------------------------------------------------------------
  # Version 2: same procedure using the "intact" (score>16, complete)
  # gene table
  # ---------------------------------------------------------------------
  symbols.version2 <- read.table(paste0("mHetGla.", i, ".intact.tsv"),
                                 header = TRUE, as.is = TRUE, sep = "\t")
  symbols.version2 <- symbols.version2[!is.na(symbols.version2$tid), ]
  rownames(symbols.version2) <- symbols.version2$tid
  inter_version2 <- intersect(names(mHetGla), rownames(symbols.version2))
  mHetGla.version2 <- mHetGla[inter_version2]
  symbols.version2.inter <- symbols.version2[inter_version2, ]
  
  names(mHetGla.version2) <- symbols.version2.inter$mikado_id
  Biostrings::writeXStringSet(mHetGla.version2,
                              file = paste0("mHetGla.", i, ".mikadoID.version2.faa"))
  
  names(mHetGla.version2) <- symbols.version2.inter$unique_gene
  Biostrings::writeXStringSet(mHetGla.version2,
                              file = paste0("mHetGla.", i, ".version2.faa"))
  
  # ---------------------------------------------------------------------
  # Version 3: same procedure using the "intact.stringtie"
  # (transcription-supported) gene table
  # ---------------------------------------------------------------------
  symbols.version3 <- read.table(paste0("mHetGla.", i, ".intact.stringtie.tsv"),
                                 header = TRUE, as.is = TRUE, sep = "\t")
  symbols.version3 <- symbols.version3[!is.na(symbols.version3$tid), ]
  rownames(symbols.version3) <- symbols.version3$tid
  inter_version3 <- intersect(names(mHetGla), rownames(symbols.version3))
  mHetGla.version3 <- mHetGla[inter_version3]
  symbols.version3.inter <- symbols.version3[inter_version3, ]
  
  names(mHetGla.version3) <- symbols.version3.inter$mikado_id
  Biostrings::writeXStringSet(mHetGla.version3,
                              file = paste0("mHetGla.", i, ".mikadoID.version3.faa"))
  
  names(mHetGla.version3) <- symbols.version3.inter$unique_gene
  Biostrings::writeXStringSet(mHetGla.version3,
                              file = paste0("mHetGla.", i, ".version3.faa"))
}


##
##### CDMR
##
### Version 1 - still get rid of copies, just don't have the cutoff in version 2
# ---------------------------------------------------------------------
# Same Version 1/2 protein FASTA curation, applied to the BatSui
# assembly. Only two versions here (no stringtie-supported Version 3),
# consistent with only Versions 1-2 having been built for BatSui above.
# ---------------------------------------------------------------------
BatSui <- readAAStringSet("BatSuiPri.protein.faa")
symbols <- read.table("BatSui.version1.tsv",
                      header = TRUE, as.is = TRUE, sep = "\t")
rownames(symbols) <- symbols$tid
inter_version1 <- intersect(names(BatSui), rownames(symbols))
# rownames(symbols.version2) <- BatSui[inter_version1]
symbols.version1 <- symbols[inter_version1, ]
BatSui.version1 <- BatSui[inter_version1]
names(BatSui.version1) <- symbols.version1$unique_gene
Biostrings::writeXStringSet(BatSui.version1,
                            file = "BatSui.version1.faa")

# Version 2
symbols.version2 <- read.table("BatSui.intact.tsv",
                               header = TRUE, as.is = TRUE, sep = "\t")
rownames(symbols.version2) <- symbols.version2$tid
inter_version2 <- intersect(names(BatSui), rownames(symbols.version2))
BatSui.version2 <- BatSui[inter_version2]
symbols.version2.inter <- symbols.version2[inter_version2, ]
names(BatSui.version2) <- symbols.version2.inter$unique_gene
Biostrings::writeXStringSet(BatSui.version2,
                            file = "BatSui.version2.faa")
