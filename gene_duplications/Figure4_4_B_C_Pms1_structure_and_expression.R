library(rtracklayer)

# ===========================================================================
#
# This is simplified version of the Figure4_3 script with the finalized
# post-processed objects
#
# Figure 4 (Pms1 gene family + DDR gene duplications) support script:
#   1. Build a color-coded BED track of NMR gene models (protein-coding,
#      pseudogene, ncRNA) for genome browser display (repeat of the
#      earlier "gene color BED" script seen previously in this session)
#   2. Extract Pms1 gene copies and their protein sequences
#   3. Build a StringTie-transcript-based tissue-expression presence/
#      absence heatmap specifically for Pms1 gene copies
#   4. Cross-check a curated list of gene copy-number calls against the
#      annotation metrics table
#   5. Build rRNA BED tracks for BatSui and NMR
#   6. RNA-seq QC analysis (RIN/DV200) across species, with ANOVA/Tukey
#      post-hoc tests
#   7. Generalize the Pms1-style tissue-expression heatmap approach to a
#      broader set of DNA-damage-repair (DDR) gene duplications
#
# ===========================================================================

# ---------------------------------------------------------------------
# Load per-transcript annotation metrics and gene-level GFF features,
# and build the color-coded gene BED track (identical logic to the
# "mHetGlaPri_gene_color.bed" script from earlier in this session)
# ---------------------------------------------------------------------
metrics <- read.table("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/mHetGla1.pri_annotation_metrics_with_pseudogenes_v120426.tsv.gz",
                      header = TRUE, as.is = TRUE, sep = "\t")
gff <- readGFF("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/mHetGla1.pri.MT.v020326.gff.gz")
gff <- as.data.frame(gff[gff$type %in% c("gene", "ncRNA_gene"), ])

metrics_gene <- metrics[!duplicated(metrics$parent), ]
metrics_gene <- metrics_gene[, c("Name", "parent", "gene_type", "is_complete")]
metrics_bed <- gff[, c("seqid", "start", "end", "strand", "ID")]
merged_bed <- merge(metrics_bed, metrics_gene, by.x = "ID", by.y = "parent")
merged_bed$score <- 0
merged_bed$length <- merged_bed$end - merged_bed$start

# Make RGB: green = complete protein-coding gene, orange = pseudogene/
# incomplete, lavender = ncRNA gene
rgb <- rep("0,138,0", nrow(merged_bed))
rgb[merged_bed$gene_type == "pseudogene" | merged_bed$is_complete == "False"] <- "211,84,0"
rgb[merged_bed$gene_type == "ncRNA_gene"] <- "216,191,216"
merged_bed$rgb <- rgb

mHetGla1Pri_color_bed <- merged_bed[, c("seqid", "start", "end", "Name", "score", "strand", "start", "end", "rgb")]
mHetGla1Pri_color_bed2 <- mHetGla1Pri_color_bed
mHetGla1Pri_color_bed2$Name[mHetGla1Pri_color_bed2$rgb != "0,138,0"] <- "."

write.table(mHetGla1Pri_color_bed, file = "gene_copy_analysis/mHetGlaPri_gene_color.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Clean gene names (strip mikado prefix, shorten Ensembl NMR IDs), write
# a version keeping "-copy" suffixes, then a version collapsing them
mHetGla1Pri_color_bed2$Name <- gsub("mikado.", "", mHetGla1Pri_color_bed2$Name)
mHetGla1Pri_color_bed2$Name <- gsub("ENSHGLG000000", "Ens",
                                    mHetGla1Pri_color_bed2$Name)
write.table(mHetGla1Pri_color_bed2, file = "gene_copy_analysis/mHetGlaPri_PcGName_withcopy.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

mHetGla1Pri_color_bed2$Name <- unlist(lapply(strsplit(mHetGla1Pri_color_bed2$Name, "-copy"), function(x) return(x[1])))
write.table(mHetGla1Pri_color_bed2, file = "gene_copy_analysis/mHetGlaPri_PcGName.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


##
##### Same bed file but for BatSui
##
# ---------------------------------------------------------------------

alias <- metrics$alias
alias <- alias[grep("stringtie", alias)]
alias <- gsub("stringtie_", "", alias)
before <- sub("\\..*", "", alias)

metrics <- read.table("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/mBatSui1.1.pri/BatSui_gene_symbols_metrics.tsv",
                      header = TRUE, as.is = TRUE, sep = "\t")
alias <- metrics$alias
alias <- alias[grep("stringtie", alias)]
alias <- gsub("stringtie_", "", alias)

## Figure 4C
## Pms1 copy locations

# ---------------------------------------------------------------------
# Extract Pms1 protein copies from the "version2" NMR pri protein set
# ---------------------------------------------------------------------
library(Biostrings)
v2 <- readAAStringSet("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure4/Figure4A_Duplications/gene_copy_orthofinder/referece_proteins/mHetGla.pri.version2.faa")
writeXStringSet(v2[grep("Pms1", names(v2))], file = "mHetGla1.1.pri.version2.pms1.faa")


# Figure 4C gene expression
# ---------------------------------------------------------------------
# Load a table mapping Pms1 gene copies to their StringTie transcript
# group IDs, then load the full StringTie transcript GFF and match
# transcripts (via their Parent gene ID) to these Pms1 copies
# ---------------------------------------------------------------------
pms1ID <- read.table("../Figure4BC_Pms1/pms1_gene_group_id.txt", header = TRUE, as.is = TRUE)
library(rtracklayer)
string <- readGFF("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/mHetGla1.pri.stringtie.s.gff.gz")
library(pbapply)
stringDF <- as.data.frame(string)

# The "Parent" column from a GFF read via rtracklayer is typically a
# CharacterList (each row can have 0 or more parents) - this coerces it
# to a flat character vector, using "" for rows with no parent, so it
# can be used for simple string matching below
result <- pbvapply(stringDF$Parent, function(z) if (length(z)) z else "", character(1))
stringDF$ParentChar <- result

# Restrict to StringTie records whose parent gene ID matches one of the
# known Pms1 copy IDs, then further restrict to transcript-level rows
# and extract tissue of origin from the transcript ID (text before the
# first ".")
stringDFpms1 <- stringDF[stringDF$ParentChar %in% pms1ID$StringTieID, ]
stringDFpms1.transcript <- stringDFpms1[stringDFpms1$type == "transcript", ]
stringDFpms1.transcript$tissue <- sub("\\..*", "", stringDFpms1.transcript$ID)

# ---------------------------------------------------------------------
# Build a Pms1-copy x tissue presence/absence matrix (1 = at least one
# transcript detected in that tissue for that copy), then plot as a
# black/white heatmap
# ---------------------------------------------------------------------
library(pheatmap)
colors <- colorRampPalette(c("white", "black"))
tables <- table(stringDFpms1.transcript$ParentChar, stringDFpms1.transcript$tissue)


tables[tables > 0] <- 1

# Manually relabel rows (Pms1 copy StringTie group IDs) with more
# readable copy names
rownames(tables) <- c("copy8",
                      "copy5.6",
                      "copy1.7",
                      "copy2.3",
                      "copy0",
                      "copy4")


p1 <- pheatmap(tables, color = colors(2), cellwidth = 10, cellheight = 10,
               cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45, legend = FALSE)

pdf(file = "Figure4C_Pms1_expression.pdf", width = 8, height = 11)
print(p1)
dev.off()

pdf(file = "Figure4C_Pms1_expression.pdf", width = 8, height = 11)
print(p1)
dev.off()

ggsave("Pms1_Isoforms.png", p1)

pullString <- pblapply(string$Parent, function(x) return(x %in% pms1ID$StringTieID))
SP <- lengths(string$Parent)
parent1 <- 1:length(SP)

stringPms <- string$Parent[pullString, ]

## Figure 4D recheck copies
#
# ---------------------------------------------------------------------
# Cross-check a curated gene copy-number table ("mHetGla_copies.csv")
# against the annotation metrics, by grep-matching gene symbols and
# counting how many "complete, protein-coding" gene records match each
# curated gene name - useful as a sanity check that manually curated
# copy counts agree with what's in the annotation
# ---------------------------------------------------------------------
copies <- read.csv("gene_copy_analysis/mHetGla_copies.csv", header = T, as.is = T)
copies <- copies[complete.cases(copies), ]
copies$Gene <- tolower(copies$Gene)
library(stringr)
copies$Gene <- str_to_title(copies$Gene)

rows <- c()
for (x in 1:length(copies$Gene)) {
  i = copies$Gene[x]
  print(i)
  ge <- metrics_gene[grep(i, metrics_gene$Name), ]
  if (nrow(ge) == 0) next
  ge$fin <- FALSE
  ge$fin[ge$is_complete == "True" & ge$gene_type == "gene"] <- TRUE
  nrow(ge[ge$fin, ])  # (autoprint - not stored)
  row <- c(i, copies$HetGla.Copies[x], nrow(ge[ge$fin, ]))
  rows <- rbind(rows, row)
}
rows <- as.data.frame(rows)
colnames(rows) <- c("gene", "copies", "grepped")

# ---------------------------------------------------------------------
# Build rRNA BED tracks for BatSui and NMR, from separate rRNA
# annotation GFFs
# ---------------------------------------------------------------------
library(rtracklayer)
rRNA <- readGFF("Figure1E_annotations/mBatSui1.1.pri/mBatSui.1.1.primary.rrna.s.gff")
table(rRNA$seqid, rRNA$Name)
rRNA_bed <- rRNA[, c("seqid", "start", "end", "Name", "score", "strand")]
write.table(rRNA_bed, file = "Figure1E_annotations/mBatSui1.1.pri/mBatSui.1.1.primary.rrna.s.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


mHetGla <- readGFF("Figure1E_annotations/mHetGla1.1.pri.rrna.gff")
rRNA_bed <- mHetGla[, c("seqid", "start", "end", "Name", "score", "strand")]
write.table(rRNA_bed, file = "Figure1E_annotations/mHetGla.1.1.pri.rrna.s.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

##
##### Figure4C associated - Other DDR genes
##
# ---------------------------------------------------------------------
# Generalize the Pms1 tissue-expression heatmap approach to a broader
# panel of DNA-damage-repair (DDR) gene duplications: for each DDR gene
# family, build a copy x tissue presence/absence matrix from StringTie
# transcript evidence, plot as a heatmap, and export a combined long-
# format expression matrix across all gene families
# ---------------------------------------------------------------------
library(rtracklayer)
library(pbapply)
setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure4/Figure4BC_Pms1/")

DDR_genes <- read.table("DDR_specific_duplications_gene_and_Stringtie.txt", header = TRUE, as.is = TRUE, sep = "\t")
string <- readGFF("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/mHetGla1.pri.stringtie.s.gff.gz")

## Lets start by getting a list of all possible tissues
stringDF <- as.data.frame(string)
result <- pbvapply(stringDF$Parent, function(z) if (length(z)) z else "", character(1))
result  # (autoprint)
stringDF$ParentChar <- result
stringDF$tissue <- sub("\\..*", "", stringDF$ID)
transcriptDF <- stringDF[stringDF$type == "transcript", ]
tissues <- unique(stringDF$tissue)[!is.na(unique(stringDF$tissue))]
head(DDR_genes)

# ---------------------------------------------------------------------
# For each DDR gene family (grouped by "Main.Symbol"), build a gene-
# copy x tissue binary expression matrix based on whether any
# transcript for that copy's StringTie ID was detected in each tissue
# ---------------------------------------------------------------------
# now we populate all the stringtie tissues for each gene
tissue_mats <- list()
for (i in unique(DDR_genes$Main.Symbol)) {
  gene <- DDR_genes[DDR_genes$Main.Symbol == i, ]
  rownames(gene) <- gene$Gene.name
  mat <- matrix(0, nrow = length(gene$Gene.name), ncol = length(tissues))
  rownames(mat) <- gene$Gene.name
  colnames(mat) <- tissues
  for (j in rownames(gene)) {
    stringID <- gene[j, "StringTie.ID"]
    gene_tissue <- transcriptDF$tissue[transcriptDF$ParentChar %in% stringID]
    mat[j, gene_tissue] <- 1
  }
  tissue_mats[[i]] <- mat
}

# ---------------------------------------------------------------------
# Plot a heatmap per DDR gene family; special-case Pms1 to enforce a
# specific, meaningful row order (matching copy naming used earlier in
# this script) rather than the default alphabetical order
# ---------------------------------------------------------------------
library(pheatmap)
library(ggplot2)
colors <- colorRampPalette(c("white", "black"))
dir.create("DDR_segdups")

for (i in names(tissue_mats)) {
  tables <- tissue_mats[[i]]
  if (i == "Pms1") {
    tables <- tables[c("Pms1-c5", "Pms1-c1",
                       "Pms1-c2", "Pms1-c0",
                       "Pms1-c4", "Pms1-c8"), ]
  }
  p1 <- pheatmap(tables, color = colors(2), cellwidth = 10, cellheight = 10,
                 cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 45, legend = FALSE)
  
  ggsave(paste0("DDR_segdups/", i, "_gene_models.png"), p1, width = 8.5, height = 11.5)
}

# ---------------------------------------------------------------------
# Combine all per-gene-family matrices into one long-format table
# (each row = one gene copy in one gene family, columns = tissues),
# tag each row with its gene family and a cleaned individual gene name,
# and export as a single combined CSV
# ---------------------------------------------------------------------
DDR_expression <- do.call(rbind, Map(function(x, nm) {
  cbind(as.data.frame(x), tissue = nm)
}, tissue_mats, names(tissue_mats)))
colnames(DDR_expression)[ncol(DDR_expression)] <- "gene"
DDR_expression$gene_name <- unlist(lapply(strsplit(rownames(DDR_expression), split = "\\."), function(x) return(x[2])))

write.csv(DDR_expression, file = "DDR_segdups/DDR_SegDup_Expression_Binary_matrix.csv",
          quote = FALSE, row.names = TRUE, col.names = TRUE)
