library(rtracklayer)
setwd("~/Desktop/NMR_multi_assembly/v090226/annotations/")

# This script annotates and analyzes PGRE pseudogenes against the rest of the annotations.
# This script was repeated without chnanges for mHetGla1.1.alt and mBatsui1.pri. Just swap out the gff and metrics.tsv files in 
# Figure1E_annotations for Figure1E_annotations/mHetGla1.1.alt and Figure1E_annotations/mBatSui1.1.pri

# ---------------------------------------------------------------------
# Load the primary mitochondrial (MT) annotation GFF and extract
# transcript-level features (mRNA and ncRNA) as a data frame
# ---------------------------------------------------------------------
annotation_gff <- readGFF("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/mHetGla1.pri.MT.v020326.gff.gz")
annotation_transcripts <- as.data.frame(annotation_gff[annotation_gff$type %in% c("ncRNA", "mRNA"), ])
rownames(annotation_transcripts) <- annotation_transcripts$ID
annotation_transcripts <- annotation_transcripts[, c("ID", "type")]

# Load the processed/unprocessed pseudogene GFF (from a separate pipeline output)
pseudo_gff <- readGFF("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/mHetGla1.1.pri.PGRe.gff")

# ---------------------------------------------------------------------
# Load per-transcript annotation quality metrics (e.g., score,
# completeness) and restrict to transcripts present in the MT GFF
# ---------------------------------------------------------------------
annotation_metrics <- read.delim("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/mHetGla1.pri.metrics.tsv", header = TRUE, as.is = TRUE, sep = "\t")
annotation_metrics <- annotation_metrics[annotation_metrics$tid %in% annotation_gff$ID, ]
head(annotation_metrics)

# Reduce pseudogene GFF to one row per parent gene (dedupe) and build
# a lookup table of pseudogene type by parent gene ID
pseudo_gff_simple <- pseudo_gff[!duplicated(pseudo_gff$parent_gene), ]
pseudo_df <- as.data.frame((pseudo_gff_simple[, c("type", "parent_gene")]))
rownames(pseudo_df) <- pseudo_df$parent_gene

# ---------------------------------------------------------------------
# Assign a "type" to each transcript in annotation_metrics: use the
# pseudogene type if the transcript's ID is a known pseudogene,
# otherwise fall back to the mRNA/ncRNA type from the GFF
# ---------------------------------------------------------------------
annotation_metrics$type <- ""
for (i in 1:length(annotation_metrics$tid)) {
  tid <- annotation_metrics$tid[i]
  isPseudo <- as.character(pseudo_df[tid, "type"])
  if (!is.na(isPseudo)) {
    annotation_metrics$type[i] <- isPseudo
  } else {
    annotation_metrics$type[i] <- as.character(annotation_transcripts[tid, "type"])
  }
}

# Flag transcripts passing a minimum annotation score threshold (>16)
annotation_metrics$minScore <- FALSE
annotation_metrics$minScore[annotation_metrics$score > 16] <- TRUE

# ---------------------------------------------------------------------
# Plot: distribution of annotation scores, colored by transcript type
# ---------------------------------------------------------------------
library(ggplot2)
ggplot(annotation_metrics, aes(x = score, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  theme_minimal()

# ---------------------------------------------------------------------
# Summarize, per transcript type, the proportion passing the score
# threshold and the proportion flagged as "complete"
# ---------------------------------------------------------------------
gene_completion <- annotation_metrics[, c("type", "is_complete", "minScore")]
gene_completion$is_complete[gene_completion$is_complete == "True"] <- TRUE
gene_completion$is_complete[gene_completion$is_complete == "False"] <- FALSE

tab <- cbind(table(gene_completion$type, gene_completion$minScore),
             table(gene_completion$type, gene_completion$is_complete))
colnames(tab) <- c("score.lt.16", "score.gt.16", "complete.F", "complete.T")
tab <- as.data.frame(tab)

# Compute proportions: fraction scoring >16 and fraction complete, per type
tab$score.P <- tab$score.gt.16 / (tab$score.gt.16 + tab$score.lt.16)
tab$complete.P <- tab$complete.T / (tab$complete.T + tab$complete.F)
df_plot <- tab

# ---------------------------------------------------------------------
# Reshape the proportions to long format and plot as a grouped bar chart
# (score-pass proportion vs completeness proportion, per gene type)
# ---------------------------------------------------------------------
library(tidyr)
library(dplyr)

df_plot$total <- with(df_plot,
                      score.gt.16 + score.lt.16   # or complete.T + complete.F (same totals)
)
df_plot$type <- rownames(df_plot)

df_long <- df_plot %>%
  select(type, score.P, complete.P) %>%
  pivot_longer(
    cols = c(score.P, complete.P),
    names_to = "metric",
    values_to = "proportion"
  )

df_long <- df_long %>%
  left_join(df_plot[, c("type", "total")], by = "type")

library(ggplot2)
ggplot(df_long, aes(x = type, y = proportion, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # thinner bars
  theme_minimal(base_size = 10) +  # smaller overall base
  labs(x = "Type", y = "Proportion", fill = "") +
  ylim(0, 1) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  # bigger x labels
    axis.text.y = element_text(size = 14),                         # bigger y labels
    axis.title = element_text(size = 16),                          # bigger axis titles
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    plot.margin = margin(5, 5, 5, 5)  # tighter margins
  )

# Quick inspection of the raw GFF object
annotation_gff

# ---------------------------------------------------------------------
# Reclassify incomplete mRNAs (no ORF) as a distinct pseudogene category
# ---------------------------------------------------------------------
annotation_metrics$type[annotation_metrics$type == "mRNA" & annotation_metrics$is_complete == "False"] <- "pseudogene_noORF"

# annotation_metrics$gene_type <- ""

# ---------------------------------------------------------------------
# For each parent gene, decide a gene-level classification based on
# the types of its child transcripts: protein-coding gene, ncRNA gene,
# or pseudogene
# ---------------------------------------------------------------------
parents <- unique(annotation_metrics$parent)
parent_call <- c()
for (g in parents) {
  allType <- annotation_metrics$type[annotation_metrics$parent == g]
  if ("mRNA" %in% allType) {
    parent_call <- c(parent_call, "gene")
    next
  }
  if ("ncRNA" %in% allType) {
    parent_call <- c(parent_call, "ncRNA_gene")
    next
  }
  if (any(grepl("pseudogene", allType))) {
    parent_call <- c(parent_call, "pseudogene")
  }
}

# ---------------------------------------------------------------------
# Barplot: counts of specific ncRNA subtypes (tRNA, snoRNA, etc.) plus
# counts of gene-level classifications (gene/ncRNA_gene/pseudogene),
# restricted to categories with >100 features for readability
# ---------------------------------------------------------------------
each_feature <- c(table(annotation_gff$type)[c("tRNA", "snoRNA", "snRNA", "rRNA", "microRNA", "scaRNA", "sRNA")],
                  table(parent_call))
each_feature <- each_feature[each_feature > 100]
barplot(each_feature, horiz = T, las = 2)

# Sanity check: transcript type counts, and parents of processed pseudogenes
table(annotation_metrics$type)
annotation_metrics$parent[annotation_metrics$type %in% "processed_pseudogene"]

# ---------------------------------------------------------------------
# Build a gene ID -> gene Name lookup table from "gene" features
# ---------------------------------------------------------------------
gene <- annotation_gff[annotation_gff$type == "gene", ]
geneDf <- as.data.frame(gene)
geneDf <- geneDf[c("ID", "Name")]
rownames(geneDf) <- geneDf$ID

# For each transcript, look up its parent gene's Name and gene-level
# classification (gene/ncRNA_gene/pseudogene) determined above
geneName <- c()
geneType <- c()
names(parent_call) <- parents
for (i in annotation_metrics$parent) {
  geneName <- c(geneName, geneDf[i, "Name"])
  geneType <- c(geneType, parent_call[i])
}
annotation_metrics$Name <- geneName
annotation_metrics$gene_type <- geneType

# ---------------------------------------------------------------------
# Subset to a shareable set of columns and write out the final
# per-transcript annotation metrics table (including pseudogenes)
# ---------------------------------------------------------------------
annotation_metrics_share <- annotation_metrics[c("tid", "Name", "parent", "gene_type", "type",
                                                 "alias", "score", "type", "is_complete", "cdna_length", "exon_num", "minScore")]
# Fall back to parent ID as the Name when no gene Name is available
annotation_metrics_share$Name[is.na(annotation_metrics_share$Name)] <- annotation_metrics_share$parent[is.na(annotation_metrics_share$Name)]

write.table(annotation_metrics_share, file = "mHetGla1.pri_annotation_metrics_with_pseudogenes_v120426.tsv", quote = FALSE, row.names = FALSE,
            col.names = TRUE, sep = "\t")

# ---------------------------------------------------------------------
# Processed pseudogenes: strip the "-copyN" suffix to collapse multiple
# copies back to their parent gene name, then tally occurrences
# (i.e., how many processed pseudogene copies exist per source gene)
# ---------------------------------------------------------------------
processed_pseudo <- annotation_metrics$Name[annotation_metrics$type == "processed_pseudogene"]
processed_pseudo <- unlist(lapply(strsplit(processed_pseudo, "-copy"), function(x) {
  return(x[1])
}))
processed_pseudo <- sort(table(processed_pseudo), decreasing = TRUE)
processed_pseudo <- processed_pseudo[!is.na(names(processed_pseudo))]

write.table(processed_pseudo, file = "mHetGla1.Pri.processed.speduogenes.tsv", quote = FALSE,
            row.names = FALSE, col.names = FALSE, sep = "\t")

# processed_pseudogp <- gprofiler2::gost(names(processed_pseudo), organism = "mmusculus", evcodes = TRUE)
# (GO enrichment on processed pseudogene source genes - currently disabled)

# ---------------------------------------------------------------------
# Load a manually curated table of candidate NMR-specific gene copy
# events (used to spot-check gene duplications identified above)
# ---------------------------------------------------------------------
NMR_copy_pseudoCheck <- read.delim("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure4/Figure4A_Duplications/v20260516 NNMR gene copy events -sheet1.txt", header = T, as.is = T, sep = "\t")

# 8086 - new method
# 6355 - mikado only


# Split the semicolon-delimited list of HetGla gene IDs for the first
# 128 candidate gene-copy events into per-event ID vectors
sep <- strsplit(NMR_copy_pseudoCheck$HetGla.IDs, ";")[1:128]
unlist(strsplit(NMR_copy_pseudoCheck$HetGla.IDs, ";")[1:128])

# ---------------------------------------------------------------------
# Subset annotation_metrics to transcripts whose gene Name matches one
# of the curated candidate gene-copy IDs, and isolate the subset
# classified as pseudogenes among those copies
# ---------------------------------------------------------------------
annotation_metrics[annotation_metrics$Name %in% unlist(strsplit(NMR_copy_pseudoCheck$HetGla.IDs, ";")[1:128]), ]
annotation_metrics_copy <- annotation_metrics[annotation_metrics$Name %in% unlist(strsplit(NMR_copy_pseudoCheck$HetGla.IDs, ";")[1:128]), ]
duplicated_pseudos <- annotation_metrics_copy[annotation_metrics_copy$gene_type == "pseudogene", c("tid", "alias", "score", "type", "minScore", "gene_type", "Name")]

# Write out a fuller table of the gene-copy candidate transcripts
annotation_metrics_copy1 <- annotation_metrics_copy[, c("parent", "tid", "alias", "score", "type", "is_complete", "cdna_length", "exon_num", "minScore", "gene_type", "Name")]
write.table(annotation_metrics_copy1, file = "mHetGla1.1.pri.gene_copy_analysis.tsv",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# Reduce to one row per parent gene (i.e., gene-level view of the copy set)
annotation_metrics_copy_geneonly <- annotation_metrics_copy[!duplicated(annotation_metrics_copy$parent), ]




# ---------------------------------------------------------------------
# Extract pseudogene GFF records corresponding to the duplicated
# pseudogenes identified above, and export as a standalone GFF3
# ---------------------------------------------------------------------
pseudo_loss <- pseudo_gff[pseudo_gff$parent_gene %in% duplicated_pseudos$tid, ]
export(pseudo_loss, con = "mHetGla1.1.pri_genecopy_event_pseudogenes.gff", format = "GFF3")

# ---------------------------------------------------------------------
# Reload the gene-copy analysis table (written earlier) without a
# header, to use as a row-index lookup for extracting cDNA sequences
# ---------------------------------------------------------------------
annotation_metrics_copy1 <- read.table("mHetGla1.1.pri.gene_copy_analysis.tsv", header = FALSE, as.is = T, sep = "\t")

library(Biostrings)

# Load the full cDNA FASTA and clean sequence names to just the transcript ID
mHetGla1.1.pri.cdna <- readDNAStringSet("mHetGla1.1.pri.cdna.fa.gz")
names(mHetGla1.1.pri.cdna) <- unlist(lapply(strsplit(names(mHetGla1.1.pri.cdna), " "), function(x) return(x[1])))

# Subset the cDNA set to just the gene-copy transcripts (column V2 = tid),
# rename sequences as "<tid>-<gene_type>" (column V11), and write to FASTA
rownames(annotation_metrics_copy1) <- annotation_metrics_copy1$V2
mHetGla1.1.pri.cdna1 <- mHetGla1.1.pri.cdna[rownames(annotation_metrics_copy1)]
names(mHetGla1.1.pri.cdna1) <- paste0(rownames(annotation_metrics_copy1), "-", annotation_metrics_copy1$V11)
writeXStringSet(mHetGla1.1.pri.cdna1, file = "mHetGla1.1.pri.genecopies.cDNA.fa")


