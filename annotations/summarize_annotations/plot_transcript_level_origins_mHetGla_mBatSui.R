setwd("~/Desktop/NMR_multi_assembly/v090226/annotations/MT_annotations/")
library(ggplot2)
library(rtracklayer)
library(pheatmap)

# ---------------------------------------------------------------------
# For each species/assembly (naked mole-rat and Cape dune mole-rat), load the gene
# annotation metrics table and plot a heatmap of gene counts by
# annotation source ("origin") x gene type
# ---------------------------------------------------------------------
for (i in c("mHetGla1.pri", "BatSuiPri")) {
  
  # Select the correct input file depending on species
  if (i == "mHetGla1.pri") {
    file <- "~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure4/Figure4BC_Pms1/mHetGla1.pri_annotation_metrics_with_pseudogenes_v120426.tsv.gz"
  } else {
    file <- "~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure4/Figure4BC_Pms1/BatSuiPri_annotation_metrics_with_pseudogenes_v20260605.tsv"
  }
  
  # Load annotation metrics for this species
  mhetgla_genes <- read.delim(file)
  
  # Extract the "origin" (annotation source/tool) from the alias field,
  # assuming it's encoded as the first underscore-delimited token
  mhetgla_genes$origin <- unlist(lapply(strsplit(mhetgla_genes$alias, "_"), function(x) {
    return(x[1])
  }))
  
  # Cross-tabulate: number of genes per origin x gene type combination
  gene_ratio <- table(mhetgla_genes$origin, mhetgla_genes$type)
  
  # 1. Create the palette function
  white_to_red <- colorRampPalette(c("white", "orange"))
  # 2. Generate 100 colors along the gradient
  my_colors <- white_to_red(100)
  
  # Build heatmap: origin x type counts, with counts displayed on tiles,
  # no clustering (preserve natural row/column order)
  p1 <- pheatmap(gene_ratio, scale = "none", color = my_colors, display_numbers = gene_ratio,
                 cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 12, cellwidth = 50, cellheight = 50)
  
  # Save heatmap to a species-specific PDF
  pdf(file = paste0("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Supplementary Figures/Transcript_Summaries/", i, "_transcript_by_source.pdf"))
  print(p1)
  dev.off()
}

# ---------------------------------------------------------------------
# Summarize StringTie-assembled transcripts by tissue of origin
# ---------------------------------------------------------------------

# Load the merged StringTie GFF (per-tissue transcript assemblies merged)
Transcript_by_tissue <- readGFF("~/Desktop/NMR_multi_assembly/v090226/annotations/mHetGla1.pri.stringtie.s.gff.gz")

# Keep only "transcript"-level features (drop exon/other sub-features)
Transcript_by_tissue <- Transcript_by_tissue[Transcript_by_tissue$type == "transcript", ]

# Extract tissue of origin from the transcript ID, assuming StringTie's
# merged ID format encodes the source tissue before ".merged."
Transcript_by_tissue$origin <- unlist(lapply(strsplit(Transcript_by_tissue$ID, ".merged."), function(x) return(x[1])))

# Count transcripts per tissue
transcriptTab <- as.data.frame(table(Transcript_by_tissue$origin))
colnames(transcriptTab) <- c("tissue", "transcripts")

# Quick unsorted bar plot of transcript counts by tissue
ggplot(transcriptTab, aes(x = transcripts, y = tissue)) +
  geom_bar(stat = "identity") + theme_classic()

# Final plot: horizontal bar chart, tissues sorted by transcript count
p1 <- ggplot(transcriptTab, aes(x = reorder(tissue, transcripts), y = transcripts)) +
  geom_col() + coord_flip() + theme_classic()

# Save the sorted tissue-transcript plot to PDF
pdf(file = "~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Supplementary Figures/Transcript_Summaries/mHetGla1.pri_transcript_by_tissue.pdf")
print(p1)
dev.off()