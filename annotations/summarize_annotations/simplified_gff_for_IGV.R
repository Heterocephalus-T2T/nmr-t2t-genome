library(rtracklayer)

# ---------------------------------------------------------------------
# Build a BED12-style color-coded gene track for genome browser display,
# with colors indicating gene type (protein-coding / pseudogene / ncRNA)
# ---------------------------------------------------------------------

# Load per-transcript annotation metrics (includes pseudogenes)
metrics <- read.table("~/Desktop/NMR_multi_assembly/v090226/annotations/mHetGla1.pri_annotation_metrics_with_pseudogenes_v120426.tsv.gz",
                      header = TRUE, as.is = TRUE, sep = "\t")

# Load the MT/nuclear annotation GFF and restrict to gene-level features
# (protein-coding genes and ncRNA genes only)
gff <- readGFF("~/Desktop/NMR_multi_assembly/v090226/annotations/MT_annotations/mHetGla1.pri.MT.v020326.gff")
gff <- as.data.frame(gff[gff$type %in% c("gene", "ncRNA_gene"), ])

# Reduce transcript-level metrics to one row per gene (parent), keeping
# just the columns needed for the BED track (name, gene type, completeness)
metrics_gene <- metrics[!duplicated(metrics$parent), ]
metrics_gene <- metrics_gene[, c("Name", "parent", "gene_type", "is_complete")]

# Extract basic coordinate/strand info for each gene from the GFF
metrics_bed <- gff[, c("seqid", "start", "end", "strand", "ID")]

# Join gene coordinates with gene-level metrics on gene ID
merged_bed <- merge(metrics_bed, metrics_gene, by.x = "ID", by.y = "parent")

# BED score column (unused, set to 0) and feature length
merged_bed$score <- 0
merged_bed$length <- merged_bed$end - merged_bed$start

# ---------------------------------------------------------------------
# Assign RGB color per gene for browser display:
# - green (0,138,0)   = complete protein-coding gene (default)
# - orange (211,84,0) = pseudogene OR incomplete gene
# - lavender (216,191,216) = ncRNA gene
# ---------------------------------------------------------------------
rgb <- rep("0,138,0", nrow(merged_bed))
rgb[merged_bed$gene_type == "pseudogene" | merged_bed$is_complete == "False"] <- "211,84,0"
rgb[merged_bed$gene_type == "ncRNA_gene"] <- "216,191,216"
merged_bed$rgb <- rgb

# Assemble BED-style columns: chrom, start, end, name, score, strand,
# thickStart, thickEnd, itemRgb (using start/end again for thickStart/thickEnd
# since no separate CDS boundaries are being tracked here)
mHetGla1Pri_color_bed <- merged_bed[, c("seqid", "start", "end", "Name", "score", "strand", "start", "end", "rgb")]

# Write out the full color-coded BED (labels included for all genes,
# including pseudogenes/ncRNAs)
write.table(mHetGla1Pri_color_bed, file = "mHetGlaPri_gene_color.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# ---------------------------------------------------------------------
# Make a second version where only protein-coding gene names are shown
# (pseudogenes/ncRNAs get their name blanked to "." to reduce clutter)
# ---------------------------------------------------------------------
mHetGla1Pri_color_bed2 <- mHetGla1Pri_color_bed
mHetGla1Pri_color_bed2$Name[mHetGla1Pri_color_bed2$rgb != "0,138,0"] <- "."

# Clean up gene names: strip "mikado." prefix and shorten Ensembl-style
# NMR gene IDs (ENSHGLG000000...) to a shorter "Ens..." format
mHetGla1Pri_color_bed2$Name <- gsub("mikado.", "", mHetGla1Pri_color_bed2$Name)
mHetGla1Pri_color_bed2$Name <- gsub("ENSHGLG000000", "Ens",
                                    mHetGla1Pri_color_bed2$Name)

# Write BED with cleaned names, still retaining "-copyN" suffixes
# (i.e., distinguishing individual gene-copy duplicates by name)
write.table(mHetGla1Pri_color_bed2, file = "mHetGlaPri_PcGName_withcopy.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Strip "-copyN" suffixes entirely, collapsing gene-copy duplicates back
# to their base gene name, then write a final simplified BED
mHetGla1Pri_color_bed2$Name <- unlist(lapply(strsplit(mHetGla1Pri_color_bed2$Name, "-copy"), function(x) return(x[1])))
write.table(mHetGla1Pri_color_bed2, file = "mHetGlaPri_PcGName.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
