library(rtracklayer)

# ---------------------------------------------------------------------
# File paths: per-transcript annotation metrics table (includes
# pseudogenes) and the full genome GFF (includes MT annotations)
# ---------------------------------------------------------------------
gene <- "~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/mHetGla1.pri_annotation_metrics_with_pseudogenes_v120426.tsv.gz"
gff <- readGFF("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/mHetGla1.pri.MT.v020326.gff.gz")

# ---------------------------------------------------------------------
# Load gene-level annotation metrics and reduce to one row per gene
# (parent), since the table originally has one row per transcript
# ---------------------------------------------------------------------
gene_list <- read.table(gene, header = T, as.is = T, sep = "\t")
head(gene_list)
gene_list <- gene_list[!duplicated(gene_list$parent), ]
table(gene_list$type)

# ---------------------------------------------------------------------
# Extract non-coding RNA annotations predicted by cmscan (Rfam/Infernal)
# or MirMachine (miRNA prediction), keeping only feature types with
# "RNA" in the name (e.g., tRNA, snoRNA, miRNA, etc.)
# ---------------------------------------------------------------------
cmscan <- gff[gff$source == "cmscan" | gff$source == "MirMachine", ]
cmscan <- cmscan[grep("RNA", cmscan$type), ]
cmscan_tab <- cmscan[, c("ID", "type")]

# Simplify the gene-level metrics table to just gene ID and biotype/gene_type
gene_list_simple <- gene_list[, c("parent", "gene_type")]

# ---------------------------------------------------------------------
# Extract mitochondrial gene annotations from MitoFinder, excluding
# CDS features so only gene/exon-level records remain
# ---------------------------------------------------------------------
mito <- gff[gff$source == "mitofinder", ]
mito <- mito[mito$type != "CDS", ]
mito_simple <- mito[, c("ID", "type")]

# ---------------------------------------------------------------------
# Standardize column names across all three sources (mitochondrial,
# ncRNA, nuclear gene metrics) so they can be row-bound into one table
# ---------------------------------------------------------------------
colnames(mito_simple) <- colnames(gene_list_simple) <- colnames(cmscan_tab) <- c("geneID", "biotype")

# Combine mitochondrial, ncRNA (cmscan/MirMachine), and nuclear gene
# biotype tables into a single genome-wide gene ID -> biotype table
biotype_total <- as.data.frame(rbind(mito_simple, cmscan_tab, gene_list_simple), stringsAsFactors = FALSE)

# Write out the combined biotype table for downstream use (e.g., figure generation)
write.table(biotype_total, file = "~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1E_annotations/mHetGla1pri_geneID_by_biotype.txt",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")