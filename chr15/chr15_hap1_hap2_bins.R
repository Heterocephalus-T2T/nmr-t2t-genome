library(dplyr)
library(readr)
library(tidyr)
library(purrr)
setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1B_synteny/")

#
## This script takes the chr15 hap1 segments and the hap1-hap2 alignments to identify the segment regions in hap2
#

# ---------------------------------------------------------------------
# Load haplotype 1 vs haplotype 2 whole-genome alignment (PAF format,
# restricted to chrX/chrY-containing scaffolds per filename) and keep
# only the standard 12 PAF columns (drop any extra tags)
# ---------------------------------------------------------------------
mHetGla1_h1_h2 <- read.table("mHetGla1_h1_vs_h2.chrXY.paf.gz", header = F, as.is = T, sep = "\t")
mHetGla1_h1_h2 <- mHetGla1_h1_h2[, 1:12]
colnames(mHetGla1_h1_h2) <- c(
  "qname",      # query sequence name
  "qlen",       # query sequence length
  "qstart",     # query start (0-based)
  "qend",       # query end (0-based, exclusive)
  "strand",     # "+" or "-"
  "tname",      # target sequence name
  "tlen",       # target sequence length
  "tstart",     # target start
  "tend",       # target end
  "nmatch",     # number of matching bases
  "alen",       # alignment block length
  "mapq"        # mapping quality
)

# Load manually curated chr15 hap1 segment boundaries (used to split
# chr15 into named blocks for cross-haplotype comparison)
seg_h1 <- read.delim("chr15_hap1_segment_locations.bed", header = F, as.is = T, sep = "\t")
colnames(seg_h1) <- c("chr", "start", "end", "name")

# Restrict alignments to chr15-vs-chr15 (both query and target), and to
# the region below 30 Mb where the segments of interest are located
mHetGla1_h1_h2 <- mHetGla1_h1_h2[grepl("-chr15", mHetGla1_h1_h2$qname) & grepl("-chr15", mHetGla1_h1_h2$tname) &
                                   mHetGla1_h1_h2$qend < 3e7 & mHetGla1_h1_h2$tend < 3e7, ]

write.table(mHetGla1_h1_h2, file = "chr15_hap1_hap2.paf",
            quote = FALSE,row.names=FALSE,col.names = FALSE,sep="\t")

# ---------------------------------------------------------------------
# For each hap1-defined segment, extract the alignments whose target
# (hap1) coordinates fall within that segment and have a substantial
# match length (>10kb), storing each subset in a named list
# ---------------------------------------------------------------------
segs <- list()
for (i in 1:length(seg_h1$name)) {
  
  start <- seg_h1$start[i]
  end <- seg_h1$end[i]
  name <- seg_h1$name[i]
  
  mHetGla1_h1_h2_seg <- mHetGla1_h1_h2[mHetGla1_h1_h2$tstart > start & mHetGla1_h1_h2$tend < end & mHetGla1_h1_h2$nmatch > 1e4, ]
  
  #plot(mHetGla1_h1_h2_seg$qstart,mHetGla1_h1_h2_seg$qend)
  
  #plot(mHetGla1_h1_h2_seg$tstart,mHetGla1_h1_h2_seg$tend)
  
  segs[[name]] <- mHetGla1_h1_h2_seg
}

# ---------------------------------------------------------------------
# Manually determine, for each hap1 segment, the corresponding hap2
# (query) coordinate range by inspecting each segment's alignments.
# This is an interactive/exploratory process, not a fully automated one.
# ---------------------------------------------------------------------

# segA: straightforward, use full min/max of query coordinates
hap2_segA <- c("chr15", min(segs$segA$qstart), max(segs$segA$qend), "segA")

## Manually check segB
segs$segB # mapQ is always 0 or 1 -- i.e. segB alignments are low-confidence/repetitive, so segB is not used further

# Manually check segC
plot(segs$segC$qstart, segs$segC$alen)

hap2_segA <- c("chr15", min(segs$segA$qstart), max(segs$segA$qend), "segA")

# SegC: Start > end of segA,
hap2_segC <- c("chr15", segs$segC$qstart, max(segs$segC$qend), "segC")

hap2_segC <- c("chr15", segs$segD$qstart, max(segs$segD$qend), "segD")

plot(segs$segE$qstart, segs$segE$qend)

# Restrict segD to alignments past 15 Mb, for closer inspection
segs$segD[segs$segD$qstart > 1.5e7, ]
plot(segs$segD$qstart, segs$segD$qend)

# ---------------------------------------------------------------------
# Clean up sequence names in the full PAF table for JBrowse display:
# strip assembly-version prefixes from query (alt/hap2) and target
# (pri/hap1) sequence names
#
# ---------------------------------------------------------------------
head(mHetGla1_h1_h2)
mHetGla1_h1_h2$V1 <- gsub("v090226.alt-", "", mHetGla1_h1_h2$qname)
mHetGla1_h1_h2$V6 <- gsub("v090226.pri-", "", mHetGla1_h1_h2$tname)

# Write out the full cleaned PAF file for JBrowse visualization
write.table(mHetGla1_h1_h2, file = "mHetGla1_h1_vs_h2.chrXY.jBrowse.paf",
            quote = FALSE, row.names = FALSE, col.names = FALSE)


mHetGla1_h1_h2_chr15 <- mHetGla1_h1_h2[mHetGla1_h1_h2$V1 == "chr15" & mHetGla1_h1_h2$V6 == "chr15", ]
write.table(mHetGla1_h1_h2_chr15[, 1:12], file = "mHetGla1_h1_vs_h2.chr15.jBrowse.paf",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# seg C/E hap2: chr15:3,795,013-19123885

# ---------------------------------------------------------------------
# Refine segD: restrict to alignments past 15 Mb with mapq > 1 (i.e.,
# reasonably confident mappings), and report the resulting coordinate span
# ---------------------------------------------------------------------
segD_match <- segs$segD[segs$segD$qstart > 1.5e7 & segs$segD$mapq > 1, ]
min(segD_match$qstart)
max(segD_match$qend)

# ---------------------------------------------------------------------
# Final hap2 segment coordinates, determined via the manual inspection
# above. segCE and segD use hardcoded coordinates derived from the
# preceding interactive exploration rather than being recomputed live.
# ---------------------------------------------------------------------
library(dplyr)
library(GenomicRanges)

## ---- 1. Read PAF -------------------------------------------------
paf_file <- "/Users/dsokolowski/Library/CloudStorage/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1B_synteny/map_chr15_hap2_segCE_to_mHetGla1.1.pri.segE.paf"

paf_cols <- c("qname","qlen","qstart","qend","strand",
              "tname","tlen","tstart","tend",
              "nmatch","alnlen","mapq")

raw <- read.table(paf_file, sep = "\t", header = FALSE,
                  stringsAsFactors = FALSE, fill = TRUE,
                  col.names = paste0("V", 1:30), quote = "")

paf <- raw[, 1:12]
names(paf) <- paf_cols

tag_cols <- raw[, 13:ncol(raw)]
paf$tags <- apply(tag_cols, 1, function(x) paste(na.omit(x[x != ""]), collapse = ";"))
paf$dv <- as.numeric(sub(".*dv:f:([0-9.eE+-]+).*", "\\1", paf$tags))

## ---- 2. Keep only near-full-length hits (the "green" alignments) ----
# target (mHetGla1.1 segE) is 48270 bp; the big bars cover almost all of it
paf_full <- paf %>% filter(alnlen >= 40000)

## ---- 3. Cluster by query position, merging across gaps --------------
gr <- GRanges(seqnames = paf_full$qname,
              ranges   = IRanges(start = paf_full$qstart + 1,
                                 end   = paf_full$qend))

# Increase min.gapwidth until you land on 5 clusters.
# Start big and back off if it over-merges.
gap_size <- 500000   # <-- "net size": tune this

merged <- reduce(gr, min.gapwidth = gap_size)
merged_df <- as.data.frame(merged)
print(merged_df)
cat("Number of blocks:", length(merged), "\n")

gap_sizes <- c(1e4, 5e4, 1e5, 2e5, 3e5, 5e5, 7e5, 1e6)
sapply(gap_sizes, function(g) length(reduce(gr, min.gapwidth = g)))

paf_full$block <- findOverlaps(gr, merged, select = "first")

block_summary <- paf_full %>%
  group_by(block) %>%
  summarise(
    qstart_min = min(qstart), qend_max = max(qend),
    n_alignments = n(),
    mean_dv = mean(dv, na.rm = TRUE),
    tstart_min = min(tstart), tend_max = max(tend)
  )

block_summary

merged_df$start <- merged_df$start + 4012549
merged_df$end <- merged_df$end + 4012549

merged_df$seqnames <- "chr15"

write.table(merged_df[,1:3], file = "~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1B_synteny/chr15_hap2_segE.bed",
            quote=FALSE,row.names = FALSE,col.names = FALSE,
            sep="\t")

# bedtools subtract -a chr15_hap2_segCE.bed -b chr15_hap2_segE.bed > chr15_hap2_segC.bed

setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1B_synteny")

locations_CE <- read.table("chr15_hap2_segment_locations.bed",header=F,as.is=T,sep="\t")
segC <- read.table("chr15_hap2_segC.bed",header=F,as.is=T,sep="\t")
segE <- read.table("chr15_hap2_segE.bed",header=F,as.is=T,sep="\t")
segC$V4 <- "segC"
segE$V4 <- "segE"
chr15_hap2_segment_locations <- as.data.frame(rbind(locations_CE[c(1,3,4),],segC,segE))
chr15_hap2_segment_locations <- chr15_hap2_segment_locations[order(chr15_hap2_segment_locations$V2),]

write.table(chr15_hap2_segment_locations,
            file = "chr15_hap2_segment_locations_v260916.bed",
            quote=F,row.names=F,col.names=F,sep="\t")           
