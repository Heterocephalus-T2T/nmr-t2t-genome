# ==========================================================================
# Purpose: Map a given haplotype-specific contig back to its constituent
#          diploid unitigs, using shared read (A-line) records across GFA
#          assembly graphs from hifiasm-style output.
# ==========================================================================

setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1A_overview/gfa/")

# Load the GFA files (no-sequence versions) as tab-delimited tables.
# Each file's rows include "A" lines, which record which reads support
# each contig/unitig (used later to link haplotype contigs to dip unitigs).
h1 <- read.delim("Pup3-NMR.asm.dip.hap1.p_ctg.noseq.gfa", header = FALSE, as.is = T, sep = "\t")   # haplotype 1 primary contigs
h2 <- read.delim("Pup3-NMR.asm.dip.hap2.p_ctg.noseq.gfa", header = FALSE, as.is = T, sep = "\t")   # haplotype 2 primary contigs
dip <- read.delim("Pup3-NMR.asm.dip.p_utg.noseq.gfa", header = FALSE, as.is = T, sep = "\t")       # diploid (unphased) unitig graph

# Keep only "A" (read alignment) lines from the diploid graph, and retain
# just the unitig name (col 2) and read name (col 5) columns.
dip <- dip[dip$V1 == "A", c(2, 5)]
head(dip)
colnames(dip) <- c("unitig", "reads")

# Store both haplotype contig graphs in a named list for easy lookup
# by the get_unitigs() function below.
contig_graphs <- list(h1 = h1,
                      h2 = h2)

# ---------------------------------------------------------------------
# get_unitigs()
# For a given haplotype ("h1" or "h2") and contig name, find all the
# diploid (dip) unitigs that share reads with that haplotype contig.
# Returns a table of unitig names -> count of shared reads, which is
# useful for identifying which dip unitigs make up / overlap a given
# phased contig.
# ---------------------------------------------------------------------
get_unitigs <- function(haplotype = "h2", contig = "h2tg000476l") {
  
  # Pull out the requested haplotype's contig graph
  hap <- contig_graphs[[haplotype]]
  
  # Filter to "A" lines only, keep contig name (col 2) and read name (col 5)
  hap <- hap[hap$V1 == "A", c(2, 5)]
  colnames(hap) <- c("contig", "reads")
  
  # Get the set of reads belonging to the contig of interest
  contig_reads <- hap$reads[hap$contig == contig]
  
  # Find which diploid unitigs contain those same reads
  dip_contig <- dip[dip$reads %in% contig_reads, ]
  
  # Return a frequency table of matching dip unitigs
  return(table(dip_contig$unitig))
  
  
}

# Example usage: investigate the unitig composition of a known
# problematic contig (h1tg000001l) in haplotype 1
unis <- get_unitigs("h1", "h1tg000001l")

# Print the unitig names as a single comma-separated string (for quick viewing)
paste0(names(unis), collapse = ",")

# Write the list of matching unitig names to file, one per line
write.table(names(unis), file = "unis.txt", quote = F, row.names = F, col.names = F, sep = "\t")

# Print unitig names to console without quotes, for a quick sanity check
print(names(unis), quote = F)