# ---------------------------------------------------------------------
# Replot FukAns + NMR
# Plot RELAX selection intensity (K values) per gene as a grouped bar chart
# ---------------------------------------------------------------------

setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure3/Figure3EF_relax_analysis/")

# Load merged RELAX results across all tested selection schemes/lineages
relax_old <- read.delim("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure3/Figure3EF_relax_analysis/relax_all_schemes_merged.tsv")
relax <- read.delim("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure3/Figure3EF_relax_analysis/sperm_summary_fdr.tsv")
gene_all <- read.delim("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure3/Figure3EF_relax_analysis/sperm_shared_genes.txt",header=F)



# ---------------------------------------------------------------------
# Main gene set used for plotting
# ---------------------------------------------------------------------
# Filter for genes that are not artifats
# relax <- relax[relax$Category_FDR %in% c("INTENSIFIED", "NOT_SIGNIFICANT","RELAXED"),]
relax <- relax[relax$Gene %in% gene_all$V1,]

relax_nmr <- relax[relax$Scheme == "nmr_vs_nonmolerat",]
relax_fans <- relax[relax$Scheme == "amr_vs_nonmolerat",]
relax_bsui <- relax[relax$Scheme == "bats_vs_nonmolerat",]

rownames(relax_nmr) <- relax_nmr$Gene
rownames(relax_bsui) <- relax_bsui$Gene
rownames(relax_fans) <- relax_fans$Gene

relax_nmr <- relax_nmr[gene_all$V1,]
relax_bsui <- relax_bsui[gene_all$V1,]
relax_fans <- relax_fans[gene_all$V1,]

relax_nmr$K_nmr_only <- relax_nmr$K
relax_nmr$Category_nmr_only <- relax_nmr$Category_FDR

relax_fans$K_amr_only <- relax_fans$K
relax_fans$Category_amr_only <- relax_fans$Category_FDR

relax_bsui$K_bats_only <- relax_bsui$K
relax_bsui$Category_bats_only <- relax_bsui$Category_FDR

relax_simple <- as.data.frame(cbind(relax_nmr[,c("K_nmr_only","Category_nmr_only")],
                                    relax_bsui[,c("K_bats_only","Category_bats_only")],
                                    relax_fans[,c("K_amr_only","Category_amr_only")]))

relax_simple$Gene <- relax_nmr$Gene

relax <- relax_simple

relax_union <- relax[relax$Category_nmr_only == "RELAXED" & relax$Category_amr_only %in% c("NOT_SIGNIFICANT", "INTENSIFIED") & relax$Category_bats_only %in%  c("NOT_SIGNIFICANT", "INTENSIFIED"), ]

# Standardize gene name capitalization (e.g., "GENE1" -> "Gene1") and
# use gene names as row identifiers
relax_union$Gene <- stringr::str_to_title(relax_union$Gene)
rownames(relax_union) <- relax_union$Gene

# Keep only the K (selection intensity parameter) values per lineage,
# rename to friendly species labels, and treat missing K as 1 (i.e.,
# no evidence of intensified/relaxed selection - neutral baseline)
relax_union <- relax_union[, c("K_nmr_only", "K_amr_only", "K_bats_only")]
colnames(relax_union) <- c("NMR", "AnMR", "CDMR")
relax_union[is.na(relax_union)] <- 1

# Restore Gene as an explicit column (needed for reshaping to long format)
relax_union$Gene <- rownames(relax_union)
df <- relax_union



##
# Filter for top 15 genes based on k-value for NMR

k_min <- apply(df[,c("AnMR","CDMR")],1,min)

k_diff <- k_min-df$NMR



df <- df[df$Gene %in% names(sort(k_diff,decreasing = T)[1:20]),]

##

# ---------------------------------------------------------------------
# Reshape from wide (one column per species) to long format
# (one row per gene-species-K value combination) for grouped bar plotting
# ---------------------------------------------------------------------
library(reshape2)
df_long <- melt(df)
colnames(df_long) <- c("Gene", "Species", "Value")

df_long$Value[df_long$Value > 2] <- 2

# Match colors from the Venn diagram
species_colors <- c(
  "K_nmr_only"  = "#F27C59",  # orange/red
  "K_anmr_only" = "#6FA8DC",  # blue
  "K_cdmr_only" = "#4DB6AC"   # teal
)


# ---------------------------------------------------------------------
# Plot: grouped bar chart of K value per gene, one bar per species/lineage
# ---------------------------------------------------------------------
p <- ggplot(df_long, aes(x = Gene, y = Value, fill = Species)) +
  #coord_flip() +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  scale_fill_manual(
    values = c(
      "NMR" = "#F27C59",
      "AnMR" = "#6FA8DC",
      "CDMR" = "#4DB6AC"
    )
  ) +
  labs(
    x = NULL,
    y = "K value",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "top"
  )
ggsave("NMR_relax_final_top20.png", plot = p, dpi = 300)
ggsave("NMR_relax_final_top20.pdf", plot = p, dpi = 300)


forVenn <- list(nmr=relax_nmr$Gene[relax_nmr$Category_nmr_only == "RELAXED"],
                bsui=relax_bsui$Gene[relax_bsui$Category_bats_only == "RELAXED"],
                amr=relax_fans$Gene[relax_fans$Category_amr_only == "RELAXED"])


# install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)

ggVennDiagram(forVenn, 
              category.names = c("NMR", "B. sui", "F. anselli"),
              label = "count") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "none")

#
## Make updated GFF for RNAseq
# 

relax <- read.delim("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure3/Figure3EF_relax_analysis/sperm_summary_fdr.tsv")
relax2 <- relax[relax$Scheme %in% c("nmr_vs_nonmolerat","allmolerat_vs_nonmolerat"),]
relax2 <- relax2[relax2$Category_FDR %in% c("RELAXED","INTENSIFIED","NOT_SIGNIFICANT"),]

gene_all <- read.delim("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure3/Figure3EF_relax_analysis/sperm_shared_genes.txt",header=F)
relax2 <- relax2[relax2$Gene %in% gene_all$V1,]

relax2$category_both <- paste0(relax2$Scheme,"-",relax2$Category_FDR)

relaxed_gmt <- list()
uni <- unique(relax2$category_both)
library(stringr)
relax2$Gene_lower <- stringr::str_to_title(relax2$Gene)


for(i in uni) {
  relaxed_gmt[[i]]$name <- relaxed_gmt[[i]]$id <- i
  relaxed_gmt[[i]]$genes <- relax2$Gene_lower[relax2$category_both == i]
  
}

relaxed_gmt[["total"]]$name <- relaxed_gmt[["total"]]$id <- "Total"
relaxed_gmt[["total"]]$genes <- relax_analysis_in$Gene_mmus


for (i in uni) {
  for (j in status) {
    relaxed_gmt[[paste0(i, "_", j)]]$genes <-
      relax_analysis_in$Gene_mmus[relax_analysis_in$Category == j & relax_analysis_in$Scheme == i]
    relaxed_gmt[[paste0(i, "_", j)]]$id <-
      relaxed_gmt[[paste0(i, "_", j)]]$name <- paste0(i, "_", j)
  }
}
# Add a "total" background set covering all RELAX-tested genes
relaxed_gmt[["total"]]$genes <- unique(relax2$Gene_lower)
relaxed_gmt[["total"]]$name <- relaxed_gmt[["total"]]$id <- "Total"
class(relaxed_gmt) <- "GMT"
library(ActivePathways)
is.GMT(relaxed_gmt)
saveRDS(relaxed_gmt, file = "relaxed_sperm_genes_gmt.rds")
write.GMT(relaxed_gmt, file = "relaxed_sperm_genes_Sept10_26.gmt")

# NMR vs F.ans
fisher.test(matrix(c(83,385-83,
              44,385-44),nrow = 2,ncol = 2))

# NMR vs B.sui
fisher.test(matrix(c(83,385-83,
                     21,385-21),nrow = 2,ncol = 2))
# F.ans vs B.sui
fisher.test(matrix(c(44,385-44,
                     21 ,385-21 ),nrow = 2,ncol = 2))


            