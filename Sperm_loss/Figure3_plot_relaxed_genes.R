# ---------------------------------------------------------------------
# Replot FukAns + NMR
# Plot RELAX selection intensity (K values) per gene as a grouped bar chart
# ---------------------------------------------------------------------

# Load merged RELAX results across all tested selection schemes/lineages
relax <- read.delim("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure3/Figure3C_GeneLoss/relax_results_all_schemes_merged.tsv")


# ---------------------------------------------------------------------
# Main gene set used for plotting
# ---------------------------------------------------------------------
# Filter for genes that are not artifats
relax <- relax[relax]
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

# ---------------------------------------------------------------------
# Reshape from wide (one column per species) to long format
# (one row per gene-species-K value combination) for grouped bar plotting
# ---------------------------------------------------------------------
library(reshape2)
df_long <- melt(df)
colnames(df_long) <- c("Gene", "Species", "Value")


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
ggsave("plot.png", plot = p, dpi = 300)