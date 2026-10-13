library(stringr)
setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Supplementary Figures/OR_Ligand_analysis/")

# ---------------------------------------------------------------------
# Load M2OR (molecule-to-olfactory-receptor) response summary data
# ---------------------------------------------------------------------
M2OR <- read.table("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Supplementary Figures/OR_Ligand_analysis/NMR_M2OR_summaries.txt",
                   header = T, as.is = T, sep = "\t")

# Standardize ligand names to lowercase for consistent matching
M2OR$Molecule.Name <- tolower(M2OR$Molecule.Name)
head(M2OR)

# Build a unique OR-ligand pair identifier (receptor-ligand combination)
LR_pair <- paste0(M2OR$OR, "-", M2OR$Molecule.Name)
M2OR$pair <- LR_pair

# Remove duplicate OR-ligand pairs, keeping first occurrence
M2OR_paired <- M2OR[!duplicated(LR_pair), ]

# Keep only pairs where the receptor showed a positive response
M2OR_responsive <- M2OR_paired[M2OR_paired$Responsive > 0, ]

# Standardize OR gene naming convention (e.g., "OR1A1" -> "Or1a1" style)
M2OR_responsive$OR <- stringr::str_to_title(M2OR_responsive$OR)

# Manual correction: reassign clade for Or6c to "Mus"
M2OR_responsive$Clade[M2OR_responsive$OR == "Or6c"] <- "Mus"

# Get the unique sets of responsive ORs and their associated ligands
ORS <- unique(M2OR_responsive$OR)
ligands <- unique(M2OR_responsive$Molecule.Name)

# ---------------------------------------------------------------------
# Load OR gene copy number data across species and subset to the
# responsive ORs identified above
# ---------------------------------------------------------------------
setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1B_synteny/")
OR <- read.csv("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure1/Figure1B_synteny/OR gene copy number.csv",
               header = TRUE, as.is = TRUE)
colnames(OR) <- c("name", "mouse", "hsap", "rat", "DeerMouse", "NMR", "CDMR", "mmus_minus_NMR")
rownames(OR) <- OR$name
OR1 <- OR[, 2:7]

# Subset copy-number table to only the ORs found responsive in M2OR data
OR1_in <- OR1[rownames(OR1) %in% ORS, ]


# ---------------------------------------------------------------------
# Load ligand summary/annotation table (functional group, odor tone,
# volatility, etc.) and restrict it to ligands relevant to this analysis
# ---------------------------------------------------------------------
Ligand_summary <- read.table("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Supplementary Figures/OR_Ligand_analysis/Ligand_Summary2_20_Sarah_20-_20Ligand_Summary2_20.tsv",
                             header = TRUE, as.is = TRUE, sep = "\t")
rownames(Ligand_summary) <- Ligand_summary$Ligand

# Keep only ligands present in both the summary table and the M2OR data
ligand2 <- intersect(Ligand_summary$Ligand, ligands)
Ligand_summary <- Ligand_summary[ligand2, ]

# ---------------------------------------------------------------------
# Build a one-row-per-OR clade lookup table (deduplicated by OR-clade
# combination) for use as heatmap column annotation later
# ---------------------------------------------------------------------
OR_clade <- paste0(M2OR_responsive$OR, "-", M2OR_responsive$Clade)
OR_clade_df <- M2OR_responsive[!duplicated(OR_clade), ]

# ---------------------------------------------------------------------
# Build a binary OR x ligand response matrix: 1 = OR responds to ligand
# ---------------------------------------------------------------------
mat <- matrix(0, nrow = length(ORS), length(ligand2))
rownames(mat) <- ORS
colnames(mat) <- ligand2

for (i in ORS) {
  for (j in ligand2) {
    pair <- paste0(i, "-", j)
    if (pair %in% M2OR_responsive$pair) mat[i, j] <- 1
  }
}

# Transpose so ligands are rows and ORs are columns (for heatmap layout)
mat <- t(mat)

# ---------------------------------------------------------------------
# Prepare annotation data frames for the heatmap: OR clade (columns)
# and ligand functional group / odor descriptors (rows)
# ---------------------------------------------------------------------
annotation_col <- as.data.frame(OR_clade_df$Clade)
colnames(annotation_col) <- "clade"

annotation_row <- Ligand_summary[, c(2:3, 6)]
rownames(annotation_row) <- rownames(Ligand_summary)
rownames(annotation_col) <- OR_clade_df$OR

# Reorder columns (ORs) by clade grouping
group_order <- order(annotation_col$clade)
mat2 <- mat[, group_order]
annotation_col2 <- annotation_col[group_order, , drop = FALSE]

# Reorder rows (ligands) by primary functional group
group_order <- order(annotation_row$Primary_functional_group)
mat3 <- mat2[group_order, ]
annotation_row2 <- annotation_row[group_order, , drop = FALSE]

# Binary color scale: white = no response, red = responsive
my_colors <- c("white", "red")

# ---------------------------------------------------------------------
# Heatmap 1: full OR x ligand response matrix, annotated by clade
# (columns) and functional group/odor tone (rows)
# ---------------------------------------------------------------------
pdf(file = "M2OR_ligand_analysis_scent.pdf", width = 16, height = 22)
print(pheatmap::pheatmap(mat3, color = my_colors,
                         annotation_col = annotation_col2, annotation_row = annotation_row2,
                         cluster_cols = FALSE, cluster_rows = FALSE))
dev.off()

# ---------------------------------------------------------------------
# Collapse the response matrix by summing responses within each clade,
# producing a ligand x clade summary matrix
# ---------------------------------------------------------------------
clades <- unique(annotation_col2$clade)
eachClade <- c()
for (i in clades) {
  toSum <- mat3[, rownames(annotation_col2)[annotation_col2$clade == i]]
  eachClade <- cbind(eachClade, rowSums(toSum))
}
colnames(eachClade) <- clades

# Heatmap 2: ligand x clade summed response counts (row-scaled)
print(pheatmap::pheatmap(eachClade, annotation_row = annotation_row2,
                         cluster_cols = FALSE, cluster_rows = FALSE, scale = "row",
                         color = my_colors, cellwidth = 10))

# ---------------------------------------------------------------------
# Further collapse by ligand functional group, summing clade responses
# within each functional group
# ---------------------------------------------------------------------
Groups <- unique(annotation_row2$Primary_functional_group)
eachGroup <- c()
for (i in Groups) {
  toSum <- eachClade[rownames(annotation_row2)[annotation_row2$Primary_functional_group == i], ]
  if (is.null(nrow(toSum))) {
    # Only one ligand in this group -> toSum is a vector, not a matrix
    eachGroup <- rbind(eachGroup, toSum)
  } else {
    eachGroup <- rbind(eachGroup, colSums(toSum))
  }
}
rownames(eachGroup) <- Groups

# Heatmap 3: functional group x clade summed response counts,
# with values displayed as numbers on the tiles
pdf(file = "Clade_by_Group.pdf")
print(pheatmap::pheatmap(eachGroup, scale = "row", display_numbers = eachGroup, color = my_colors,
                         cellwidth = 20, cluster_rows = FALSE, cluster_cols = FALSE))
dev.off()

# ---------------------------------------------------------------------
# Same collapsing approach, but grouped by odor tone instead of
# functional group
# ---------------------------------------------------------------------
Smells <- unique(annotation_row2$Odor_tone)
eachSmell <- c()
for (i in Smells) {
  toSum <- eachClade[rownames(annotation_row2)[annotation_row2$Odor_tone == i], ]
  if (is.null(nrow(toSum))) {
    eachSmell <- rbind(eachSmell, toSum)
  } else {
    eachSmell <- rbind(eachSmell, colSums(toSum))
  }
}
rownames(eachSmell) <- Smells

# Heatmap 4: odor tone x clade summed response counts
pdf(file = "clade_by_odour.pdf")
print(pheatmap::pheatmap(eachSmell, scale = "row", display_numbers = eachSmell,
                         color = my_colors, cellwidth = 20, cluster_rows = FALSE, cluster_cols = FALSE))
dev.off()

# ---------------------------------------------------------------------
# Identify ligands for which the 3rd clade column
# has the highest summed response count, and save that ligand list
# ---------------------------------------------------------------------
NMR_top <- rownames(eachClade)[apply(eachClade, 1, which.max) == 3]
write.table(NMR_top, file = "NMR_top_ligands.tsv", quote = F, row.names = F,
            col.names = F, sep = "\t")