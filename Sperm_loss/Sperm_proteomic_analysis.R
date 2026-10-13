###############################################################################
# Sperm proteomics: mouse vs. naked mole-rat (NMR)
#
# Purpose
#   1. Run GO enrichment (g:Profiler) on all proteins detected in mouse and
#      NMR sperm, and plot the top terms for each species.
#   2. Test one-to-one orthologous proteins for differential abundance
#      between species (per-protein Welch t-test on quantity ratios, FDR).
#   3. Annotate DE proteins with NMR gene type (gene vs. pseudogene) from the
#      mHetGla1 T2T annotation, and draw volcano plots.
#   4. Run GO enrichment on mouse-biased and NMR-biased proteins.
#   5. Compare proteomic fold changes with RNA-seq fold changes from the
#      mouse vs. NMR differential expression analysis.
#
# Inputs
#   - SCPC087-01_mouse_sperm_globalproteomics_raw__top_protein_quantities.csv
#       Mouse sperm protein quantities; column "Gene" plus 3 replicate
#       quantity columns (columns 3:5) and "Ratio" columns.
#   - sub_bre_sperm_proteomics.csv
#       NMR sperm protein quantities; column "Gene" (mouse orthologue
#       symbols) plus replicate columns from column 3 onward.
#   - mHetGla1.pri_annotation_metrics_with_pseudogenes_v120426.tsv.gz
#       NMR gene annotation with columns Name, parent, score, gene_type.
#   - mouse_vs_nmr_DE_allsamples.rds
#       RNA-seq DE table with Row.names, logFC, significance.
#
# Outputs
#   - NMR_mmus_one2one.csv          merged ratio table for shared proteins
#   - mouse_biased_proteomics.pdf   top GO terms, mouse-higher proteins
#   - NMR_biased_proteomics.pdf     top GO terms, NMR-higher proteins
#   - Plot objects gp2-gp6, p1 and two volcano plots (printed, not saved)
#
# Fold-change convention: log2FC = log2(mouse / NMR)
#   positive = higher in mouse, negative = higher in NMR
###############################################################################

library(gprofiler2)
library(ggplot2)
library(dplyr)
  setwd("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure3/Figure3_KLM_mouse_vs_NMR/")


# =============================================================================
# 1. Load and order the proteomics tables
# =============================================================================

# ---- Mouse sperm ----
mouse_sperm <- read.csv("SCPC087-01_mouse_sperm_globalproteomics_raw__top_protein_quantities.csv",header=T,as.is=T)

# Drop proteins with any missing value
mouse_sperm <- mouse_sperm[complete.cases(mouse_sperm),]

# Sort by mean abundance across the 3 replicates (columns 3:5), highest first
mouse_sperm <- mouse_sperm[order(rowMeans(mouse_sperm[,3:5]),decreasing = T),]

# Named vector of mean abundance per gene
mouse_means <- rowMeans(mouse_sperm[,3:5])
names(mouse_means) <- mouse_sperm$Gene

# ---- NMR sperm ----
sub_bre_proteomics_counts <- read.csv("sub_bre_sperm_proteomics.csv",header=T,as.is=T)
sub_bre_proteomics_counts <- sub_bre_proteomics_counts[complete.cases(sub_bre_proteomics_counts),]

# Sort by mean abundance across all replicate columns (3 to end), highest first
sub_bre_proteomics_counts <- sub_bre_proteomics_counts[order(rowMeans(sub_bre_proteomics_counts[3:ncol(sub_bre_proteomics_counts)]),decreasing = T),]


# =============================================================================
# 2. GO enrichment of all detected sperm proteins, per species
# =============================================================================
# Both queries use mouse gene symbols and the mouse annotation. No custom
# background is given, so g:Profiler tests against all annotated mouse genes
# (hence "nobg").

mouse <- gost(mouse_sperm$Gene,organism = "mmusculus",evcodes = T)
sperm_protein_nobg <- gost(sub_bre_proteomics_counts$Gene,organism = "mmusculus",evcodes = T)

# Keep GO terms only (GO:BP, GO:MF, GO:CC) and drop very broad terms
# (more than 3,000 genes)
sperm_protein_nobg1 <- sperm_protein_nobg$result[sperm_protein_nobg$result$term_size < 3e3 & grepl("GO",sperm_protein_nobg$result$source ),]
mouse1 <- mouse$result[mouse$result$term_size < 3e3 & grepl("GO",mouse$result$source ),]

# Replace the "source" column with the species label used for plot colours
mouse1$source <- "mouse"
sperm_protein_nobg1$source <- "NMR"

# ---- Top 5 GO terms: all mouse sperm proteins ----
gp2 <- ggplot(mouse1 %>% head(5),
              aes(x = reorder(term_name, -p_value),
                  y = -log10(p_value),color=source)) +
  geom_col() +
  coord_flip() +
  theme_bw() + theme(legend.position = "none") + 
  labs(x = "", y = "-log10(FDR)")+ scale_color_manual(values=c(
    "mouse"="purple",
    "NMR" = "gold"
  ))

# ---- Top 5 GO terms: all NMR sperm proteins ----
gp3 <- ggplot(sperm_protein_nobg1 %>% head(5),
              aes(x = reorder(term_name, -p_value),
                  y = -log10(p_value),color=source)) +
  geom_col() +
  coord_flip() +
  theme_bw() + theme(legend.position = "none") + 
  labs(x = "", y = "-log10(FDR)")+ scale_color_manual(values=c(
    "mouse"="purple",
    "NMR" = "gold"
  ))

# ---- Top 5 sperm-related GO terms in NMR ----
# "perm" matches "sperm", "spermatid", "spermatogenesis", etc.
gp4 <- ggplot(sperm_protein_nobg1[grep("perm",sperm_protein_nobg1$term_name),] %>% head(5),
              aes(x = reorder(term_name, -p_value),
                  y = -log10(p_value),color=source)) +
  geom_col() +
  coord_flip() +
  theme_bw() + theme(legend.position = "none") + 
  labs(x = "", y = "-log10(FDR)")+ scale_color_manual(values=c(
    "mouse"="purple",
    "NMR" = "gold"
  ))


# =============================================================================
# 3. Spot check: abundance of Irgc in NMR sperm
# =============================================================================

nmr_means <- rowMeans(sub_bre_proteomics_counts[3:ncol(sub_bre_proteomics_counts)])
names(nmr_means) <- sub_bre_proteomics_counts$Gene

# Abundance percentile of Irgc among all NMR proteins (100 = most abundant)
(rank(nmr_means)/length(nmr_means) * 100)["Irgc"]

# Mean abundance of Irgc
nmr_means["Irgc"]


# =============================================================================
# 4. Differential abundance of shared proteins
# =============================================================================

# Inner join on gene symbol: keeps only proteins detected in both species
# (treated as one-to-one orthologues)
merge_xspecies <- as.data.frame(merge(mouse_sperm,sub_bre_proteomics_counts,by="Gene"))

# Keep the first row per gene if a symbol appears more than once
merge_xspecies <- merge_xspecies[!duplicated(merge_xspecies$Gene),]
rownames(merge_xspecies) <- merge_xspecies$Gene

# Keep only the "Ratio" (normalized quantity) columns.
# Column order after merge: mouse columns first, then NMR columns.
merge_xspecies1 <- merge_xspecies[,grepl("Ratio",colnames(merge_xspecies))]

write.csv(merge_xspecies1,
          file = "NMR_mmus_one2one.csv",
          quote = FALSE, row.names = TRUE,col.names = TRUE)

# Per-protein Welch two-sample t-test.
# ASSUMES columns 1:3 are the mouse replicates and 4:end are NMR.
# Returns p-value, mouse mean ratio, NMR mean ratio.
DE <- apply(merge_xspecies1,1,function(x) {
  tt <- t.test(x[1:3],x[4:length(x)])
  p <- tt$p.value
  mouse <- mean(x[1:3])
  nmr <- mean(x[4:length(x)])
  return(c(p,mouse,nmr))
})
DE <- as.data.frame(t(DE))
colnames(DE) <- c("p","mouse","nmr")

# Benjamini-Hochberg FDR correction
DE$fdr <- p.adjust(DE$p,method = "fdr")

# Positive = higher in mouse, negative = higher in NMR
DE$log2FC <- log2(DE$mouse/DE$nmr)

library(ggplot2)

df <- DE
# Assume your df has columns: log2FC and pvalue (or padj)
df$negLog10P <- -log10(df$fdr)

# Define significance thresholds
fc_cutoff <- 1        # log2FC threshold (e.g., 2-fold change)
p_cutoff  <- 0.05      # p-value threshold

# Classify points
df$sig <- "Not significant"
df$sig[df$log2FC > fc_cutoff & df$fdr < p_cutoff]  <- "Mouse"
df$sig[df$log2FC < -fc_cutoff & df$fdr < p_cutoff] <- "NMR"

ggplot(df, aes(x = log2FC, y = negLog10P, color = sig)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_color_manual(values = c("Mouse" = "purple", 
                                "NMR" = "gold", 
                                "Not significant" = "grey70")) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +
  labs(x = "log2 Fold Change", 
       y = "-log10(p-value)", 
       color = "",
       title = "Volcano Plot") +
  theme_minimal() +
  theme(legend.position = "top")

df <- df[order(df$fdr),]
mouseDE <- gost(rownames(df)[df$sig == "Mouse"],evcodes = TRUE)
mouseDE1 <- mouseDE$result[mouseDE$result$term_size < 3e3 & grepl("GO",mouseDE$result$source),]
mouseDE1$source <- "mouse"
nmrDE <- gost(rownames(df)[df$sig == "NMR"],evcodes = TRUE)
nmrDE1 <- nmrDE$result[nmrDE$result$term_size < 3e3 & grepl("GO",nmrDE$result$source),]

nmrDE1$source <- "NMR"


gp5 <- ggplot(mouseDE1 %>% head(5),
              aes(x = reorder(term_name, -p_value),
                  y = -log10(p_value),color=source)) +
  geom_col() +
  coord_flip() +
  theme_bw() + theme(legend.position = "none") + 
  labs(x = "", y = "-log10(FDR)")+ scale_color_manual(values=c(
    "mouse"="purple",
    "NMR" = "gold"
  ))

gp6 <- ggplot(nmrDE1 %>% head(5),
              aes(x = reorder(term_name, -p_value),
                  y = -log10(p_value),color=source)) +
  geom_col() +
  coord_flip() +
  theme_bw() + theme(legend.position = "none") + 
  labs(x = "", y = "-log10(FDR)")+ scale_color_manual(values=c(
    "mouse"="purple",
    "NMR" = "gold"
  ))

# Last part, compare mouse vs. NMR RNAseq vs proteomics sperm 
rownames(mouse_nmr) <- mouse_nmr$Row.names
merge_xspecies <- as.data.frame(merge(df,mouse_nmr,by="row.names"))
rownames(merge_xspecies) <- merge_xspecies$Row.names

# Significance is from RNAseq
# Sig is proteomics

# NoDE
merge_xspecies$color <- "grey"

# Opposite DE
merge_xspecies$color[merge_xspecies$sig == "Mouse" & merge_xspecies$significance == "NMR"] <- "gold"
merge_xspecies$color[merge_xspecies$sig == "NMR" & merge_xspecies$significance == "Mouse"] <- "gold"

# same DE
merge_xspecies$color[merge_xspecies$sig == "Mouse" & merge_xspecies$significance == "Mouse"] <- "purple"
merge_xspecies$color[merge_xspecies$sig == "NMR" & merge_xspecies$significance == "NMR"] <- "purple"

# Proteomics Only
merge_xspecies$color[merge_xspecies$sig == "Mouse" & merge_xspecies$significance == "NS"] <- "pink"
merge_xspecies$color[merge_xspecies$sig == "NMR" & merge_xspecies$significance == "NS"] <- "pink"

# RNAseq only
merge_xspecies$color[merge_xspecies$sig == "Not significant" & merge_xspecies$significance == "Mouse"] <- "lightblue"
merge_xspecies$color[merge_xspecies$sig == "Not significant" & merge_xspecies$significance == "NMR"] <- "lightblue"



merge_xspecies$label <- ""
# merge_xspecies$label[merge_xspecies$col == "purple"] <- rownames(merge_xspecies)[merge_xspecies$col == "purple"]
df <- merge_xspecies
p1 <- ggplot(df, aes(x = log2FC, y = logFC, color = color,label = label)) +
  geom_point() +
  scale_color_identity() + theme_classic() +    geom_text_repel(
    data = subset(df, label != ""),
    aes(label = label),
    size = 6,
    max.overlaps = Inf
  ) + xlab("Proteomics Log2FC") + ylab("RNAseq Log2FC") + ggtitle(paste0("NMR vs. mouse proteomics vs. rnaseq: R^2=", round(cor(df$log2FC,df$logFC),3))) +   
  geom_smooth(method = "lm", se = FALSE, aes(group = 1), color = "black")
ggsave("Proteomics and RNAseq corrplot.png", p1,width = 4,height = 4)
##
### Check in with relaxed genes
##

mouseDE <- gost(rownames(df)[df$sig == "Mouse"],evcodes = TRUE)
mouseDE1 <- mouseDE$result[mouseDE$result$term_size < 3e3 & grepl("GO",mouseDE$result$source),]
mouseDE1$source <- "mouse"
nmrDE <- gost(rownames(df)[df$sig == "NMR"],evcodes = TRUE)
nmrDE1 <- nmrDE$result[nmrDE$result$term_size < 3e3 & grepl("GO",nmrDE$result$source),]

library(ActivePathways)
library(stringr)
relaxed_gmt <- readRDS("../sperm_rnaseq/relaxed_sperm_genes_gmt.rds")
apPlts <- list()

for(i in c("Mouse", "NMR")) {
  
  MouseDE <- merge_xspecies[merge_xspecies$sig == i,]
  head(as.data.frame(MouseDE$FDR,row.names = MouseDE$Row.names,col.names="FDR"))
  
  mat <- as.matrix(as.data.frame(MouseDE$FDR,row.names = MouseDE$Row.names,col.names="FDR"))
  colnames(mat) <- "FDR"
  rownames(merge_xspecies) <- merge_xspecies$Row.names
  bg <-rownames(merge_xspecies)
  bg2 <- union(makeBackground(relaxed_gmt), bg)
  bg3 <- intersect(makeBackground(relaxed_gmt), bg)
  
  #write.GMT(relaxed_gmt, "sperm_relax.gmt")
  
  AP <- ActivePathways(scores = mat,gmt = "../sperm_rnaseq/sperm_relax.gmt",background =bg3,significant = 1)
  
  AP <- ActivePathways(scores = mat,gmt = "../sperm_rnaseq/sperm_relax.gmt",background =bg2,significant = 1)
  
  
  AP$species <- unlist(lapply(strsplit(AP$term_id,"_"), function(x) x[1]))
  AP$state <- do.call("rbind",str_split(AP$term_id, "_(?=[A-Z])", n = 2))[,2]
  AP$state_ol <- paste0(lengths(AP$overlap),"/",AP$term_size)
  
  
  
  AP$species <- factor(AP$species,levels=c("Total","allmolerat","nmr"))
  
  AP <- AP[AP$species %in% c("Total","allmolerat","nmr"),]
  
  apPlts[[i]] <- AP_plt <- ggplot(AP, aes(x = term_id, y = -log10(adjusted_p_val))) +
    
    geom_col(aes(fill = species), color = NA) +
    
    geom_col(aes(color = state), fill = NA, linewidth = 1) +
    
    geom_text(aes(y=2,label = state_ol), vjust = 0.5, size = 6) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_manual(values = c(
      "RELAXED" = "blue",
      "INTENSIFIED" = "red",
      "NOT_SIGNIFICANT" = "grey"
    )) +
    
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12)) + coord_flip()
  
  pdf(file = paste0(i,"_DE_relaxed_olap_v0619.pdf"))
  print(AP_plt)
  dev.off()
}

ggsave("~/Dropbox-SickKidsResearch/Dustin Sokolowski/NMR_T2T/Figure_Data/Figure3/Figure3_KLM_mouse_vs_NMR/Mouse_DE_relaxed_olap_fdr_20260714.pdf",
       AP_plt)
ggsave("relaxed_mouse_DEGs.png",apPlts$Mouse)