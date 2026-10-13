#!/usr/bin/env Rscript
# =============================================================================
# part5b_sperm_gene_list.R
#
# Extracts human gene symbols from sperm-related GO terms (including all
# descendant terms via GO2ALLEGS), deduplicates, and writes a clean gene
# list for the downstream 600-gene pipeline.
#
# GO terms:
#   GO:0007283  spermatogenesis
#   GO:0030317  flagellated sperm motility
#   GO:0048232  male gamete generation
#
# OUTPUT (results/10_extended_sperm_genes/):
#   go_sperm_genes.tsv         SYMBOL, ENTREZID, GENENAME, GO_TERMS, GO_IDS
#   go_sperm_gene_symbols.txt  one SYMBOL per line (sorted, unique)
#
# Symbols are HUMAN, UPPERCASE (e.g. ZAN, MROH2B). Downstream focal-ortholog
# matching against the mole-rat TSVs (which use Title-case, e.g. Zan) MUST be
# case-insensitive
#
# Usage:  Rscript part5b_sperm_gene_list.R
# Requires: BiocManager::install(c("org.Hs.eg.db", "GO.db"))
# =============================================================================

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(GO.db)
})

OUTPUT_DIR <- "results/10_extended_sperm_genes"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

GO_TERMS <- list(
  "GO:0007283" = "spermatogenesis",
  "GO:0030317" = "flagellated_sperm_motility",
  "GO:0048232" = "male_gamete_generation"
)

cat("Part 10a: extract sperm-related gene lists from GO\n\n")

all_results <- data.frame()

for (go_id in names(GO_TERMS)) {
  go_name <- GO_TERMS[[go_id]]
  cat(sprintf("Processing %s (%s)...\n", go_id, go_name))

  entrez_ids <- tryCatch({
    unique(get(go_id, org.Hs.egGO2ALLEGS))
  }, error = function(e) {
    res <- select(org.Hs.eg.db, keys = go_id,
                  columns = "ENTREZID", keytype = "GOALL")
    unique(res$ENTREZID)
  })
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  cat(sprintf("  %d Entrez IDs (incl. descendant terms)\n", length(entrez_ids)))

  mapping <- select(org.Hs.eg.db,
                    keys = entrez_ids,
                    columns = c("SYMBOL", "GENENAME"),
                    keytype = "ENTREZID")
  mapping <- unique(mapping[!is.na(mapping$SYMBOL),
                            c("ENTREZID", "SYMBOL", "GENENAME")])
  mapping$GO_ID   <- go_id
  mapping$GO_TERM <- go_name

  cat(sprintf("  %d unique gene symbols\n", length(unique(mapping$SYMBOL))))
  all_results <- rbind(all_results, mapping)
}

cat("\nCombining + deduplicating...\n")

go_membership <- aggregate(GO_TERM ~ SYMBOL, data = all_results,
                           FUN = function(x) paste(sort(unique(x)), collapse = ";"))
names(go_membership)[2] <- "GO_TERMS"

go_ids <- aggregate(GO_ID ~ SYMBOL, data = all_results,
                    FUN = function(x) paste(sort(unique(x)), collapse = ";"))
names(go_ids)[2] <- "GO_IDS"

dedup <- unique(all_results[, c("SYMBOL", "ENTREZID", "GENENAME")])
dedup <- dedup[!duplicated(dedup$SYMBOL), ]
dedup <- merge(dedup, go_membership, by = "SYMBOL", all.x = TRUE)
dedup <- merge(dedup, go_ids,       by = "SYMBOL", all.x = TRUE)
# Uppercase symbols for consistency with human convention (defensive)
dedup$SYMBOL <- toupper(dedup$SYMBOL)

# The 17 targeted genes are members of this set by design; several lack the
# three GO terms (e.g. NDC80: cell division; ZAN: fertilization), so they are
# added explicitly rather than assumed.
targets17 <- c("ADIG","ATP1A4","CABS1","DNHD1","DPY19L2","GARIN3","GARIN1A",
               "FSIP2","MROH2B","NDC80","NPM2","RIMBP3","SPATA33","SPEM3",
               "TMPRSS12","TTLL8","ZAN")
missing17 <- setdiff(targets17, dedup$SYMBOL)
if (length(missing17) > 0) {
  dedup <- rbind(dedup, data.frame(
    SYMBOL = missing17, ENTREZID = NA,
    GENENAME = "targeted gene (17-gene study)",
    GO_TERMS = "targeted_17", GO_IDS = NA))
  cat(sprintf("\nAdded %d targeted genes not covered by the GO terms: %s\n",
              length(missing17), paste(missing17, collapse = ", ")))
}
dedup <- dedup[order(dedup$SYMBOL), ]

cat(sprintf("\nTotal unique genes: %d\n", nrow(dedup)))

cat("\nPer-GO-term counts:\n")
for (go_id in names(GO_TERMS)) {
  n <- length(unique(all_results$SYMBOL[all_results$GO_ID == go_id]))
  cat(sprintf("  %s %-28s : %d\n", go_id, GO_TERMS[[go_id]], n))
}

in_sperm  <- unique(all_results$SYMBOL[all_results$GO_ID == "GO:0007283"])
in_motil  <- unique(all_results$SYMBOL[all_results$GO_ID == "GO:0030317"])
in_gamete <- unique(all_results$SYMBOL[all_results$GO_ID == "GO:0048232"])
cat(sprintf("\n  union (total unique): %d\n",
            length(unique(c(in_sperm, in_motil, in_gamete)))))

out_genes <- file.path(OUTPUT_DIR, "go_sperm_genes.tsv")
write.table(dedup, out_genes, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\nSaved: %s (%d genes)\n", out_genes, nrow(dedup)))

out_symbols <- file.path(OUTPUT_DIR, "go_sperm_gene_symbols.txt")
writeLines(sort(unique(dedup$SYMBOL)), out_symbols)
cat(sprintf("Saved: %s\n", out_symbols))

cat("\nNext: python3 part5d_select_orthogroups.py --gene_list results/10_extended_sperm_genes/go_sperm_genes.tsv --outdir results/05_sets/sperm --min_species 15 ...\n")
