
swap_symbols <- function(gene_symbols , symbol_colname, gene_key, matching_column = "transcriptID" ) {
  
  # This function swaps the transcript ID with the gene symbol for a specific row of the "gene_symbols" table.
  
  # gene_symbols is the "gene_symbols.tsv" file outputted in GenAnT and loaded into R with "gene_symbols <- read.table(gene_symbols.tsv",header=T,as.is=T,sep="\t")
  
  # symbol_colname is the name of the column whose gene symbols you're transferring. e.g., "liftoff_gene"
  
  # gene_key is the $prefix.clean.genekkey.tsv file read into R. An example of for mouse is below "genekey_mouse <-  read.table("GCF_000001635.27_GRCm39_genomic.genekkey.tsv",header=T,as.is=T,sep="\t")
  
  # matching_column is the colname of the genekey file that matches "symbol_colname"
  
  transcriptID <- gene_symbols[,symbol_colname]
  
  if(matching_column == "geneID") {
    
    gene_key <- gene_key[!duplicated(gene_key$geneID),]
    rownames(gene_key) <- gene_key[,matching_column]
    gene_key$geneSymbol[is.na(gene_key$geneSymbol)] <- gene_key$geneID[is.na(gene_key$geneSymbol)]
    gene_key$geneSymbol[gene_key$geneSymbol == ""] <- gene_key$geneID[gene_key$geneSymbol == ""]
    
    gene_key1 <- gene_key[transcriptID,]
    
  } else {
    
    rownames(gene_key) <- gene_key[,matching_column]
    gene_key$geneSymbol[is.na(gene_key$geneSymbol)] <- gene_key$geneID[is.na(gene_key$geneSymbol)]
    gene_key$geneSymbol[gene_key$geneSymbol == ""] <- gene_key$geneID[gene_key$geneSymbol == ""]
    
    gene_key1 <- gene_key[transcriptID,]
    
  }
  
  
  
  return(gene_key1$geneSymbol)
  
}

refDir="~/Desktop/Heliophobius/"

genekey_nmr <- read.table(paste0(refDir,"mHetGla1.1.pri.cleaned.genekkey.tsv"),header=TRUE,as.is=TRUE,sep="\t")
genekey_mouse <- read.table(paste0(refDir,"GCF_000001635.27_GRCm39_genomic.genekkey.tsv"),header=T,as.is=T,sep="\t")

setwd("~/Desktop/NMR_multi_assembly/v090226/annotations/alt_annos/")
annoDir="~/Desktop/NMR_multi_assembly/v090226/annotations/alt_annos/"
gene_symbols <- read.table(paste0(annoDir,"mHetGla1.alt.gene_symbols_noCopies.tsv"),header=TRUE,as.is=TRUE,sep="\t")

# fix liftoff_gene
gene_symbols$liftoff_gene <- swap_symbols(gene_symbols = gene_symbols,symbol_colname = "liftoff_gene",
                                          gene_key = genekey_nmr,matching_column = "transcriptID")

# fix togar1_gene
gene_symbols$togar1_gene <- swap_symbols(gene_symbols = gene_symbols,symbol_colname = "togar1_gene",
                                         gene_key = genekey_mouse,matching_column = "transcriptID")

gene_symbols <- gene_symbols[,(colnames(gene_symbols) != "toga_gene")]



load("~/Downloads/bioMart_ortholog_human_mouse.rda")

# As an example, here is a very simple script to convert "liftoff_gene" and "togar2_gene" symbols to mouse.



for(i in c("liftoff_gene")) {
  genes <- gene_symbols[,i]
  gene_symbols[,i] <- ifelse(
    !is.na(match(genes, bioMart_orthologs$human)), 
    bioMart_orthologs$mouse[match(genes, bioMart_orthologs$human)], 
    genes
  )
}


# NOTE. This stage will be quite different outside of mammals but the general principal will follow. Find the most closesly related model organism and have your gene symbols match that. It'll also be helpful for gene-set enrichment later on.

for(j in 2:ncol(gene_symbols)) {
  gene_symbols[,j] <- sub("-copy[0-9]+$", "", gene_symbols[,j])
}


gene_symbols$liftoff_gene <- gsub("mikado.","Hg1pri.",gene_symbols$liftoff_gene)

write.table(gene_symbols, file = "postprocess_gene_symbols_noCopies.tsv", quote=FALSE, 
            row.names=FALSE, col.names=TRUE, sep = "\t")



