setwd("~/Desktop/NMR_multi_assembly/v090226/")

pri_repeat <- read.table("v090226.pri.repeat.bed",
                         header=FALSE, as.is=TRUE,
                         sep="\t")

alt_repeat <- read.table("v090226.alt.repeat.bed",
                         header=FALSE, as.is=TRUE,
                         sep="\t")

# Active ORFs

repeats <- list(pri=pri_repeat,
                alt=alt_repeat)

for(i in c("pri","alt")) {
  
  rep <- repeats[[i]]
  
  rep <- rep[(rep$V3-rep$V2) >= 2e3,]
  rep <- rep[grep("LINE/L1",rep$V4),]
  
 # rep_pos <- rep[rep$V6 == "+",]
#  rep_neg <- rep[rep$V6 == "-",]
  
  
  write.table(rep_pos,file = paste0("annotations/",i,".2kbp.L1.bed"),
              quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  write.table(rep_neg,file = paste0("annotations/",i,".2kbp.L1.neg.bed"),
              quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  
}

##
### After interproscan
##

setwd("/Users/dsokolowski/Desktop/NMR_multi_assembly/Repeats/ActiveLines/")


library(stringr)

for(p in c("pri","alt")) {
  
  Interpro <- read.delim(paste0("mHetGla1.1.",p,".L1.2Kb.ORF.fa.tsv"),
                         header=FALSE,as.is = TRUE)
  
  Interpro <- Interpro[Interpro$V4 == "Pfam",]
  x <- sub(".*_(chr[0-9XY]+:[0-9]+-[0-9]+).*", "\\1", Interpro$V1)
  
  
  df <- str_match(x, "(chr[0-9XYM]+):(\\d+)-(\\d+)")
  
  df <- data.frame(
    region = x,
    chr = df[,2],
    start = as.integer(df[,3]),
    end = as.integer(df[,4]),
    stringsAsFactors = FALSE
  )
  
  # Define genome order (like UCSC/Ensembl style)
  chr_levels <- c(paste0("chr", 1:29), "chrX", "chrY", "chrM")
  
  df$chr <- factor(df$chr, levels = chr_levels, ordered = TRUE)
  
  Interpro$chr <- df$chr
  Interpro$start <- df$start
  Interpro$end <- df$end
  Interpro$base_coord <- df$region
  Interpro_sort <- Interpro[order(df$chr, df$start, df$end),]
  
  keys <- paste0("k", match(Interpro_sort$base_coord, unique(Interpro_sort$base_coord)))
  
  
  
  
  Interpro_sort$key <- keys
  
  k2 <- names(table(Interpro_sort$key))[table(Interpro_sort$key) > 1]
  
  Interpro_sort_k2 <- Interpro_sort[Interpro_sort$key %in% k2,]
  
  
  Interpro_Active <- list()
  Interpro_Assign <- c()
  
  all_domains <- c("Reverse transcriptase domain","L1 transposable element, RRM domain","Endonuclease/exonuclease/phosphatase")

  
  for(i in unique(Interpro_sort_k2$base_coord)) {
    inter_each <- Interpro_sort_k2[Interpro_sort_k2$base_coord == i,]
    
    domain_inter <- intersect(inter_each$V13,all_domains)
    if(length(domain_inter) == 3) {
      Interpro_Active[[i]] <- inter_each[1,c("chr","start","end")]
      Interpro_Assign <- c(Interpro_Assign,"Active")
    } 
    if(all(domain_inter %in% all_domains[1:2])) {
      Interpro_Active[[i]] <- inter_each[1,c("chr","start","end")]
      Interpro_Assign <- c(Interpro_Assign,"Active-noEN")
    }
    
  }
  
  Active_bed <- do.call("rbind",Interpro_Active)
  Active_bed$Assign <- Interpro_Assign
  Active_bed <- Active_bed[complete.cases(Active_bed),]
  write.table(Active_bed, file = paste0("mHetGla1.1.",p,"_activeL1s.bed"),
              quote=FALSE,row.names = FALSE,col.names=FALSE,sep="\t")

}

Active_bed <- read.table(paste0("mHetGla1.1.","pri","_activeL1s.bed"),header=F,as.is=T,sep="\t")

pri2Kbp <- read.table("pri.2kbp.L1.bed",header=F,as.is=T,sep="\t")

library(GenomicFeatures)
txdb <- makeTxDbFromGFF("~/Desktop/NMR_multi_assembly/v090226/annotations/MT_annotations/mHetGla1.pri.MT.v020326.gff")
library(rtracklayer)
gff <- readGFF("~/Desktop/NMR_multi_assembly/v090226/annotations/MT_annotations/mHetGla1.pri.MT.v020326.gff")
gff <- gff[gff$strand != "*",]
export(gff,con = "~/Desktop/NMR_multi_assembly/v090226/annotations/MT_annotations/mHetGla1.pri.MT.v020326.strandfilter.gff",format = "gff3")
txdb <- makeTxDbFromGFF("~/Desktop/NMR_multi_assembly/v090226/annotations/MT_annotations/mHetGla1.pri.MT.v020326.strandfilter.gff")
colnames(Active_bed) <- c("chr","start","end","name")
df <- makeGRangesFromDataFrame(Active_bed)
library(Chipseeker)

ensembltest <- annotatePeak(df, TxDb = txdb)
plotAnnoPie(ensembltest)
anno <- ensembltest@anno

promoter <- anno[abs(anno$distanceToTSS) < 1e3,]

# Add gene symbol

mHetGla1.1.pri <- read.table("~/Desktop/NMR_multi_assembly/v090226/annotations/mHetGla1.1.pri_gene_symbols_metrics.tsv",header=T,as.is=T,sep="\t")
rownames(mHetGla1.1.pri) <- mHetGla1.1.pri$mikado_id
mHetGla1.1.pri1 <- mHetGla1.1.pri[promoter$geneId,]
promoter$hierarchical_gene <- mHetGla1.1.pri1$hierarchical_gene

gp2 <- gprofiler2::gost(unique(promoter$hierarchical_gene),organism = "mmusculus",correction_method = "fdr",evcodes = T)
uni <- unlist(lapply(strsplit(promoter$uniqueGene,"-copy"), function(x) return(x[1])))

write.table(unique(promoter$hierarchical_gene), file = "mHetGla1.1.ActiveLine.txt",
            quote = FALSE, row.names=FALSE, col.names = FALSE, sep = "\t")
library(Biostrings)
setwd("~/Desktop/NMR_multi_assembly/v090226/")

fa <- readDNAStringSet("mHetGla1.1.pri.cleaned.fa")
chrX <- fa["chrX"]
chrY <- fa["chrY"]

writeXStringSet(chrX, filepath = "mHetGla1.1.chrX.fa")
writeXStringSet(chrY, filepath = "mHetGla1.1.chrY.fa")

## Install required packages
install.packages("devtools")
library(devtools)

## Install from GitHub repository
devtools::install_github("daewoooo/SVbyEye", branch="master")
library(SVbyEye)

paf <- readPaf("~/Downloads/map_mHetGla1.1.chrY_to_mHetGla1.1.chrX.paf")

setwd("~/Desktop/NMR_multi_assembly/Cenpa/")
fimo_cenbp <- read.table("Cenbp_cenpa_pos.fimo.tsv",header=FALSE,as.is=TRUE,sep="\t")

cenpb <- table(fimo_cenbp$V3)
cenpb <- cenpb[paste0("chr",c(1:29,"X","Y"))]
write.csv(cenpb, file = "fimo_cenbp_motif.csv",
          quote = FALSE, row.names = FALSE, col.names = FALSE)

annotations <- read.delim("mHetGla Centromere Locations - Centromere Heatmap Annotations.tsv")
for(i in 7:14) {
  annotations[,i] <- toupper(annotations[,i])
}

annotation_heat <- annotations[,7:14]
rownames(annotation_heat) <- annotations$Chromosome

for(i in 1:ncol(annotation_heat)) {
  annotation_heat[,i] <- as.numeric(factor(annotation_heat[,i]))-1
}
apply(annotation_heat,2,factor)
annotation_col <- data.frame(
  cenpb=annotations$Cenpb.motifs,
  KBp= round(annotations$Length..bp./1e3)
)
rownames(annotation_col) <- rownames(annotation_heat)
color <- c("white","black")
anno_transpose <- t(annotation_heat)
rownames(anno_transpose) <- c("SFC","Gene","CDR","Individual Variation", "Dual Peak",)
p1 <- pheatmap::pheatmap(t(annotation_heat),color = color, cluster_rows = FALSE,cluster_cols = FALSE,cellwidth = 10, cellheight = 10,
                   annotation_col = annotation_col,legend = FALSE)

pdf(file = "mHetGla.centromere.annotation.pdf")
print(p1)
dev.off()

table(annotation_heat$SFC.logical, annotation_heat$Active.L1)

fisher.test(table( annotation_heat$Gene.logical,annotation_heat$SFC.logical)[c("1","0"),c("1","0")])

tab <- table(annotation_heat$SFC.logical,annotation_heat$Looping)


rownames(tab) <- paste0("SFC",rownames(tab))
colnames(tab) <- paste0("looping",colnames(tab))

# SFCs are generally more Cenpb motif rich (p_wilcoxon = 0.02398)
# SFCs are on average shorter (p_wilcoxon = 0.000184)
# wilcox.test(annotation_col$cenpb ~ annotation_heat$SFC.logical)
wilcox.test(annotation_col$KBp ~ annotation_heat$SFC.logical)
t.test(annotation_col$KBp ~ annotation_heat$SFC.logical)

# Codon counts

musT2T <- read.csv("~/Downloads/musT2T.counts.plotly.csv",
                   header=TRUE,as.is=TRUE)

musT2T_cut <- musT2T[musT2T$repeat_length >= 19,]

hsapT2T <- read.csv("~/Downloads/chm13v2.0_maskedY_rCRS.fa.PanSN.counts.plotly.csv",
                    header=TRUE,as.is=TRUE)

hsapT2T_cut <- hsapT2T[hsapT2T$repeat_length >= 19,]

nmrT2T <- read.csv("~/Downloads/mHetGla1.1.pri.cleaned.counts.plotly.csv",
                    header=TRUE,as.is=TRUE)
nmrT2T_cut <- nmrT2T[nmrT2T$repeat_length >=19,]


nmrT2T_cut$BP <- nmrT2T_cut$repeat_length*nmrT2T_cut$count
musT2T_cut$BP <- musT2T_cut$repeat_length*musT2T_cut$count
hsapT2T_cut$BP <- hsapT2T_cut$repeat_length*hsapT2T_cut$count

hsap_20copy_bp <- aggregate(hsapT2T_cut$BP~hsapT2T_cut$codon,FUN = sum)
mmus_20copy_bp <- aggregate(musT2T_cut$BP~musT2T_cut$codon,FUN = sum)
nmr_20copy_bp <- aggregate(nmrT2T_cut$BP~nmrT2T_cut$codon,FUN = sum)

colnames(hsap_20copy_bp) <- colnames(mmus_20copy_bp) <- colnames(nmr_20copy_bp) <- c("codon","BP")
merged.hs.mm <- merge(hsap_20copy_bp,mmus_20copy_bp,"codon")

colnames(merged.hs.mm) <- c("codon","hsap","mmus")
merged.all <- merge(merged.hs.mm,nmr_20copy_bp,"codon")
colnames(merged.all) <-  c("codon","hsap","mmus","nmr")
rownames(merged.all) <- merged.all$codon
merged.mat <- merged.all[,colnames(merged.all) != "codon"]

col_fun <- colorRampPalette(c("white", "red"))

pdf(file = "TNR_repeat_number.pdf",height = 11,width = 8)
print(pheatmap::pheatmap(log10(merged.mat),display_numbers = TRUE,
                   color = col_fun(10),fontsize = 10,
                   cellheight = 25, cellwidth = 25,
                   cluster_rows = TRUE,cluster_cols = TRUE))
dev.off()

pdf(file = "TNR_repeat_number_width.pdf",height = 11,width = 8)
print(pheatmap::pheatmap(t(log10(merged.mat)),display_numbers = TRUE,
                         color = col_fun(10),fontsize = 8,
                         cellheight = 20, cellwidth = 20,
                         cluster_rows = TRUE,cluster_cols = TRUE))
dev.off()


pheatmap::pheatmap(t(merged.mat/1e3),display_numbers = TRUE,
                   color = col_fun(10),fontsize = 12,
                   cellheight = 30, cellwidth = 30,
                   cluster_rows = TRUE,cluster_cols = TRUE)

pheatmap::pheatmap(log2(merged.mat),display_numbers = TRUE,
                   color = col_fun(10))
