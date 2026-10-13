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
  
  rep_pos <- rep[rep$V6 == "+",]
  rep_neg <- rep[rep$V6 == "-",]
  
  
  write.table(rep_pos,file = paste0("annotations/",i,".2kbp.L1.pos.bed"),
              quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  write.table(rep_neg,file = paste0("annotations/",i,".2kbp.L1.neg.bed"),
              quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  
}
