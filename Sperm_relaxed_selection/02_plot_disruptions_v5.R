#!/usr/bin/env Rscript
# =============================================================================
# 02_plot_disruptions_v5.R  (v3 layout; honest absence labelling)
# Protein-level disruption plots from MACSE codon-aware NT sub-alignments.
# Layout matches v3: title "{GENE} protein alignment", subtitle with aa/species,
# bottom 7-marker legend in 2 rows; HetGla/BatSui no-sequence = "Absent".
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(optparse) })

option_list <- list(
  make_option("--input_dir", type="character",
              default="results/focal_alignments_17genes"),
  make_option("--output_dir", type="character",
              default="results/focal_alignments_17genes/figures"),
  make_option("--conservation_threshold", type="double", default=0.7),
  make_option("--format", type="character", default="both"),  # pdf | png | both
  make_option("--dpi", type="double", default=600),  # PNG resolution
  make_option("--width", type="double", default=13),
  make_option("--height_per_species", type="double", default=0.55),
  make_option("--base_height", type="double", default=2.6),
  make_option("--genes", type="character", default="")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)

FS_TITLE<-16; FS_SUBTITLE<-11; FS_YLAB<-12; FS_XLAB<-11; FS_AXIS_TIT<-13; FS_LEGEND<-11

COL_SEQ<-"#B3D9FF"; COL_GAP<-"#D0D0D0"
# Two DIFFERENT reasons a species has no sequence. They are not the same claim.
#   unannotated      gene not in that species' proteome/orthogroup. tBLASTn
#                    recovers the locus from the genome for several of these,
#                    so the gene very likely exists but is unannotated.
#   not_single_copy  NMR/BatSui come from curated single-copy TSVs; a gene
#                    excluded from that call yields no sequence. A filtering
#                    decision, NOT evidence of gene loss.
COL_UNANNOT<-"#FFE8CC"   # pale orange
COL_NOTSC  <-"#E4E4F0"   # pale lilac
COL_SUBST<-"#E67E22"; COL_STOP<-"#C0392B"; COL_FRAME<-"#2980B9"
GROUP_COLORS <- c(Bathyergidae="#B2182B", Hystricomorpha="#E08214",
                  Muridae="#1B7837", Primates="#2166AC",
                  Laurasiatheria="#5E3C99", Afrotheria="#01665E")

codon_table <- c(
  TTT="F",TTC="F",TTA="L",TTG="L",CTT="L",CTC="L",CTA="L",CTG="L",
  ATT="I",ATC="I",ATA="I",ATG="M",GTT="V",GTC="V",GTA="V",GTG="V",
  TCT="S",TCC="S",TCA="S",TCG="S",CCT="P",CCC="P",CCA="P",CCG="P",
  ACT="T",ACC="T",ACA="T",ACG="T",GCT="A",GCC="A",GCA="A",GCG="A",
  TAT="Y",TAC="Y",TAA="*",TAG="*",CAT="H",CAC="H",CAA="Q",CAG="Q",
  AAT="N",AAC="N",AAA="K",AAG="K",GAT="D",GAC="D",GAA="E",GAG="E",
  TGT="C",TGC="C",TGA="*",TGG="W",CGT="R",CGC="R",CGA="R",CGG="R",
  AGT="S",AGC="S",AGA="R",AGG="R",GGT="G",GGC="G",GGA="G",GGG="G")
translate_codon <- function(cod){
  cod<-toupper(cod)
  if(grepl("!",cod,fixed=TRUE)) return("!")
  if(grepl("-",cod,fixed=TRUE)) return("-")
  if(nchar(cod)<3) return("X")
  aa<-codon_table[[cod]]; if(is.null(aa)) return("X"); aa }
read_fasta <- function(path){
  lines<-readLines(path); seqs<-list(); cur<-NULL; buf<-c()
  for(ln in lines){
    if(startsWith(ln,">")){
      if(!is.null(cur)) seqs[[cur]]<-paste0(buf,collapse="")
      cur<-sub("^>","",strsplit(ln,"\\s+")[[1]][1]); buf<-c()
    } else buf<-c(buf,trimws(ln)) }
  if(!is.null(cur)) seqs[[cur]]<-paste0(buf,collapse=""); seqs }

build_gene_df <- function(gene, input_dir, cons_thr){
  aln_path<-file.path(input_dir,paste0(gene,"_focal_NT.aln"))
  pres_path<-file.path(input_dir,paste0(gene,"_presence.tsv"))
  info_path<-file.path(input_dir,"species_info.tsv")
  if(!file.exists(aln_path)||!file.exists(info_path)) return(NULL)
  info<-read.delim(info_path,stringsAsFactors=FALSE)
  # controls extractor writes one gene-keyed species_info; subset to this gene
  if("gene" %in% names(info)) info<-info[info$gene==gene,,drop=FALSE]
  presence<-if(file.exists(pres_path)) read.delim(pres_path,stringsAsFactors=FALSE) else data.frame(species=info$species,present=1,reason="")
  if(!"reason" %in% names(presence)) presence$reason<-""
  seqs<-read_fasta(aln_path); if(length(seqs)==0) return(NULL)
  n_codons<-floor(nchar(seqs[[1]])/3); if(n_codons<1) return(NULL)
  present_species<-intersect(info$species,names(seqs))
  aa_mat<-matrix("-",nrow=length(present_species),ncol=n_codons,dimnames=list(present_species,NULL))
  for(sp in present_species){
    s<-seqs[[sp]]
    for(j in seq_len(n_codons)) aa_mat[sp,j]<-translate_codon(substr(s,(j-1)*3+1,(j-1)*3+3)) }
  consensus<-character(n_codons); conserved<-logical(n_codons)
  for(j in seq_len(n_codons)){
    col<-aa_mat[,j]; valid<-col[!(col %in% c("-","*","!","X"))]
    if(length(valid)==0){consensus[j]<-NA; conserved[j]<-FALSE; next}
    tb<-sort(table(valid),decreasing=TRUE); consensus[j]<-names(tb)[1]
    conserved[j]<-(tb[1]/length(present_species))>=cons_thr }
  order_map<-setNames(info$order,info$species); rows<-list(); ri<-1
  for(sp in info$species){
    is_present<-presence$present[match(sp,presence$species)]; is_present<-ifelse(is.na(is_present),0,is_present)
    yorder<-order_map[[sp]]
    if(is_present==0){
      why<-presence$reason[match(sp,presence$species)]; if(is.na(why)||why=="") why<-"unannotated"
      tt<-if(why=="not_single_copy") "absent_notsc" else "absent_unannot"
      rows[[ri]]<-data.frame(species=sp,yorder=yorder,pos=seq_len(n_codons),type=tt,stringsAsFactors=FALSE); ri<-ri+1; next }
    for(j in seq_len(n_codons)){
      aa<-aa_mat[sp,j]
      type<-if(aa=="-")"gap" else if(aa=="*")"stop" else if(aa=="!")"frameshift" else if(conserved[j]&&!is.na(consensus[j])&&aa!=consensus[j])"subst" else "present"
      rows[[ri]]<-data.frame(species=sp,yorder=yorder,pos=j,type=type,stringsAsFactors=FALSE); ri<-ri+1 } }
  df<-do.call(rbind,rows)
  attr(df,"n_codons")<-n_codons; attr(df,"info")<-info; attr(df,"presence")<-presence; df }

plot_gene <- function(gene, df, out_dir, fmt, width, hps, base_h){
  info<-attr(df,"info"); n_codons<-attr(df,"n_codons"); presence<-attr(df,"presence")
  info<-info[order(info$order),]; n_species<-nrow(info)
  absent_sp<-presence$species[presence$present==0]
  absent_labels<-info$label[match(absent_sp,info$species)]
  n_present<-n_species-length(absent_sp)
  df$y<-max(info$order)-df$yorder
  ylabs<-info$label; ybreaks<-max(info$order)-info$order; ycols<-GROUP_COLORS[info$group]
  base_df<-df
  base_df$fill<-ifelse(base_df$type=="gap","gap",
                ifelse(base_df$type=="absent_unannot","absent_unannot",
                ifelse(base_df$type=="absent_notsc","absent_notsc","present")))
  subst<-df[df$type=="subst",]; stops<-df[df$type=="stop",]; frames<-df[df$type=="frameshift",]
  lv<-c("Sequence","Gap","Substitution at conserved site","Stop codon","Frameshift")
  legend_df<-data.frame(pos=rep(-Inf,5),y=rep(-Inf,5),label=factor(lv,levels=lv))

  p<-ggplot()+
    geom_tile(data=base_df,aes(x=pos,y=y,fill=fill),height=0.82)+
    scale_fill_manual(values=c(present=COL_SEQ,gap=COL_GAP,absent_unannot=COL_UNANNOT,absent_notsc=COL_NOTSC),guide="none")
  if(nrow(stops)) p<-p+geom_point(data=stops,aes(x=pos,y=y),shape=17,color=COL_STOP,size=2.0)
  if(nrow(frames)) p<-p+geom_point(data=frames,aes(x=pos,y=y),shape=18,color=COL_FRAME,size=2.3)
  if(nrow(subst)) p<-p+geom_point(data=subst,aes(x=pos,y=y),shape=16,color=COL_SUBST,size=1.5)
  if(length(absent_sp)>0){
    # NEVER "GENE ABSENT": tBLASTn finds the locus in the genome for several of
    # these species (SPEM3/CABS1 in gpig, chin, degu at 51-59% identity), so
    # asserting loss would be false. State only what is known.
    why_v<-presence$reason[match(absent_sp,presence$species)]
    why_v[is.na(why_v)|why_v==""]<-"unannotated"
    lab_v<-ifelse(why_v=="not_single_copy","absent","unannotated / not detected")
    col_v<-ifelse(why_v=="not_single_copy","#6C5B9E","#B8651B")
    aden<-data.frame(pos=n_codons/2,
                     y=max(info$order)-info$order[match(absent_sp,info$species)],
                     lab=lab_v,col=col_v,stringsAsFactors=FALSE)
    p<-p+geom_text(data=aden,aes(x=pos,y=y,label=lab),colour=aden$col,
                   fontface="bold.italic",size=3.6) }
  p<-p+
    geom_point(data=legend_df,aes(x=pos,y=y,shape=label,color=label),inherit.aes=FALSE)+
    scale_shape_manual(name=NULL,values=c("Sequence"=15,"Gap"=15,"Substitution at conserved site"=16,"Stop codon"=17,"Frameshift"=18))+
    scale_color_manual(name=NULL,values=c("Sequence"=COL_SEQ,"Gap"=COL_GAP,"Substitution at conserved site"=COL_SUBST,"Stop codon"=COL_STOP,"Frameshift"=COL_FRAME))+
    guides(shape=guide_legend(nrow=1,override.aes=list(size=3.5)),color=guide_legend(nrow=1))
  subtitle<-paste0(n_codons," aa | ",n_present," species with sequence",
                   if(length(absent_sp)>0) paste0(" | no sequence: ",paste(absent_labels,collapse=", ")) else "")
  p<-p+
    scale_y_continuous(breaks=ybreaks,labels=ylabs,expand=expansion(add=0.6))+
    scale_x_continuous(expand=expansion(mult=c(0.01,0.01)))+
    labs(x="Alignment position (aa)",y=NULL,title=paste0(gene," protein alignment"),subtitle=subtitle)+
    theme_minimal(base_size=12)+
    theme(plot.title=element_text(face="bold",size=FS_TITLE),
          plot.subtitle=element_text(size=FS_SUBTITLE,color="grey35"),
          axis.text.y=element_text(color=ycols,face="bold",size=FS_YLAB),
          axis.text.x=element_text(size=FS_XLAB),
          axis.title.x=element_text(size=FS_AXIS_TIT),
          panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
          panel.grid.major.x=element_line(color="grey90",linewidth=0.3),
          legend.position="bottom",legend.text=element_text(size=FS_LEGEND),
          legend.margin=margin(2,2,2,2),plot.margin=margin(12,16,10,12))
  h<-base_h+n_species*hps
  fmts <- if(fmt=="both") c("pdf","png") else fmt
  for(ff in fmts){
    out<-file.path(out_dir,paste0(gene,"_disruptions.",ff))
    if(ff=="png") ggsave(out,p,width=width,height=h,dpi=opt$dpi,limitsize=FALSE)
    else ggsave(out,p,width=width,height=h,limitsize=FALSE)
    cat(sprintf("  saved %s  (%d aa, %d species, %d absent)\n",out,n_codons,n_species,length(absent_sp)))
  } }

main <- function(){
  if(nzchar(opt$genes)) genes<-trimws(strsplit(opt$genes,",")[[1]])
  else { fs<-list.files(opt$input_dir,pattern="_focal_NT\\.aln$"); genes<-sub("_focal_NT\\.aln$","",fs) }
  cat("======================================================\n")
  cat("02_plot_disruptions_v5 (v3 layout):",length(genes),"genes\n")
  cat("  conservation threshold:",opt$conservation_threshold,"\n")
  cat("======================================================\n")
  for(g in genes){
    df<-tryCatch(build_gene_df(g,opt$input_dir,opt$conservation_threshold),error=function(e){cat("  ",g,"ERROR:",conditionMessage(e),"\n");NULL})
    if(is.null(df)){cat("  ",g,"-> skipped (no data)\n");next}
    plot_gene(g,df,opt$output_dir,opt$format,opt$width,opt$height_per_species,opt$base_height) }
  cat("======================================================\nDONE\n") }
main()
