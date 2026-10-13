#!/bin/bash
#$ -l h_vmem=16G,h_rt=150:00:00,h_stack=32M
#$ -pe smp 16

###########################
##
### This script uses RagTag and bedtools to identify the amount of sequences that is missing from mHetGla1.1.pri
##
###########################

source /.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate RagTag

cd /.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/PreviousAssemblyMissingBP

refFA="/.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/v090226/mHetGla1.1.pri.cleaned.fa" # directory in /data -- adding species can be done with scripts in /utils
refbed="/.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/v090226/mHetGla1.1.pri.cleaned.chrom.sizes.bed" # directory in /data -- adding species can be done with scripts in /utils

fa=$I
prefix=`basename ${fa} .fa`

mkdir -p ${prefix}.align 

# ragtag.py correct -t 16 -o ${prefix}.ragtag  $refFA $fa
# ragtag.py scaffold -t 16 -r -m 10000000 -o ${prefix}.ragtag $refFA ${prefix}.ragtag/ragtag.correct.fasta 
# HetGla.tor.d5d682.dip.hap2.scaffold.fasta

cd ${prefix}.align


minimap2 -t 4 -x asm5  $refFA ../${fa} > ${prefix}.paf

awk -F'\t' '$1 == $6 && $12 > 0' ${prefix}.paf > ${prefix}.filtered.paf
cut -f 6,8,9 ${prefix}.filtered.paf > ${prefix}.filtered.bed
bedtools subtract -a ${refbed} -b ${prefix}.filtered.bed > mHetGla1.no.${prefix}.bed
bedtools sort -i mHetGla1.no.${prefix}.bed > ../mHetGla1.no.${prefix}.s.bed


