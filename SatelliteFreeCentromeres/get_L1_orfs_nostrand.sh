#!/bin/bash
#$ -l h_vmem=8G,h_rt=60:00:00,h_stack=32M
#$ -pe smp 8

source /.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate annotation_tutorial

cd /.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/v090226/active_line_analysis

# Bring in assembly and repeats

# Rscript --vanilla fifilter_repeat_long_line.R # Done on local desktop but can be done on cluster by changing directories

i=$I

# for i in pri alt ; do 

bedtools getfasta -fi /.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/v090226/mHetGla1.1.${i}.cleaned.fa -bed ${i}.2kbp.L1.bed -fo mHetGla1.1.${i}.L1.2Kb.fa -s 
# bedtools getfasta -fi /.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/v090226/mHetGla1.1.${i}.cleaned.fa -bed ${i}.2kbp.L1.neg.bed -fo mHetGla1.1.${i}.L1.2Kb.minus.fa -s

		
${PWD}/ORFfinder -in mHetGla1.1.${i}.L1.2Kb.fa -ml 300 -out mHetGla1.1.${i}.L1.2Kb.ORF.fa
# ${PWD}/ORFfinder -in mHetGla1.1.${i}.L1.2Kb.minus.fa -ml 300 -out mHetGla1.1.${i}.L1.2Kb.ORF.minus.fa -strand minus


