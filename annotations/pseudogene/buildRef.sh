#!/bin/bash
#$ -l h_vmem=24G,h_rt=20:00:00,h_stack=32M
#$ -pe smp 8


source /.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate annotation_tutorial

cd /.mounts/labs/simpsonlab/users/dsokolowski/projects/GenAnT/data/references/mHetGla1.1.pri

wd=/.mounts/labs/simpsonlab/users/dsokolowski/projects/GenAnT/data/references/mHetGla1.1.pri
GenAnT=/.mounts/labs/simpsonlab/users/dsokolowski/projects/GenAnT
fa=mHetGla1.1.pri.cleaned.fa # GCF_000001635.27_GRCm39_genomic.fna
gff=mHetGla1.1.pri.cleaned.gff # GCF_000001635.27_GRCm39_genomic.gff

####
# script to set up the reference genome directory, reference for GenAnT and also the structure used to call pseudogenes
####

bash $GenAnT/setup/reference_directory_mikado.sh $wd $GenAnT $fa $gff

