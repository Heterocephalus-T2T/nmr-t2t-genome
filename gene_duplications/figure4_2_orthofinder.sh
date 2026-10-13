#!/bin/bash
#$ -l h_vmem=16G,h_rt=150:00:00,h_stack=32M
#$ -pe smp 16

cd /.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/Figure4_Orthofinder

i=$I

cd ${i}

source /.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate annotation_tutorial

orthofinder -t 16 -a 16 -o $i"_orthofinder_v090326" -f protein_seqs
