#!/bin/bash
#$ -l h_vmem=8G,h_rt=150:00:00,h_stack=32M
#$ -pe smp 28



cd /.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/mHetGlaT2T_repeats

source /.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate tidecluster

faDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/mHetGlaT2T_repeats

##
# Liner to call tandem repeats on each haplotype using TideCluster
##

for i in pri alt ; do

	mkdir -p mHetGla1.1.${i}_tidecluster ; mHetGla1.1.${i}_tidecluster

	TideCluster.py run_all -f $faDir/mHetGla1.1.${i}.cleaned.fa -c 28 -pr mHetGla1.1.${i}

done
