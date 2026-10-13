#!/bin/bash
#$ -l h_vmem=8G,h_rt=48:00:00,h_stack=32M
#$ -pe smp 1

outputDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/v090226/annotations
inputDir=/scratch3/groups/simpsonlab/dsokolowski/mHetGla1.pri.v060226
prefix=mHetGla1.pri

#############
# Annotations for each assembly was completed using GenAnT, aligning the RNAseq and Isoseq data with Hisat and minimap2 respectively, before merging bam files across each tissue
# The input for GenAnT is the yaml files shown. This script denotes the GenAnT outputs that are continued to be used in subsequent analyses
#############

# Final GFF
mv $inputDir/full_annotation.geneSymbols.gff $outputDir/$prefix.geneSymbols.gff
# GFF without gene symbols, for manual gene symbols to add
mv $inputDir/ncRNA_analysis/full_annotation.gff $outputDir/$prefix.gff

# gene symbol table from each source
mv $inputDir/gene_symbol/gene_symbols_noCopies.tsv $outputDir/$prefix.gene_symbols_noCopies.tsv
mv $inputDir/gene_symbol/gene_symbols_full.tsv $outputDir/$prefix.gene_symbols_full.tsv

# Orthofinder reference-to-target
mv $inputDir/orthofinder/orthofinder_protein.tsv $outputDir/$prefix.orthofinder_protein.tsv

# Mikado gene scoring algorithm
mv $inputDir/transcript_selection/mikado_lenient.metrics.tsv $outputDir/$prefix.metrics.tsv
mv $inputDir/transcript_selection/mikado_lenient.scores.tsv $outputDir/$prefix.scores.tsv

# Input GFFs
mv $inputDir/transcript_selection/stringtie.gffread.gff $outputDir/$prefix.stringtie.gff

#Toga outs
mv $inputDir/toga_outr1/*_toga/inact_mut_data.txt $outputDir/$prefix.inact_mut_data.txt
mv $inputDir/toga_outr1/*_toga/loss_summ_data.tsv $outputDir/$prefix.loss_summ_data.tsv
