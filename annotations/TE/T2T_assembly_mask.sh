#!/bin/bash
#$ -l h_vmem=8G,h_rt=60:00:00,h_stack=32M
#$ -pe smp 16

source /.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate annotation_tutorial

cd /.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/mHetGlaT2T_repeats

# Bring in assembly and repeats

##
# Run and analyze repeat masker with a custom curated library
## 

fa=$I

# fa=`basename $i .fa`

threads=16
maskDir="/.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/final/primary_alt/repeat_mask"

# RepeatMasker -pa ${threads} -dir ./${fa}"_mask" -lib ${maskDir}/HetGlaHap1_rm2.0_23Oct_withTRC.fa -xsmall -s $fa  -a -gff

prefix=`basename $fa .fa`

# cp ${fa}.masked ${prefix}.sm.fa
samtools faidx $fa

size=$(awk -F'\t' '{sum += $2} END {print sum}' ${fa}.fai)

export PERL5LIB=$CONDA_PREFIX/share/RepeatMasker:$CONDA_PREFIX/share/RepeatMasker/util:$PERL5LIB

cd ${fa}"_mask"

# calcDivergenceFromAlign.pl -s ${fa}.align.divsum -a ${fa}.div.align  ${fa}.align
createRepeatLandscape.pl -div ${fa}.align.divsum -g $size > ../${prefix}.landscape.html
# buildSummary.pl ${fa}.out > ../${prefix}.summary.txt
