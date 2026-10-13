#!/bin/bash
#$ -l h_vmem=12G,h_rt=230:00:00,h_stack=32M
#$ -pe smp 8

# /.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/P-GRe

source /.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate annotation_tutorial

cd /.mounts/labs/simpsonlab/users/dsokolowski/projects/GenAnT/data/references/mHetGla1.1.pri

genomeDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/GenAnT/data/references/mHetGla1.1.pri
pGreDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/annotation_pipeline/P-GRe
prefix=mHetGla1.1.pri


####
# Uses our genome annotation and assembly + mmus pseuogene and gene structures to call pseudogenes with PGre
####

${pGreDir}/PGRe.sh -f ${genomeDir}/${prefix}.cleaned.fa \
		-g ${genomeDir}/${prefix}.cleaned.gff \
		-p ${genomeDir}/${prefix}.cleaned.protein.faa \
		-o ${genomeDir} \
		-d ${pGreDir}
		-u /.mounts/labs/simpsonlab/users/dsokolowski/projects/GenAnT/data/references/mmus_GRC39/GCF_000001635.27_GRCm39_genomic.protein.faa \
		-d ${pGreDir} -A -t 8

