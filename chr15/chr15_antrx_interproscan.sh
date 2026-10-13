#!/bin/bash
#$ -l h_vmem=8G,h_rt=20:00:00,h_stack=32M
#$ -pe smp 8

source /.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate annotation_tutorial

cd /.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/v090226/chr15/antxrl


i=$I

# Run interproscan on a faa file, in this case it's for the antrx orthogenes

externalDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/GenAnT/external

export PATH="$externalDir/my_interproscan/interproscan-5.69-101.0:$PATH"

interproDir=${externalDir}/my_interproscan/interproscan-5.69-101.0

mkdir -p interpro.${i} ; cd interpro.${i} 

${interproDir}/interproscan.sh -i ../${i} -f tsv -dp --cpu 8
