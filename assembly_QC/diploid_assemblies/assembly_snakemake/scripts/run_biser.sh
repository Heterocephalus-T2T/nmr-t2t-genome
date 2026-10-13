#!/bin/bash

outDir=$1
hap=$2
threads=$3

mkdir -p ${outDir}/biser.${hap} ; cd ${outDir}/biser.${hap}

export PATH="/.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/envs/NMR_T2T/bin:$PATH" 

mkdir -p "./temp"

biser -o biser_${hap} -t ${threads} -T "temp" --keep-temp --keep-contigs --gc-heap 4G ${outDir}/assembly/${hap}.sm.fa
