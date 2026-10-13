#!/bin/bash

outDir=$1
sifDir=$2
compleasm_odb=$3
fa=$4
wildcard=$5
threads=$6


cd ${outDir}

# /.mounts/labs/simpsonlab/users/dsokolowski/projects/AMR_genome_assembly/tools

${sifDir}/compleasm_kit/compleasm.py run -a ${fa} -o ${wildcard}"_compleasm" -l rodentia -t $threads -L ${compleasm_odb}

