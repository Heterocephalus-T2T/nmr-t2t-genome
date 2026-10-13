#!/bin/bash
#$ -l h_vmem=16G,h_rt=150:00:00,h_stack=32M
#$ -pe smp 8

# cd /.mounts/labs/simpsonlab/users/dsokolowski/projects/AMR_genome_assembly/Daniel_M/HelEmini/HelEmi_repeat_libary
cd /.mounts/labs/simpsonlab/users/dsokolowski/projects/AMR_genome_assembly/Bathyergus/BatSui_repeat_library

export sourceDir=/.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate
source $sourceDir NMR_T2T


##
# Combine libraries from bultiple sources in order of priority. Clean redundant repeat elements with cd-hit-est, where the tidecluster called repeat gets kept > repeatmodeler > annotatios from NMR
## 

cat BatSui_tidecluster.fa BatSui_repeatmodeler.fa HetGlaHap1_rm2.0_23Oct_withTRC.fa > BatSui_combined.fa


cd-hit-est \
  -i BatSui_combined.fa \
  -o BatSui_combined.nr.fa \
  -c 0.9 \
  -n 8 \
  -M 0 \
  -T 8

