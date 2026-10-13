#!/bin/bash
#SBATCH -t 20:0:0
#SBATCH --mem=32000
#SBATCH -N 1
#SBATCH -c 20

##################################################
##
### MACS2 peak calling of PE reads with effective NMR genome size
##
##################################################

module load MACS2/2.2.9.1

i=$I

set -e

cd /home/dsokolo3/mdwilson/dustin/NMR_new_genome/T2T_assemblies/CENPA_ChIPseq/${i}/bam
 
mkdir -p tmp1

macs2 callpeak -t *CENPA_liver*.bam -c *Input_liver*.bam -f BAMPE -n "CENPA_liver.${i}.bam.macs2" -g 1870000000 --broad --tempdir tmp1
