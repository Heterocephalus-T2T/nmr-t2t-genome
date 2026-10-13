#!/bin/bash

#SBATCH -t 200:0:0
#SBATCH --mem=100000
#SBATCH -N 1
#SBATCH -c 50

source /hpf/largeprojects/mdwilson/external/miniconda3/etc/profile.d/conda.sh

#module purge

##################################################
##
### Liner to call topologically associated domains
##
##################################################

cd /home/dsokolo3/mdwilson/dustin/NMR_general/ABC_enhancer/HiC

module load juicer/1.6-encode-CPU
#module load bwa
module load bwa/0.7.17
module load coreutils
module load java
#module load samtools
module load samtools/1.9

module load bedtools
conda activate /home/dsokolo3/mdwilson/dustin/NMR_general/ABC_enhancer/ABC/

HiC=/home/dsokolo3/mdwilson/dustin/NMR_general/ABC_enhancer/HiC

java -jar juicer_tools_1.19.02.jar arrowhead --cpu --threads 48 -m 1024 -c "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29,chrX" -k KR mHetGlaPri_hypothalamus.hic mHetGlaPri_hypothalamus_arrowhead
