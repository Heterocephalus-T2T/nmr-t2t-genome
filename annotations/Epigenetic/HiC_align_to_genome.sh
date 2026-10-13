#!/bin/bash

#SBATCH -t 200:0:0
#SBATCH --mem=100000
#SBATCH -N 1
#SBATCH -c 50
#SBATCH --tmp=7000000

module load juicer/1.6-encode-CPU
module load bwa/0.7.17
module load coreutils
module load java
module load samtools/1.9

##################################################
##
### Juicer scripts to align and process Hi-C data against mHetGla1.1.pri
##
##################################################


assembly=/home/dsokolo3/mdwilson/genomes/mHetGla1.1.pri/bwa_index/mHetGla1.1.pri.cleaned.fa
outsDir=/home/dsokolo3/mdwilson/dustin/NMR_general/ABC_enhancer/HiC

HiCDir="/home/dsokolo3/mdwilson/dustin/NMR_new_genome/Hypo_HiC/"

cd $SCRATCH

mkdir fastq

cd fastq 

ln -s "/home/dsokolo3/mdwilson/dustin/NMR_new_genome/Hypo_HiC/HiC1/"*.gz ./
ln -s "/home/dsokolo3/mdwilson/dustin/NMR_new_genome/Hypo_HiC/HiC2/"*.gz ./

cat WL2949_none_hypothal_NA_none_hglaber_NA_Basal_TCAG1_R1_001.fastq.gz WL2950_none_hypothal_NA_none_hglaber_NA_Basal_TCAG1_R1_001.fastq.gz > NMR_HiC_R1.fastq.gz
cat WL2949_none_hypothal_NA_none_hglaber_NA_Basal_TCAG1_R2_001.fastq.gz WL2950_none_hypothal_NA_none_hglaber_NA_Basal_TCAG1_R2_001.fastq.gz > NMR_HiC_R2.fastq.gz

unlink WL2949_none_hypothal_NA_none_hglaber_NA_Basal_TCAG1_R1_001.fastq.gz
unlink WL2950_none_hypothal_NA_none_hglaber_NA_Basal_TCAG1_R1_001.fastq.gz
unlink WL2949_none_hypothal_NA_none_hglaber_NA_Basal_TCAG1_R2_001.fastq.gz
unlink WL2950_none_hypothal_NA_none_hglaber_NA_Basal_TCAG1_R2_001.fastq.gz

cd $SCRATCH

juicer.sh -d $SCRATCH -z $assembly -a 'Fibroblast cellline alignment to longread NMR_male genome' -p $assembly -D "/hpf/tools/centos7/juicer/2.0-encode-CPU/"

cd aligned
samtools view -O SAM -F 1024 merged_dedup.sam | awk -v mnd=1 -f /hpf/tools/centos7/juicer/1.6-encode-CPU/CPU/common/sam_to_pre.awk > merged_nodups.txt
java -Xmx64g -jar /home/dsokolo3/mdwilson/dustin/NMR_general/ABC_enhancer/HiC/juicer_tools_1.19.02.jar pre merged_nodups.txt "maternal_ktype.hic" /home/dsokolo3/mdwilson/genomes/mHetGla1.1.pri/chrom.sizes
mv maternal_ktype.hic $outsDir"/mHetGlaPri_hypothalamus.hic"


