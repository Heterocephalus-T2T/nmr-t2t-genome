#!/bin/bash
#SBATCH -t 200:0:0
#SBATCH --mem=100000
#SBATCH -N 1
#SBATCH -c 50
#SBATCH --tmp=1000000

############################
# Alignment and post-processing of ChIP-seq data, using parameters for centromere annotations in Arora, Sullivan and Dumont (2023, cel rep. PMID: 37742188)
############################


module load bedtools
module load bwa
module load samtools
module load python/3.6.3

j=$J
i=$I

f=`basename $i _paired_R1.fastq.gz`

inputDir=/home/dsokolo3/mdwilson/dustin/NMR_new_genome/T2T_assemblies/CENPA_ChIPseq # directory with files are
outputDir=${inputDir}/${j}

genomeDir=/home/dsokolo3/mdwilson/genomes/${j}/bwa_index

# cd $SCRATCH

cd ${outputDir}

# Alignment + sorting
bwa mem \
	-k 19 \
	-w 100 \
	-d 100 \
	-r 1.5 \
	-c 10000 \
	-A 1 \
	-B 4 \
	-O 6 \
	-E 1 \
	-L 5 \
	-U 9 \
	-T 30 \
	-v 3 \
	-t 48 \
  ${genomeDir}/*.fa \
  ${inputDir}/$f"_paired_R1.fastq.gz" ${inputDir}/$f"_paired_R2.fastq.gz" \
| samtools view -@ 48 -b - | samtools sort -@ 48 -o ${f}.sorted.bam

# Index the sorted BAM
samtools index ${f}.sorted.bam

module load python/3.6.7_deepTools

mkdir -p $SCRATCH/${f}
export TMPDIR=$SCRATCH/${f}

bamCoverage -b ${f}.sorted.bam -o ${f}.bw --normalizeUsing CPM --binSize 10 --centerReads

mv ${f}.sorted.bam ${outputDir}/bam
mv ${f}.sorted.bam.bai ${outputDir}/bam
mv ${f}.bw ${outputDir}/bigwig
