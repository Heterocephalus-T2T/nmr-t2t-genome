#!/bin/bash

outDir=$1
sifDir=$2
mmi=$3
threads=$4

cd ${outDir}

# minimap2 -t ${threads} -ax map-hifi --MD ${mmi} ${outDir}/hifi/hifi.fq.gz | samtools view -@ ${threads} -b | samtools sort -@ ${threads} -o ${outDir}/hifi/diploid.sort.bam
# samtools index ${outDir}/hifi/diploid.sort.bam

mkdir -p flagger ; cd flagger

ln -s ${outDir}/assembly/diploid.fa ./
samtools faidx diploid.fa

cd $outDir/flagger

# Create a bed file covering whole genome

cat diploid.fa.fai | \
awk '{print $1"\t0\t"$2}' > whole_genome.bed

# Put the path to the whole-genome bed file in a json file
echo "{" > annotations_path.json
echo \"whole_genome\" : \"${PWD}/whole_genome.bed\" >> annotations_path.json
echo "}" >> annotations_path.json


singularity exec --bind ${outDir} $sifDir/flagger.sif bam2cov \
	--bam ${outDir}/hifi/diploid.sort.bam \
	--output ./coverage_file.cov.gz \
	--annotationJson ./annotations_path.json \
	--threads ${threads} \
	--baselineAnnotation whole_genome

date 
echo FlaggerRun

cp /.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/assembly_summary/7_flagger/container/alpha_optimum_trunc_exp_gaussian_w_16000_n_50.HiFi_DC_1.2_DEC_2024.v1.1.0.tsv alpha.tsv

singularity exec --bind ${outDir} ${sifDir}/flagger.sif hmm_flagger \
	--input coverage_file.cov.gz \
	--outputDir ${outDir}/flagger  \
	--alphaTsv ./alpha.tsv \
	--labelNames Err,Dup,Hap,Col \
	--threads ${threads} \
	--windowLen 16000
