#!/bin/bash

outDir=$1
sifDir=$2
mmi=$3
threads=$4
hap=$5
fa=$6
fai=$7

cd ${outDir}


# minimap2 -t ${threads} -ax map-ont --MD ${mmi} ${outDir}/ont/ont.fq.gz | samtools view -@ ${threads} -b | samtools sort -@ ${threads} -o ${outDir}/ont/${hap}.ont.sort.bam

# samtools index ${outDir}/ont/${hap}.ont.sort.bam


mkdir -p ${hap}_ont_flagger ; cd ${hap}_ont_flagger

ln -s $fa ./
samtools faidx $hap.fa

cat ${hap}.fa.fai | \
    awk '{print $1"\t0\t"$2}' > whole_genome.bed

echo "{" > annotations_path.json
echo \"whole_genome\" : \"${PWD}/whole_genome.bed\" >> annotations_path.json
echo "}" >> annotations_path.json

date


# fullout=/.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/test_NMR_assemblies/local_tests/font
# mkdir -p $fullont

# config["outDir"]+"/ont/{hap}.ont.sort.bam.bai"

singularity exec --bind ${outDir} ${sifDir}/flagger.sif bam2cov \
        --bam ${outDir}/ont/${hap}.ont.sort.bam \
        --output ./coverage_file.cov.gz \
        --annotationJson ./annotations_path.json \
        --threads $threads \
        --baselineAnnotation whole_genome

date
echo FlaggerRun

cp /.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/assembly_summary/7_flagger/container/alpha_optimum_trunc_exp_gaussian_w_16000_n_50.ONT_R941_Guppy6.3.7_DEC_2024.v1.1.0.tsv alpha.tsv

singularity exec --bind ${outDir} ${sifDir}/flagger.sif hmm_flagger \
        --input coverage_file.cov.gz \
        --outputDir ${outDir}/${hap}_ont_flagger  \
        --alphaTsv ./alpha.tsv \
        --labelNames Err,Dup,Hap,Col \
        --threads $threads \
        --windowLen 16000
