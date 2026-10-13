#!/bin/bash
#$ -l h_vmem=16G,h_rt=150:00:00,h_stack=32M
#$ -pe smp 16

cd ${outDir}

# ==============================================================================
# summarize_sr_assemblies.sh
#
# Purpose: Per-assembly QC/report-card pipeline for candidate NMR (naked mole-rat)
#          genome assemblies. Designed to be run in parallel (via qsub) across
#          many candidate assembly FASTAs to compare completeness, contiguity,
#          annotation transfer, and repeat content before selecting a final
#          assembly.
#
# Usage:   qsub -N <name> -P simpsonlab -cwd -V -v I=<assembly.fa> summarize_sr_assemblies.sh
#          (typically looped: for i in G*fa ; do qsub ... -v I=$i ... ; done)
#
# Pipeline steps:
#   1. Compleasm      - BUSCO-like genome completeness check (rodentia lineage)
#   2. QUAST          - contiguity/assembly quality metrics (N50, length, etc.)
#   3. Liftoff        - transfers gene annotations from the mHetGla1.1.pri T2T
#                       reference onto the new assembly (builds/reuses a
#                       gffutils annotation DB)
#   4. gffread+BUSCO  - extracts lifted-over protein sequences and re-checks
#                       completeness at the protein level (rodentia_odb12)
#   5. RepeatMasker   - masks repeats using a custom NMR repeat library, then
#                       computes divergence, builds a repeat landscape plot,
#                       summarizes repeat content by family, and outputs a
#                       sorted, cleaned BED file of repeat annotations
#                       (for downstream genome browser / JBrowse visualization)
#
# Outputs (per assembly, prefix = basename of input .fa):
#   - ${fa}_compleasm/          : genome-level BUSCO-style completeness
#   - ${fa}_quast/               : assembly contiguity stats
#   - liftoff.${fa}/${fa}.gff    : lifted-over gene annotations
#   - liftoff.${fa}/${fa}.faa    : predicted protein sequences
#   - liftoff.${fa}/${fa}_buscoFaa/ : protein-level completeness check
#   - ${fa}_mask/${prefix}.landscape.html : repeat divergence landscape
#   - ${fa}_mask/${prefix}.summary.txt    : repeat content summary
#   - ${fa}_mask/${prefix}.repeat.bed     : sorted repeat annotation BED
# ==============================================================================

# /.mounts/labs/simpsonlab/users/dsokolowski/projects/AMR_genome_assembly/tools

source /.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate NMR_T2T

cd /.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/PreviousAssemblyCompleasm

outDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/PreviousAssemblyCompleasm

sifDir="/.mounts/labs/simpsonlab/users/dsokolowski/projects/AMR_genome_assembly/tools"
compleasm_db="/.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/test_NMR_assemblies/local_tests/rodentia_db"
threads=16
maskDir="/.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/final/primary_alt/repeat_mask"
refLiftOff="HetGlaT2T"
refLiftOffFa="/.mounts/labs/simpsonlab/users/dsokolowski/projects/GenAnT/data/references/mHetGla1.1.pri/mHetGla1.1.pri.cleaned.fa" # directory in /data -- adding species can be done with scripts in /utils
refLiftOffGff="/.mounts/labs/simpsonlab/users/dsokolowski/projects/GenAnT/data/references/mHetGla1.1.pri/mHetGla1.1.pri.cleaned.gff" # directory in /data -- adding species can be done with scripts in /utils
busco_odb="/.mounts/labs/simpsonlab/users/dsokolowski/projects/AMR_genome_assembly/tools/rodentia_odb12"

fa=$I


# Run compleasm
cd ${outDir}

${sifDir}/compleasm_kit/compleasm.py run -a ${fa} -o ${fa}"_compleasm" -l rodentia -t $threads -L ${compleasm_db}


# Run Quast
cd ${outDir}

quast.py --large -o ${fa}"_quast" ${fa}

# Run liftOff
cd ${outDir}

mkdir -p liftoff.${fa} ; cd liftoff.${fa}

if [ -e "referencegff_db" ]; then
	echo "referencegff_db exists. Skipping."
else
	# Create the file
	
	gffutils-cli create ${refLiftOffGff} -o referencegff_db
fi

liftoff -db referencegff_db ../${fa} ${refLiftOffFa} -o ./${fa}.gff -u ./unmapped.${fa}.txt -copies

gffread ${fa}.gff -g ../${fa} -y ${fa}.faa

busco --offline -i ${fa}.faa -o ${fa}_buscoFaa -m protein --lineage_dataset ${busco_odb} --cpu 16 -f


# Repeatmasker 

cd ${outDir}

RepeatMasker -pa ${threads} -dir ./${fa}"_mask" -lib ${maskDir}/HetGlaHap1_rm2.0_23Oct_withTRC.fa -xsmall -s $fa  -a -gff

cd ${fa}"_mask" 

prefix=`basename $fa .fa`

cp ${fa}.masked ${prefix}.sm.fa
samtools faidx $fa

size=$(awk -F'\t' '{sum += $2} END {print sum}' ${fa}.fai)

export PERL5LIB=$CONDA_PREFIX/share/RepeatMasker:$CONDA_PREFIX/share/RepeatMasker/util:$PERL5LIB

calcDivergenceFromAlign.pl -s ${fa}.align.divsum -a ${fa}.div.align  ${fa}.align
createRepeatLandscape.pl -div ${fa}.align.divsum -g $size > ${prefix}.landscape.html
buildSummary.pl ${fa}.out > ${prefix}.summary.txt

awk 'NR>3 {
	chr=$5;
	start=$6;
	end=$7;
	name=$10"_"$11; # repeat name + family
	score=$1;       # SW score
	strand=($9=="C" ? "-" : "+");
	print chr"\t"start"\t"end"\t"name"\t"score"\t"strand
}' ${fa}.out > ${fa}.mask.bed

sort -k1,1 -k2,2n -k3,3n ${fa}.mask.bed > ${fa}.mask.s.bed

awk 'BEGIN{FS=OFS="\t"} { sub(/^[^_]*_/, "", $4); print }' ${fa}.mask.s.bed > ${prefix}.repeat.bed

conda deactivate source /.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate barrnap
f=`basename $fa .fa`

barrnap --kingdom euk $i --threads 16 > $f.rrna.gff
