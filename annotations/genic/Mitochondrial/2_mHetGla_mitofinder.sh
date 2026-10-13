#!/bin/bash
#$ -l h_vmem=8G,h_rt=48:00:00,h_stack=32M
#$ -pe smp 16

cd /.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/TOR_assemblies/MT_assembly

sifDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/AMR_genome_assembly/tools

contigDir=/.mounts/labs/simpsonlab/production/mole_rat/assembly/HetGla.tor.6add52

module load singularity

mkdir -p /scratch2/groups/simpsonlab/dsokolowski/tor_hifi

hifiDir=/scratch2/groups/simpsonlab/dsokolowski/tor_hifi

SINGULARITY_CACHEDIR=$sifDir/cachedir
SINGULARITY_TMPDIR=$sifDir/tmpdir

mkdir -p $sifDir/cachedir
mkdir -p $sifDir/tmpdir
refDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/TOR_assemblies/MT_assembly/MT_reference


mkdir -p mfinder ; cd mfinder

cp ../mHetGla1.1.pri.chrMT.fa ./
# cat /.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/basecall/fastq/hifi/*gz > $hifiDir/torhifi.fastq.gz

#########
##
### This script uses mitofinder to annotate the mHetGla mitochondrial genome
##
#########

singularity exec --bind ${PWD},${refDir},${hifiDir},${sifDir} $sifDir/mitohifi.sif mitofinder -j mfinder -a ${PWD}/mHetGla1.1.pri.chrMT.fa -o 2 -r $refDir/NC_015112.1.gb -p  16

agat sp gxf2gxf -g mfinder/mfinder_MitoFinder_mitfi_Final_Results -o mfinder_mtDNA_contig.clean.gff

# If you need to do it from the assembly its
# 1- run mitohifi on contigs to get mt assembly and genbank version
# 2- rerun mitofinder because mitohifi throws out gtf file and gb --> gtf conversion is messy (https://github.com/marcelauliano/MitoHiFi/issues/126)
	# this is this script
# 3- clean with agat ( this is this script )
