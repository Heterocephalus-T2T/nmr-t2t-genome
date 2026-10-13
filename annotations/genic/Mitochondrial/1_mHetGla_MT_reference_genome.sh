#!/bin/bash
#$ -l h_vmem=8G,h_rt=60:00:00,h_stack=32M
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
refDir=/.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/TOR_assemblies/MT_assembly
mkdir -p no_Contig ; cd no_Contig
#########
##
### This script uses mitohifi/mitofinder to pull the existing naked mole-rat mitochrondiral genome as a reference to build around. Internet is required
##
#########
singularity exec --bind ${PWD},${refDir},${hifiDir},${sifDir},${contigDir} $sifDir/mitohifi.sif findMitoReference.py --species "Heterocephalus glaber" --outfolder ./MT_reference --min_length 14000

