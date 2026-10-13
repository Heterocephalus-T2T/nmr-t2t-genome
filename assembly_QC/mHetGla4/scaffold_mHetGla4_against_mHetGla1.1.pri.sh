#!/bin/bash
#$ -l h_vmem=24G,h_rt=20:00:00,h_stack=32M
#$ -pe smp 4



source  /.mounts/labs/simpsonlab/users/dsokolowski/miniconda3/bin/activate RagTag

cd /.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/TOR_assemblies

###########################
##
### This script uses RagTag and to orient mHetGla4.1.pri to mHetGla1.1.pri. It does not use patch, so regions that were not assemble-able without the CAM ultras are left as gaps.
##
###########################


i=$I

refFA=/.mounts/labs/simpsonlab/users/dsokolowski/projects/NMR_genome_assembly/final/primary_alt/mHetGla1.1.pri.cleaned.fa
f=`basename $i .p_ctg.fa`

cut -f 1-2 $i.fai > $i.sizes
(echo "$i.idx"; cat "$i.sizes") > $i.idx

mkdir -p $f

minimap2 -t 4 -x asm5  $refFA $i > $f/$f.paf

ragtag.py correct -t 4 -o $f $refFA $i

ragtag.py scaffold -t 4 -o $f -r -m 20000000 $refFA $f/ragtag.correct.fasta
# ragtag.py patch -t 4 -o $f $refFA $f/ragtag.scaffold.fasta

mv $f/ragtag.correct.fasta $f/$f.correct.fasta
mv $f/ragtag.scaffold.fasta pafs/$f.scaffold.fasta
# mv $f/ragtag.patch.fasta $f/$f.patch.fasta


