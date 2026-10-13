#!/bin/bash

# {params.outDir} {wildcards.hap} {params.maskDir} {input.fa} {output.fa}
outDir=$1
hap=$2
maskDir=$3
inputFa=$4
outputFa=$5
outputBed=$6
threads=$7

cd $outDir

RepeatMasker -pa ${threads} -dir ./$hap"_mask" -lib ${maskDir}/HetGlaHap1_rm2.0_23Oct_withTRC.fa -xsmall -s $inputFa  -a -gff -html

cp $outDir/$hap"_mask"/$hap.fa.masked $outputFa
samtools faidx $inputFa

size=$(awk -F'\t' '{sum += $2} END {print sum}' ${inputFa}.fai)

cd $hap"_mask"

export PERL5LIB=$CONDA_PREFIX/share/RepeatMasker:$CONDA_PREFIX/share/RepeatMasker/util:$PERL5LIB

# + "/assembly/${hap}.fa", 
calcDivergenceFromAlign.pl -s ${hap}.fa.align.divsum -a ${hap}.fa.div.align  ${hap}.fa.align
createRepeatLandscape.pl -div ${hap}.fa.align.divsum -g $size > $hap.landscape.html
buildSummary.pl ${hap}.fa.out> $hap.summary.txt

awk 'NR>3 {
	chr=$5;
	start=$6;
	end=$7;
	name=$10"_"$11; # repeat name + family
	score=$1;       # SW score
	strand=($9=="C" ? "-" : "+");
	print chr"\t"start"\t"end"\t"name"\t"score"\t"strand
}' ./$hap.fa.out > $hap.mask.bed

sort -k1,1 -k2,2n -k3,3n $hap.mask.bed > $hap.mask.s.bed

awk 'BEGIN{FS=OFS="\t"} { sub(/^[^_]*_/, "", $4); print }' $hap.mask.s.bed > $outputBed

samtools faidx $outputFa
