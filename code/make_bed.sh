#!/usr/bin/env bash

input_folder=
output_folder=
bam=
multiread=
pairing="paired"
reference_bed=

while getopts r:p:m opt; do
  case $opt in
  r)
      reference_bed=$OPTARG
      ;;
  p)
      pairing=$OPTARG
      ;;
  m)
      multiread=1
      ;;
  esac
done

shift $((OPTIND - 1))

bam=$1
output=$2

if [ "$pairing" == "paired" ]
then
    samtools sort -n -@ 5 -m 5G $bam ${output}_sorted
    bedtools bamtobed -split -bedpe -i ${output}_sorted.bam > $output.bed
    code/parse_bed.pl $output.bed
elif [ "$pairing" == "single" ]
then
    bedtools bamtobed -i $bam > ${output}_edited.bed
else
    echo "Unknown pairing scheme"
    exit
fi
sort -k 1,1 -k 2,2n -T ./ ${output}_edited.bed > ${output}_sorted.bed
bedtools intersect -wao -sorted -a $reference_bed -b ${output}_sorted.bed > ${output}_temp.bed
if [[ -n $multiread ]]
then
    code/multi_read_correction.pl ${output}_temp.bed > ${output}_corrected.bed
    bedtools merge -scores sum -i ${output}_corrected.bed > ${output}_coverage.bed
else
    awk 'OFS="\t" {print $4,$2,$3,$1,$13}' ${output}_temp.bed | bedtools merge -scores sum -i - > ${output}_coverage.bed
fi
awk 'OFS="\t" { print $1, $4 / ($3 - $2 + 1)}' ${output}_coverage.bed > ${output}_read_depth.bed
#if [ -f $output.bed ]
#then
#    rm ${output}_sorted.bam
#    rm $output.bed
#fi
#rm ${output}_edited.bed
#rm ${output}_temp.bed


