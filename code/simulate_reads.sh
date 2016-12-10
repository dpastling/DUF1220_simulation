#!/usr/bin/env bash

base_coverage=15
copies=2
read_length=100
replicate=
duf_ref=reference/DUF1220_canonical_domains_hg38_v2.2.bed

while getopts b:c:l:r: opt; do
  case $opt in
  b)
      base_coverage=$OPTARG
      ;;
  c)
      copies=$OPTARG
      ;;
  l)
      read_length=$OPTARG
      ;;
  r)
      replicate=$OPTARG
      ;;
  esac
done

shift $((OPTIND - 1))


fastq_folder=$1
reference_bed=$2

insert_size=350
if [[ $read_length -eq 300 ]]; then
        insert_size=700
fi
if [[ $read_length -eq 600 ]]; then
        insert_size=1400
fi


coverage=$(( $base_coverage * $copies ))
#number_of_reads=$( echo $total_size $coverage $insert_size | awk '{ printf "%.0f", ($1 * $2) / $3}' )

prefix=$fastq_folder/template_${read_length}bp_${copies}x
if [[ -n $replicate ]]
then
	prefix=${prefix}_$replicate
fi

sleep 2
#code/randomly_sample_bed.pl -n $number_of_reads -l $read_length --paired --insert $insert_size -o $prefix $reference_bed
code/randomly_sample_bed.pl -c $coverage -l $read_length --paired --insert $insert_size -o $prefix $reference_bed


read1=${prefix}_1.bed
read2=${prefix}_2.bed
paste $read1 $read2 > ${prefix}_merged.bed
awk 'OFS="\t" {min=$2;max=$3;if($8<min)min=$8;if($9>max)max=$9;print $1, min, max, $0}' ${prefix}_merged.bed > ${prefix}_append.bed
sort -k 1,1 -k 2,2n -T ./ ${prefix}_append.bed > ${prefix}_sorted.bed
rm ${prefix}_merged.bed
rm ${prefix}_append.bed

bedtools intersect -wa -sorted -a ${prefix}_sorted.bed -b $duf_ref > ${prefix}_short.bed

rm ${prefix}_sorted.bed
rm ${prefix}_1.bed
rm ${prefix}_2.bed



