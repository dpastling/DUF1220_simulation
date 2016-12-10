#!/usr/bin/env bash
#BSUB -J simulation[1-60]
#BSUB -e logs/simulate_read_lengths_%J.log
#BSUB -o logs/simulate_read_lengths_%J.out
#BSUB -R "select[mem>20] rusage[mem=20] span[hosts=1]"
#BSUB -q normal
#BSUB -P Sikela

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

base_coverage=15
fastq_folder=fastq/template
bed_ref=reference/hg38_all_regions_10Mb_merged.bed
GENOME=$HOME/genomes/bowtie2.2.5_indicies/hg38/hg38.fa

READ_LENGTH=()
COVERAGE=()
for i in 36 50 100 150 300 600
do
	for j in {1..10}
	do
		READ_LENGTH+=($i)
	done
done
for i in 36 50 100 150 300 600
do
	for j in {1..10}
	do
		COVERAGE+=($j)
	done
done

read_length=${READ_LENGTH[$(($LSB_JOBINDEX - 1))]}
copies=${COVERAGE[$(($LSB_JOBINDEX - 1))]}

code/simulate_reads.sh -b $base_coverage -c $copies -l $read_length $fastq_folder $bed_ref

prefix=$fastq_folder/template_${read_length}bp_${copies}x

cut -f 4-9 ${prefix}_short.bed > ${prefix}_1.bed
cut -f 10-15 ${prefix}_short.bed > ${prefix}_2.bed

code/simulation_spikein/split_NBPF_from_template.pl -i reference/DUF1220_canonical_domains_hg38_v2.2.bed -t $fastq_folder/template_${read_length}bp_${copies}x  -o $fastq_folder


for bed in $fastq_folder/*/*_${read_length}bp_${copies}x_1.bed
do
	prefix=`echo $bed | sed 's/_1.bed$//'`
	temp=$prefix
	for pair in 1 2
	do
	        prefix=${temp}_$pair
	        bedtools getfasta -name -s -fi $GENOME -bed $prefix.bed -fo $prefix.fa
	        code/fasta2fastq.pl -f $prefix.fa -e data/qualities_R${pair}_${read_length}bp.txt -o $prefix.fastq
	        #rm $prefix.bed
	        rm $prefix.fa
	        gzip -f $prefix.fastq
	done
done




