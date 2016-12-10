#!/usr/bin/env bash
#BSUB -J align[1-271]
#BSUB -e logs/mrsfast_%J.log
#BSUB -o logs/mrsfast_%J.out
#BSUB -R "select[mem>20] rusage[mem=20] span[hosts=1]"
#BSUB -q normal
#BSUB -n 12
#BSUB -P Sikela

# 22 genes

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

#source code/config_genes.sh
#gene=${GENES[$(($LSB_JOBINDEX - 1))]}

source code/config_domains.sh
domain=${DOMAINS[$(($LSB_JOBINDEX - 1))]}
gene=`echo $domain | perl -pe 's/^(NBPF[^_]+)_.+?$/$1/'`

GENOME=$HOME/genomes/mrsfast/hg38/hg38.fasta
fastq_folder=fastq/template
result_folder=alignments/spikein/mrsfast

for read_length in 36 100
do
#	for sample in $fastq_folder/${gene}/${gene}_*_${read_length}bp_*_1.fastq.gz
	for sample in $fastq_folder/$gene/${domain}_${read_length}bp_*_1.fastq.gz
	do
		sample=`echo $sample | perl -pe 's/^.+?\/([^\/]+?)_1.fastq.gz/$1/'`
		first_pair=$fastq_folder/$gene/${sample}_1.fastq.gz
		second_pair=$fastq_folder/$gene/${sample}_2.fastq.gz
		align_file=$result_folder/$gene/${sample}_single.bam

		if [ ! -d "$result_folder/$gene" ]
		then
			mkdir $result_folder/$gene
		fi

		crop="-c 36"
		if [ $read_length -lt 40 ]
		then
			crop="-c 0"
		fi

		code/mrsfast.sh -g $GENOME -b $align_file $crop $first_pair 

	done
done

