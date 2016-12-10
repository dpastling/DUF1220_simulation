#!/bin/bash
#BSUB -J align[1-10]
#BSUB -e logs/mrsfast_%J.log
#BSUB -o logs/mrsfast_%J.out
#BSUB -R "select[mem>10] rusage[mem=10] span[hosts=1]"
#BSUB -q normal
#BSUB -n 12
#BSUB -P Sikela

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

fastq_folder=fastq/replicates
results_folder=alignments/replicates/mrsfast

SAMPLES=( $fastq_folder/*100bp_*_1.fastq.gz )
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
sample=`echo $sample | perl -pe 's/^.+?\/([^\/]+?)_[12].fastq.gz/$1/'`
first_pair=$fastq_folder/${sample}_1.fastq.gz
second_pair=$fastq_folder/${sample}_2.fastq.gz

GENOME=$HOME/genomes/mrsfast/hg38/hg38.fasta

read_length=100

align_file=$results_folder/${sample}_single.bam

#for read_length in 36 100
#do

crop="-c 36"
if [ $read_length -lt 40 ]
then
	crop="-c 0"
fi

code/mrsfast.sh -g $GENOME -b $align_file $crop $first_pair

