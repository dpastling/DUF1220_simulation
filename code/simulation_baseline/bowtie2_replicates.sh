#!/usr/bin/env bash
#BSUB -J align[1-40]
#BSUB -e logs/bowtie2_%J.log
#BSUB -o logs/bowtie2_%J.out
#BSUB -R "select[mem>5] rusage[mem=5] span[hosts=1]"
#BSUB -q normal
#BSUB -n 12
#BSUB -P Sikela

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

genome=~astlingd/genomes/bowtie2.2.5_indicies/hg38/hg38
fastq_folder=fastq/replicates
results_folder=alignments/replicates/normal

SAMPLES=( $fastq_folder/*_1.fastq.gz )
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

sample=`echo $sample | perl -pe 's/^.+?\/([^\/]+?)_[12].fastq.gz/$1/'`
first_pair=$fastq_folder/${sample}_1.fastq.gz
second_pair=$fastq_folder/${sample}_2.fastq.gz
bam=$results_folder/${sample}_paired.bam

code/bowtie2.sh -g $genome -b $bam $first_pair $second_pair

#bowtie2 -p 12 --very-sensitive --minins 0 --maxins 2000 -x $genome -1 $first_pair -2 $second_pair | samtools view -Sb - > $results_folder/${sample}_paired.bam



