#!/usr/bin/env bash
#BSUB -J a1ign[1-40]
#BSUB -e logs/bowtie_total_%J.log
#BSUB -o logs/bowtie_total_%J.out
#BSUB -R "select[mem>5] rusage[mem=5] span[hosts=1]"
#BSUB -q normal
#BSUB -n 12
#BSUB -P Sikela

set -o nounset -o pipefail -o errexit -x

genome=~astlingd/genomes/bowtie-1.1.2/hg38/hg38
fastq_folder=fastq/replicates
result_folder=alignments/replicates/total

SAMPLES=( $fastq_folder/*_1.fastq.gz )
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

sample=`echo $sample | perl -pe 's/^.+?\/([^\/]+?)_[12].fastq.gz/$1/'`
first_pair=$fastq_folder/${sample}_1.fastq
second_pair=$fastq_folder/${sample}_2.fastq
bam=$result_folder/${sample}_paired.bam

gunzip -c $first_pair.gz > $first_pair
gunzip -c $second_pair.gz > $second_pair

code/bowtie_total.sh -g $genome -b $bam $first_pair $second_pair


bam=$result_folder/${sample}_single.bam
code/bowtie_total.sh -g $genome -b $bam $first_pair


#bowtie -p 12 -X 2000 --best --strata --all -v 2 -S --un $result_folder/${sample}_paired_unaligned $genome \
#-1 $first_pair -2 $second_pair | samtools view -Sb - > $result_folder/${sample}_paired.bam

rm $first_pair
rm $second_pair

