#!/bin/bash
#BSUB -J align[1-271]
#BSUB -e logs/bowtie2_%J.log
#BSUB -o logs/bowtie2_%J.out
#BSUB -R "select[mem>10] rusage[mem=10] span[hosts=1]"
#BSUB -q normal
#BSUB -n 12
#BSUB -P Sikela

# 271 genes

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

genome=~astlingd/genomes/bowtie2.2.5_indicies/hg38/hg38
fastq_folder=fastq/template
results=alignments/spikein/normal

source code/config_domains.sh
domain=${DOMAINS[$(($LSB_JOBINDEX - 1))]}

for sample in $fastq_folder/*/${domain}_*_1.fastq.gz
do
        sample=`echo $sample | perl -pe 's/^.+?\/([^\/]+?)_[12].fastq.gz/$1/'`
        gene=`echo $sample | perl -pe 's/^(NBPF[^_]+).+?$/$1/'`
        if [[ ! -d $results/$gene ]]; then
                mkdir -p $results/$gene
        fi

        first_pair=$fastq_folder/$gene/${sample}_1.fastq.gz
        second_pair=$fastq_folder/$gene/${sample}_2.fastq.gz
	
	code/bowtie2.sh -g $genome -b $results/$gene/${sample}_paired.bam $first_pair $second_pair

	code/bowtie2.sh -g $genome -b $results/$gene/${sample}_single.bam $first_pair

#	bowtie2 -p 12 --very-sensitive --minins 0 --maxins 2000 -x $genome -1 $first_pair -2 $second_pair | samtools view -Sb - > alignments/$folder/$gene/${sample}_paired.bam

#	bowtie2 -p 12 --very-sensitive --minins 0 --maxins 2000 -x $genome -U $first_pair | samtools view -Sb - > alignments/$folder/$gene/${sample}_single.bam
done



