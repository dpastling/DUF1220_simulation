#!/bin/bash
#BSUB -J a1ign[1-271]
#BSUB -e logs/bowtie_total_%J.log
#BSUB -o logs/bowtie_total_%J.out
#BSUB -R "select[mem>10] span[hosts=1]"
#BSUB -q normal
#BSUB -n 12
#BSUB -P Sikela

set -o nounset -o pipefail -o errexit -x

#source code/config_genes.sh
#gene=${GENES[$(($LSB_JOBINDEX - 1))]}
genome=~astlingd/genomes/bowtie-1.1.2/hg38/hg38
folder=template
result_folder=alignments/spikein/total-strata

source code/config_domains.sh
domain=${DOMAINS[$(($LSB_JOBINDEX - 1))]}

for sample in fastq/$folder/*/${domain}_*_1.fastq.gz
#for sample in fastq/$folder/${gene}_*_1.fa.gz
do
        sample=`echo $sample | perl -pe 's/^.+?\/([^\/]+?)_[12].fastq.gz/$1/'`
	gene=`echo $sample | perl -pe 's/^(NBPF[^_]+).+?$/$1/'`
        if [[ ! -d $result_folder/$gene ]]; then
                mkdir -p $result_folder/$gene
        fi
	
        first_pair=fastq/$folder/$gene/${sample}_1.fastq
        second_pair=fastq/$folder/$gene/${sample}_2.fastq

	gunzip -c $first_pair.gz > $first_pair
	gunzip -c $second_pair.gz > $second_pair
	
	code/bowtie_multiread.sh -g $genome -b $result_folder/$gene/${sample}_paired.bam $first_pair $second_pair

	code/bowtie_multiread.sh -g $genome -b $result_folder/$gene/${sample}_single.bam $first_pair

	#bowtie -p 12 -X 2000 --best --strata --all -v 2 -S --un $result_folder/${sample}_paired_unaligned $genome -1 $first_pair -2 $second_pair | samtools view -Sb - > $result_folder/${sample}_paired.bam
	#bowtie -p 12 -X 2000 --best --strata --all -v 2 -S --un $result_folder/${sample}_single_unaligned $genome $first_pair | samtools view -Sb - > $result_folder/${sample}_single.bam

	rm $first_pair
	rm $second_pair
done

