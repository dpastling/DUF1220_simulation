#!/usr/bin/env bash
#BSUB -J simulation[1-271]
#BSUB -e logs/trim_spikein_%J.log
#BSUB -o logs/trim_spikein_%J.out
#BSUB -R "select[mem>10] rusage[mem=10] span[hosts=1]"
#BSUB -q normal
#BSUB -P Sikela

set -o nounset -o pipefail -o errexit -x

source code/config_domains.sh
domain=${DOMAINS[$(($LSB_JOBINDEX - 1))]}
folder=fastq/template

for file in $folder/*/${domain}_100bp_*.fastq.gz
do
	fastq_file=`echo $file | sed 's/.gz//'`
	trimmed_file=`echo $fastq_file | sed 's/100bp/100bptrim/'`
	gunzip -c $file > $fastq_file
	code/shortenReadLength.pl -l 36 -o $trimmed_file $fastq_file
	gzip -f $trimmed_file
	rm $fastq_file
done

