#!/usr/bin/env bash
#BSUB -J simulation[1-20]
#BSUB -e logs/trim_replicates_%J.log
#BSUB -o logs/trim_replicates_%J.out
#BSUB -R "select[mem>10] rusage[mem=10] span[hosts=1]"
#BSUB -q normal
#BSUB -P Sikela

fastq_folder=fastq/replicates
result_folder=fastq/replicates

for file in $fastq_folder/template_100bp_*.fastq.gz
do
	fastq_file=`echo $file | perl -pe 's/^.+?\/([^\/]+)\.gz$/$1/'`
	fastq_file=$result_folder/$fastq_file
	trimmed_file=`echo $fastq_file | sed 's/100bp/trim/'`
	gunzip -c $file > $fastq_file
	~astlingd/code/utilities/randomlySampleFastq.pl -l 36 -o $trimmed_file $fastq_file
	gzip -f $trimmed_file
	rm $fastq_file
done


