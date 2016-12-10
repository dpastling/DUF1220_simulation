#!/usr/bin/env bash
#BSUB -J coverage[1-40]
#BSUB -e logs/coverage_%J.log
#BSUB -o logs/coverage_%J.out
#BSUB -R "select[mem>10] rusage[mem=10] span[hosts=1]"
#BSUB -q normal
#BSUB -P Sikela

set -o nounset -o pipefail -o errexit -x

#folder=replicates/mrsfast
#folder=replicates/multiread
#folder=replicates/normal
folder=replicates/total

#multiread_alignments=alignments/replicates/total
multiread_alignments=alignments/replicates/total-strata

results_folder=results/$folder
reference_bed=reference/DUF1220_canonical_domains_hg38_v2.2.bed
alignments_folder=alignments/$folder
multiread=

if [[ $results_folder =~ "multiread" ]]
then
        multiread="-m"
        alignments_folder=$multiread_alignments
fi

SAMPLES=( $alignments_folder/*.bam )
bam=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
sample=`echo $bam | perl -pe 's/^.+?\/([^\/]+?)\.bam/$1/'`

file=$results_folder/$sample

if [[ ! -d $results_folder ]]; then
	mkdir -p $results_folder
fi

pairing=
if [[ $sample =~ "paired" ]]
then
    pairing="paired"
elif [[ $sample =~ "single" ]]
then
    pairing="single"
else
    echo "not able to determine if single or paired-end reads"
    exit 1
fi

source code/make_bed.sh -r $reference_bed -p $pairing $multiread $bam $file



