#!/usr/bin/env bash
#BSUB -J coverage[1-271]
#BSUB -e logs/coverage_%J.log
#BSUB -o logs/coverage_%J.out
#BSUB -R "select[mem>10] rusage[mem=10] span[hosts=1]"
#BSUB -q normal
#BSUB -P Sikela

# catch unset variables, non-zero exits in pipes and calls, enable x-trace.
set -o nounset -o pipefail -o errexit -x

source code/config_domains.sh
domain=${DOMAINS[$(($LSB_JOBINDEX - 1))]}

folder=spikein/normal
#folder=spikein/mrsfast
#folder=spikein/total
#folder=spikein/multiread

results=results/$folder
reference_bed=reference/DUF1220_canonical_domains_hg38_v2.2.bed
alignments_folder=alignments/$folder
multiread=

if [[ $results =~ multiread ]]
then
	multiread="-m"
	alignments_folder=alignments/spikein/total-strata
fi


for bam in $alignments_folder/*/${domain}_*.bam
do
    sample=`echo $bam | perl -pe 's/^.+?\/([^\/]+?)\.bam/$1/'`
    gene=`echo $sample | perl -pe 's/^(NBPF[^_]+).+?$/$1/'`
    file=$results/$gene/$sample
    if [[ ! -d $results/$gene ]]; then
        mkdir -p $results/$gene
    fi

    pairing=
    if [[ $sample =~ paired ]]
    then
    	pairing="paired"
    else
    	pairing="single"
    fi

    code/make_bed.sh -r $reference_bed -p $pairing $multiread $bam $file
done


