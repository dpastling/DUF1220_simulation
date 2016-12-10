#!/usr/bin/env bash


genome=
first_read=
second_read=
bam_file=
crop="36"

while getopts g:b:c: opt; do
  case $opt in
  g)
      genome=$OPTARG
      ;;
  b)
      bam_file=$OPTARG
      ;;
  c)
      crop=$OPTARG
      ;;
  esac
done

shift $((OPTIND - 1))

first_read=$1
second_read=$2

if [[ -z $second_read ]]
then
        fastq_parameter="--seq $first_read"
else
        fastq_parameter="--seq1 $first_read --seq2 $second_read"
fi

if [ ! -z $crop ] && [ "$crop" != "0" ]
then
	crop="--crop $crop "
else
	crop=""
fi

align_file=`echo $bam_file | sed 's/\.bam//'`

mrsfast \
--search $genome \
$fastq_parameter \
-e 2 \
$crop\
--t 12 \
--disable-nohits \
--seqcomp \
-o $align_file.sam

samtools view -Sbh $align_file.sam > $align_file.bam
rm $align_file.sam

