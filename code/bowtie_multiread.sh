#!/usr/bin/env bash

genome=
first_read=
second_read=
bam_file=
max_insert=2000
fastq_parameter=

while getopts g:b:i: opt; do
  case $opt in
  g)
      genome=$OPTARG
      ;;
  b)
      bam_file=$OPTARG
      ;;
  i)
      max_insert=$OPTARG
      ;;
  esac
done

shift $((OPTIND - 1))

if [ ! -n "$1" ] || [ ! -n "$genome" ]
then
  echo "Usage: `basename $0` [options] -g genome_index -b output_file <read_1.fastq> <read_2.fastq if paired>"
  echo "-----------------------------------------------------------------------------------------------------"
  echo "            -i max_insert_size  (default: 2000)" 
  exit 1
fi

first_read=$1
if (( $# >= 2 ))
then
	second_read=$2
	fastq_parameter="-1 $first_read -2 $second_read"
else
	fastq_parameter="$first_read"
fi

bowtie \
  -p 12 \
  -X $max_insert \
  --best \
  --strata \
  --all \
  -v 2 \
  -S \
  $genome \
  $fastq_parameter \
  | samtools view -Sb - > $bam_file


