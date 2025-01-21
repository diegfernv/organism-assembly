#!/bin/bash

pre_filter_file=${7}
post_filter_file=${8}
file=${1}
min_length=${2}
max_length=${3}
min_quality=${4}
coverage=${5}
organism_size=${6}

# Run pre-filtering nanoq
nanoq -i $file -o "reads.fastq.gz" -vvv -r ${pre_filter_file}
nbases_pre=$(grep "Number of bases" ${pre_filter_file} | awk '{print $4}')
echo "Number of bases before filtering: $nbases_pre"

nanoq -i $file -o filtered.fastq.gz -l $min_length -m $max_length -q $min_quality -vvv -r ${post_filter_file}
