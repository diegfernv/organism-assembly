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

## Read filtered bp count
nbases_post=$(grep "Number of bases" ${post_filter_file} | awk '{print $4}')
echo "Number of bases after filtering: $nbases_post"

# Calculate required bp for $coverage coverage
required_bp=$(echo "scale=0; $coverage * $organism_size" | bc)

# Is coverage reached?
if [ $nbases_post -ge $required_bp ]; then
    echo "Coverage reached"
else
    echo "Coverage not reached"
fi
