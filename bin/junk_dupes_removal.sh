input=${1}
report_postfilter=${2}

# Remove junk_seq in header
seqkit grep -n -v -r -p "junk_seq" ${input} -o pseudomonas-aeruginosa_without_junk.fastq.gz
# Remove duplicated sequences
seqkit rmdup -s pseudomonas-aeruginosa_without_junk.fastq.gz -o pseudomonas-aeruginosa_without_junk_rmdup.fastq.gz
# Convert fastq to fasta
seqkit fq2fa pseudomonas-aeruginosa_without_junk_rmdup.fastq.gz -o cleaned.fasta
# Run nanoq
nanoq -i cleaned.fasta -o cleaned_nanoq.fasta -vvv -r post_clean_report.txt

# Get the number of reads before and after removing contaminated reads
total_before=$(grep "Number of reads" ${report_postfilter} | awk '{print $4}')
total_after=$(grep "Number of reads" post_clean_report.txt | awk '{print $4}')

echo "Number of reads before removing dups and junk: $total_before"
echo "Number of reads after removing dups and junk: $total_after"

percentage=$(echo "scale=2; $total_after * 100 / $total_before" | bc)
echo "Percentage of reads kept: $percentage%"
