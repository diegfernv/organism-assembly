data=$1
db=${2:-"db/escherichia"}

echo "Data file: $data"
echo "Database: $db"
blastn -query "$data" -db "$db" -out blast_results.txt -outfmt "6 qseqid sseqid pident length evalue bitscore" -num_threads 8  -evalue 1e-5 -perc_identity 90

echo "Getting ids from blast results"
# Get ids from results
ids=$(awk '{print $1}' "blast_results.txt")
# Delete duplicated ids
ids=$(echo "$ids" | sort | uniq)

# Count number of ids
echo "Number of ids to be deleted: $(echo "$ids" | wc -l)"

# Save ids to file
echo "$ids" > blast_results_ids.txt

# Remove contaminated reads into new file
echo "$(basename $data .fasta)_clean.fasta"
seqkit grep -v -f blast_results_ids.txt $data -o blast_$(basename $data)

