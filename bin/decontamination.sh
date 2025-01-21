data=$1
db=${2:-"db/echerichia_db"}

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
echo "blast_$(basename $data .fasta).fasta"
seqkit grep -v -f blast_results_ids.txt $data -o blast_$(basename $data)

echo "Running Nanoq with blast_$(basename $data)"
nanoq -i blast_$(basename $data) -o $(basename $data .fasta)_nanoq.fasta -vvv -r nanoq_report.txt
echo "Saving into $(basename $data .fasta)_nanoq.fasta"
