data=$1
db=$2

db_dir=$(dirname "$db")
db_file=$(basename "$db")
db_basename="${db_file%.*}"

if [[ "$db_file" == *.gz ]]; then
    echo "Descomprimiendo archivo $db_file..."
    echo "gunzip -c $db > $db_dir/${db_basename}.fasta"
    gunzip -c "$db" > "$db_dir/${db_basename}.fasta" 
    db_final="$db_dir/${db_basename}"
    echo "makeblastdb -in $db_final.fasta -dbtype nucl -out $db_final"
    makeblastdb -in "$db_final.fasta" -dbtype "nucl" -out $db_final
elif [[ "$db_file" == *.fasta ]]; then
    echo "Archivo .fasta detectado, no se necesita descompresión."
    db_final="$db" 
    makeblastdb -in "$db_final" -dbtype "nucl" -out ${db_final%.*}
else
    db_final="$db"
fi
echo "Data file: $data"
echo "Database: $db_final"
blastn -query "$data" -db "$db_final" -out blast_results.txt -outfmt "6 qseqid sseqid pident length evalue bitscore" -num_threads 8  -evalue 1e-5 -perc_identity 90

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
