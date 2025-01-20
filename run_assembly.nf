#!/usr/bin/env nextflow

process qualityControl {
    debug true
    container 'community.wave.seqera.io/library/blast_nanoq_seqkit:cac86f202ab3ebdb'

    input:
        path file
        val min_length
        val max_length
        val min_quality
        val coverage
        val organism_size

    output:
        file 'filtered.fastq.gz'
        file 'report_postfilter.txt'

    publishDir { "results/${file.name.replace('.fastq.gz','')}/quality_control" }, mode: 'copy'

    script:
    """
    run_filter.sh $file $min_length $max_length $min_quality $coverage $organism_size "report_prefilter.txt" "report_postfilter.txt"
    """
}

process junkDupesRemoval {
    debug true
    container 'community.wave.seqera.io/library/blast_nanoq_seqkit:cac86f202ab3ebdb'

    input:
        path file
        path filtered_fastq
        path report_postfilter

    output:
        file 'cleaned.fasta'

    publishDir { "results/${file.name.replace('.fastq.gz','')}/quality_control" }, mode: 'copy'
    
    script:
    """
    junk_dupes_removal.sh $filtered_fastq $report_postfilter
    """
}

process decontamination {
    debug true
    container 'community.wave.seqera.io/library/blast_nanoq_seqkit:cac86f202ab3ebdb'

    input:
        path db
        val db_name
        path cleaned_fasta

    script:
    """
    blastn -query $cleaned_fasta -db $db/echerichia -out blast_results.txt -outfmt "6 qseqid sseqid pident length evalue bitscore" -num_threads 8  -evalue 1e-5 -perc_identity 90
    """
}

workflow {
    file_ch = Channel.fromPath(params.input_file)
    db_ch = Channel.fromPath("db")
    db_name_ch = Channel.of(params.db_name)
    min_len_ch = Channel.of(params.min_length)
    max_len_ch = Channel.of(params.max_length)
    min_quality_ch = Channel.of(params.min_quality)
    coverage_ch = Channel.of(params.coverage)
    organism_size_ch = Channel.of(params.organism_size)
    
    log.info("Using min_length=${params.min_length}, max_length=${params.max_length}, min_quality=${params.min_quality}")
    qualityControl(file_ch, min_len_ch, max_len_ch, min_quality_ch, coverage_ch, organism_size_ch)
    
    junkDupesRemoval(file_ch, qualityControl.out[0], qualityControl.out[1])

    decontamination(db_ch, db_name_ch,junkDupesRemoval.out[0])
}
