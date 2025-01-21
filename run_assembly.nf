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
        file 'without_junk_rmdup.fastq.gz'

    publishDir { "results/${file.name.replace('.fastq.gz','')}/quality_control" }, mode: 'copy'
    
    script:
    """
    junk_dupes_removal.sh $filtered_fastq $report_postfilter
    """
}

process decontamination {
    debug true
    stageInMode 'copy'
    container 'community.wave.seqera.io/library/blast_nanoq_seqkit:cac86f202ab3ebdb'

    input:
        path file
        path db
        val db_name
        path cleaned_fasta

    output:
        file 'blast_cleaned.fasta'
        file 'cleaned_nanoq.fasta'
    
    publishDir { "results/${file.name.replace('.fastq.gz','')}/quality_control" }, mode: 'copy'

    script:
    """
    decontamination.sh $cleaned_fasta $db/$db_name
    """
}

process assembling {
    debug true
    container 'community.wave.seqera.io/library/flye:2.9.5--d577924c8416ccd8'

    input:
        path file
        path blast_cleaned_nanoq

    output:
        file 'assembly.fasta'

    publishDir { "results/${file.name.replace('.fastq.gz','')}/quality_control" }, mode: 'copy'
    
    script:
    """
    flye --nano-raw $blast_cleaned_nanoq --out-dir . --threads 8
    """
}

process polishing {
    debug true
    container 'community.wave.seqera.io/library/medaka:2.0.1--c15f6748e3c63d63'

    input:
        path file
        path blast_cleaned_fasta
        path flye_assembly

    output:
        file 'consensus.fasta'

    publishDir { "results/${file.name.replace('.fastq.gz','')}/quality_control" }, mode: 'copy'

    script:
    """
    medaka_consensus -o . -t 8 -i $blast_cleaned_fasta -d $flye_assembly
    """
}


process completness {
    debug true
    container 'community.wave.seqera.io/library/checkm2:1.0.2--b1ae1fcb2b8700f6'

    input:
        path file
        path medaka_consensus

    output:
        file 'checkm2_output/quality_report.tsv'

    publishDir { "results/${file.name.replace('.fastq.gz','')}/quality_control" }, mode: 'copy'

    script:
    """
    checkm2 database --download
    checkm2 predict --threads 30 --input $medaka_consensus --output-directory checkm2_output
    """
}

process depth {
    debug true
    container 'community.wave.seqera.io/library/minimap2_samtools:03e1e7cf6ec6695d'

    input:
        path file
        path junk_removed_fastq
        path medaka_consensus

    output:
        file 'output.depth'
    
    publishDir { "results/${file.name.replace('.fastq.gz','')}/quality_control" }, mode: 'copy'

    script:
    """
    minimap2 -ax map-ont $medaka_consensus $junk_removed_fastq > aln.sam

    samtools view aln.sam --threads 4 -o seq.bam --bam
    samtools sort seq.bam -o seq.sorted.bam
    samtools depth -a seq.sorted.bam -o output.depth
    """

}

process pyPlot {
    debug true
    container 'community.wave.seqera.io/library/python_pip_matplotlib_pandas_seaborn:ccd2d76e6689ee92'
    
    input:
        path file
        path depth_file

    output:
        file 'depth.png'

    publishDir { "results/${file.name.replace('.fastq.gz','')}/quality_control" }, mode: 'copy'

    script:
    """
    #!/usr/bin/env python
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns
    depth_info = pd.read_csv("output.depth", sep='\t', header=None, names=['Contig', 'Posicion', 'Profundidad'])
    depth_info["Global_Pos"] = depth_info.index + 1
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=depth_info, x='Global_Pos', y='Profundidad')
    plt.xlabel("Posición Genómica (bp)")
    plt.ylabel("Profundidad de Cobertura")
    plt.gca().set_xticklabels(['{:,}'.format(int(x)) for x in plt.gca().get_xticks() / 1e6])
    plt.grid(axis='y')
    plt.savefig("depth.png")
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

    decontamination(file_ch, db_ch, db_name_ch,junkDupesRemoval.out[0])
    
    assembling(file_ch, decontamination.out[1])

    polishing(file_ch, decontamination.out[0], assembling.out[0])

    completness(file_ch, polishing.out[0])

    depth(file_ch, junkDupesRemoval.out[1], polishing.out[0])

    pyPlot(file_ch, depth.out[0])
}
