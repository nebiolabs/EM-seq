process multiqc {
    label 'medium_cpu'
    conda "bioconda::multiqc=1.25"
    publishDir "${params.outputDir}", mode: 'copy'

    input:
        path('*')

    output:
        path("*multiqc_report.html"), emit: multiqc_report

    script:
    """
    cat <<-CONFIG > multiqc_config.yaml
    title: EM-seq Alignment Summary - ${params.flowcell}
    extra_fn_clean_exts:
        - '.md'
        - '_combined_fastp'
    custom_plot_config:
        picard_insert_size:
            xmax: 1000
    table_columns_placement:
        Samtools Stats:
            raw_total_sequences: 10
            reads_mapped_percent: 20
            reads_properly_paired_percent: 30
            reads_MQ0_percent: 35
        Picard:
            summed_median: 50
    table_columns_visible:
        Picard:
            PCT_PF_READS_ALIGNED: False
            summed_mean: False
        Samtools Stats:
            reads_mapped: False
            mapped_passed: False
            non-primary_alignments: False
            reads_MQ0_percent: True
        Samtools Flagstat:
            mapped_passed: False
        samtools_idxstats_always:
            - plasmid_puc19c
            - phage_lambda
        FastQC:
            percent_duplicates: False
            total_sequences: False
            avg_sequence_length: False
            percent_fails: False
            total_sequences: False
    disable_version_detection: true
CONFIG

    multiqc -ip . 
    """
}
