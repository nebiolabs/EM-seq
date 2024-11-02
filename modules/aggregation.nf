
process multiqc {
    tag { library }
    conda "bioconda::multiqc=1.25"
    publishDir "${params.outputDir}"

    input:
        tuple env(flowcell), val(email), val(library), val(barcodes), 
              path(nonconverted_counts_tsv), path(fastp), path(bam), path(bai), 
              path(gc_metrics), path(idxstat), path(flagstat), path(fastqc_zip), 
              path(insertsize_metrics), path(tasmanian), path(mbias), 
              path(alignment_summary_metrics_txt)

    output:
        path("multiqc_report.html"), emit: multiqc_report

    script:
    '''
    cat <<CONFIG > multiqc_config.yaml 
    title: EM-seq Alignment Summary - !{flowcell}
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
    CONFIG

    multiqc -ip . 
    '''
}

process aggregate_emseq {
    tag { library }
    conda "bioconda::samtools=1.9"
    publishDir "${params.outputDir}/ngs-agg"

    input:
        tuple env(flowcell), val(email), val(library), val(barcodes), 
              path(nonconverted_counts_tsv), path(fastp), path(bam), path(bai), 
              path(gc_metrics), path(idxstat), path(flagstat), path(fastqc_zip), 
              path(insertsize_metrics), path(tasmanian), path(mbias), 
              path(alignment_summary_metrics_txt)

    output:
        path('ngs_agg.*')

    shell:
    '''

    genome_name=$(echo !{params.genome} | awk -F"/" '{print $NF}' | sed 's/.fa|.fasta//')

    unzip *fastqc.zip

    cat !{nonconverted_counts_tsv} | awk -v l=!{library} '{print l"\t"$0}' > !{library}.nonconverted_counts.for_agg.tsv
    
    export RBENV_VERSION=$(cat !{params.path_to_ngs_agg}/.ruby-version)
    RAILS_ENV=production !{params.path_to_ngs_agg}/bin/bundle exec !{params.path_to_ngs_agg}/aggregate_results.rb \
        --bam !{bam} \
        --bai !{bai} \
        --name !{library} \
        --lane !{params.lane} \
        --contact_email !{params.email} \
        --project !{params.project} \
        --sample !{params.sample} \
        --genome !{params.genome} \
        --gc !{gc_metrics} \
        --idx_stats !{idxstat} \
        --flagstat !{flagstat} \
        --nonconverted_read_counts !{library}.nonconverted_counts.for_agg.tsv \
        --combined_mbias_records !{mbias} \
        --fastqc *_fastqc/fastqc_data.txt \
        --insert !{insertsize_metrics} \
        --tasmanian !{tasmanian} \
        --aln !{alignment_summary_metrics_txt} \
        --metadatafq_file !{metadata_fastq} \
        --fastp !{fastp} \
        --workflow !{params.workflow}_!{library} 2> ngs_agg.err 1> ngs_agg.out
    '''
}
