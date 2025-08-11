

process multiqc {
    label 'medium_cpu'
    conda "bioconda::multiqc=1.25"
    publishDir "${params.outputDir}", mode: 'copy'

    input:
        tuple val(email), path('*')

    output:
        path("*multiqc_report.html"), emit: multiqc_report

    script:
    '''
    cat <<-CONFIG > multiqc_config.yaml
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
	tuple  val(email), val(library), path(fq_or_bam), path(_read2), val(fileType),
	       val(barcodes), path(nonconverted_counts_tsv), path(fastp),
               path(bam), path(bai),
               path(gc_metrics),
               path(idxstat),
               path(flagstat),
               path(fastqc_zip),
               path(insertsize_metrics),
               path(tasmanian),
               path(mbias),
               path(alignment_summary_metrics_txt)

    output:
        path('ngs_agg.*')

    script:
    """
    path_to_ngs_agg="${params.path_to_ngs_agg}${params.revision}/"

    # bc = barcode1 + barcode2 if exists.
    if echo ${barcodes} | grep -q "+"
    then
        bc=\$(echo ${barcodes} | tr -d "][" | awk -F"+" '{bc2=""; if (length(\$2)==length(\$1)) {bc2="--barcode2 "\$2}; print \$1" "bc2;}')
    else
        bc=\$(echo ${barcodes} | tr -d "][" | awk -F"-" '{bc2=""; if (length(\$2)==length(\$1)) {bc2="--barcode2 "\$2}; print \$1" "bc2;}')
    fi

    # Validate barcodes
    if [[ ! "${barcodes}" =~ ^[+\\-ACGT]+\$ ]]; then
        echo "Warning: Invalid barcode format: ${barcodes}" >&2
    fi

    unzip *fastqc.zip

    cat ${nonconverted_counts_tsv} | awk -v l=${library} '{print l"\t"\$0}' > ${library}.nonconverted_counts.for_agg.tsv

    deploy_revision=\$(get_capistrano_release.sh) # script extracts currently executing revision of emseq from cap deploy logs
    
    export RBENV_VERSION=\$(cat \${path_to_ngs_agg}/.ruby-version)
    RAILS_ENV=production \${path_to_ngs_agg}/bin/bundle exec \${path_to_ngs_agg}/aggregate_results.rb \
    --bam ${bam} \
    --bai ${bai} \
    --barcode1 \${bc} \
    --contact_email ${params.email} \
    --genome \$(basename ${params.path_to_genome_fasta}) \
    --gc ${gc_metrics} \
    --idx_stats ${idxstat} \
    --flagstat ${flagstat} \
    --nonconverted_read_counts ${library}.nonconverted_counts.for_agg.tsv \
    --combined_mbias_records ${mbias} \
    --fastqc *_fastqc/fastqc_data.txt \
    --insert ${insertsize_metrics} \
    --tasmanian ${tasmanian} \
    --aln ${alignment_summary_metrics_txt} \
    --fastp ${fastp} \
    --metadata_bam_file ${fq_or_bam} \
    --workflow "${params.workflow} \${deploy_revision}" 2> ngs_agg.${library}.err 1> ngs_agg.${library}.out
    """
}
