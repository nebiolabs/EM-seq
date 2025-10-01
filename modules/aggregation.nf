

process multiqc {
    label 'medium_cpu'
    conda "bioconda::multiqc=1.25"
    publishDir "${params.outputDir}"

    input:
        path('*')

    output:
        path("*multiqc_report.html"), emit: multiqc_report

    script:
    """
    cat <<-CONFIG > multiqc_config.yaml
    title: EM-seq Alignment Summary - ${params.flowcell}
    extra_fn_clean_exts:
        - '.aln'
        - '.md'
        - '.fastp.json'
        - '.tmp'
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
    """
}

process aggregate_emseq {
    tag { library }
    conda "bioconda::samtools=1.9"

    input:
	tuple   val(library),
            path(aligned_bam), path(aligned_bam_bai),
            val(barcodes), val(flowcell), val(num_reads_used),
            path(fastp),
            path(nonconverted_counts_tsv),
            path(markdups_log),
            path(gc_metrics),
            path(idxstat),
            path(flagstat),
            path(fastqc_zip),
            path(mismatches),
            path(mbias),
            path(alignment_summary_metrics_txt),
            path(insertsize_metrics),
            path(file_metadata)

    output:
        path('ngs_agg.*')

    script:
    """
    path_to_ngs_agg="${params.path_to_ngs_agg}${params.revision}/"

    unzip -o *fastqc.zip

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

    echo "task attempt is ${task.attempt}"

    bc1=\$(echo \$bc | sed 's/ \\-\\-barcode2.*//')
    bc2=\$(echo \$bc | sed 's/.*\\-\\-barcode2 //')

    if [ ${task.attempt} -eq 2 ] || [ ${task.attempt} -eq 4 ]; then
        bc2=\$(echo "\$bc2" | rev | tr "[ATCG]" "[TAGC]")
    elif [ ${task.attempt} -eq 3 ] || [ ${task.attempt} -eq 4 ]; then
        bc1=\$(echo "\$bc1" | rev | tr "[ATCG]" "[TAGC]")
    fi

    if echo ${file_metadata} | grep -q "fastq\$\\|fq\$"; then
        metadata="--metadatafq_file ${file_metadata}"
    elif echo ${file_metadata} | grep -q "bam\$" ; then
        metadata="--metadata_bam_file ${file_metadata}"
    else
        echo "Error: Unsupported file type in metadata file: ${file_metadata}" >&2
        exit 1
    fi

    cat ${nonconverted_counts_tsv} | awk -v l=${library} '{print l"\t"\$0}' > ${library}.nonconverted_counts.for_agg.tsv
    export RBENV_VERSION=\$(cat \${path_to_ngs_agg}/.ruby-version)
    RAILS_ENV=production \${path_to_ngs_agg}/bin/bundle exec \${path_to_ngs_agg}/aggregate_results.rb \
    --bam ${aligned_bam} \
    --bai ${aligned_bam_bai} \
    --contact_email ${params.email} \
    --genome ${params.path_to_genome_fasta} \
    --barcode1 \${bc1} \
    --barcode2 \${bc2} \
    --flowcell ${flowcell} \
    --gc ${gc_metrics} \
    --idx_stats ${idxstat} \
    --flagstat ${flagstat} \
    --nonconverted_read_counts ${library}.nonconverted_counts.for_agg.tsv \
    --combined_mbias_records ${mbias} \
    --fastqc *_fastqc/fastqc_data.txt \
    --insert ${insertsize_metrics} \
    --tasmanian ${mismatches} \
    --aln ${alignment_summary_metrics_txt} \
    --fastp ${fastp} \
    \${metadata} \
    --workflow ${params.workflow} 2> ngs_agg.${library}.err 1> ngs_agg.${library}.out
    """
}
