process aggregate_emseq {
    tag { library }
    conda "bioconda::samtools=1.9"

    input:         
	tuple  val(library), 
            path(bam),
            path(aligned_bam), 
            path(aligned_bam_bai), 
            path(fastp),
            path(markdups_log),
            path(alignment_summary_metrics_txt),
            path(gc_metrics),
            path(idxstat),
            path(flagstat),
            path(fastqc_zip),
            path(nonconverted_counts_tsv), 
            path(tasmanian),
            path(mbias),
            path(insertsize_metrics)

    output:
        path('ngs_agg.*')

    script:
    """
    unzip -o *fastqc.zip # -f in case we need to re-run

    cat ${nonconverted_counts_tsv} | awk -v l=${library} '{print l"\t"\$0}' > ${library}.nonconverted_counts.for_agg.tsv

    export RBENV_VERSION=\$(cat ${params.path_to_ngs_agg}/.ruby-version)
    RAILS_ENV=production ${params.path_to_ngs_agg}/bin/bundle exec ${params.path_to_ngs_agg}/aggregate_results.rb \\
    --metadata_bam_file ${bam} \\
    --bam ${aligned_bam} \\
    --bai ${aligned_bam_bai} \\
    --contact_email ${params.email} \\
    --fastp ${fastp} \\
    --aln ${alignment_summary_metrics_txt} \\
    --dup ${markdups_log} \\
    --gc ${gc_metrics} \\
    --idx_stats ${idxstat} \\
    --flagstat ${flagstat} \\
    --fastqc *_fastqc/fastqc_data.txt \\
    --nonconverted_read_counts ${library}.nonconverted_counts.for_agg.tsv \\
    --combined_mbias_records ${mbias} \\
    --tasmanian ${tasmanian} \\
    --insert ${insertsize_metrics} \\
    --workflow ${params.workflow} 2> ngs_agg.${library}.err 1> ngs_agg.${library}.out
    """
}
