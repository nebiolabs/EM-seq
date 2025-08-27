process concatenate_intersections {
    label 'process_single'
    publishDir "${params.outputDir}/methylKit_intersections", mode: 'copy'

    input:
        path(intersection_files)
        path(summary_files)

    output:
        path('all_intersections_combined.tsv'), emit: combined_intersections
        path('all_summaries_combined.tsv'), emit: combined_summaries

    script:
    """
    # Combine all intersection files
    echo -e "methylkit_file\\tchr\\tstart\\tend\\tcontext\\tmethylation\\ttarget_locus\\ttarget_name" > all_intersections_combined.tsv

    # Add all intersection data (skip headers from individual files)
    for file in ${intersection_files}; do
        if [ -s "\$file" ]; then
            tail -n +2 "\$file" >> all_intersections_combined.tsv
        fi
    done

    # Combine all summary files
    echo -e "methylkit_file\\ttarget_length\\tposition\\tcontext\\tmean_methylation\\tn_loci\\tn_measurements" > all_summaries_combined.tsv

    # Add all summary data (skip headers from individual files)
    for file in ${summary_files}; do
        if [ -s "\$file" ]; then
            tail -n +2 "\$file" >> all_summaries_combined.tsv
        fi
    done

    # Report statistics
    total_intersections=\$(tail -n +2 all_intersections_combined.tsv | wc -l)
    total_summaries=\$(tail -n +2 all_summaries_combined.tsv | wc -l)

    echo "Combined \${total_intersections} intersection records from \$(echo ${intersection_files} | wc -w) files"
    echo "Combined \${total_summaries} summary records from \$(echo ${summary_files} | wc -w) files"
    """
}
