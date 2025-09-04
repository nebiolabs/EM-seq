process concatenate_files {
    label 'process_single'
    publishDir "${params.outputDir}/methylKit_intersections", mode: 'copy'

    input:
        path(files)
        val(header)
        val(output_prefix)

    output:
        path('*_combined.tsv'), emit: combined_file

    script:
    """
    echo -e "${header}" > ${output_prefix}_combined.tsv

    for file in ${files}; do
        if [ -s "\$file" ]; then
            tail -n +2 "\$file" >> ${output_prefix}_combined.tsv
        fi
    done
    """
}
