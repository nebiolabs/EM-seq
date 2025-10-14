process calculate_feature_counts {
    label 'high_cpu'
    tag { feature_type }
    conda "subread=2.1.1"
    publishDir "${params.outputDir}/feature_counts", mode: 'copy'

    input:
        tuple val(feature_type), path(feature_file)
        path(bam_with_index)
        val(feature_id_attr) // e.g., 'transcript_id', 'gene_id'

    output:
        tuple val(feature_type), path("${feature_type}_counts.tsv"), emit: counts
        tuple val("${task.process}"), val('featureCounts'), eval('featureCounts -v 2>&1 | grep "featureCounts" | sed \'s/featureCounts v//\''), topic: versions

    script:
    """
    featureCounts --primary --ignoreDup -Q 10 -M -f -O --fraction -p -P -B -C \
    -a ${feature_file} -F GTF \
    -t ${feature_type} \
    -g ${feature_id_attr} \
    --tmpDir ${params.tmp_dir ?: '/tmp'} \
    -T ${task.cpus} \
    -o ${feature_type}_counts.tsv *.bam
    """
}

process combine_counts {
    label 'low_cpu'
    publishDir "${params.outputDir}", mode: 'copy'

    input:
        path(count_files)

    output:
        path('combined_feature_counts.tsv'), emit: combined

    script:
    """
    # Construct the header
    echo -n 'File\t' > combined_feature_counts.tsv
    grep -hve '^\\s*\$' -e '^#' ${count_files[0]} | head -n 1 >> combined_feature_counts.tsv

    # Add a column containing the name of the file being processed
    for f in ${count_files}; do
        filebase=\$(basename "\${f}" _counts.tsv)
        lines=\$(wc -l < <(grep -ve '^\\s*\$' -e '^#' "\$f") | cut -f 1 -d ' ')
        if [ \$lines -gt 1 ]; then
            paste <(yes \${filebase} | head -n \$lines) <(grep -ve '^\\s*\$' -e '^#' "\$f") | tail -n +2 >> combined_feature_counts.tsv
        fi
    done
    """
}
