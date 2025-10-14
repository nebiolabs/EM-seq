process calculate_feature_methylation {
    label 'medium_cpu'
    tag { "${feature_type}" }
    conda "bedtools=2.29.2 htslib=1.9"
    publishDir "${params.outputDir}/feature_methylation", mode: 'copy'

    input:
        tuple val(feature_type), path(feature_gtf)
        path(methylation_bed)

    output:
        tuple val(feature_type), path("${feature_type}_methylation.tsv"), emit: methylation
        tuple val("${task.process}"), val('bedtools'), eval('bedtools --version | sed \'s/^bedtools v//\''), topic: versions

    script:
    """
    bedtools intersect -nonamecheck \
    -wa -wb -loj \
    -a ${feature_gtf} -b <(bgzip -d < ${methylation_bed}) \
    | awk -v FS='\\t' -v OFS='\\t' '\$14>0 {print \$10,\$11,\$12,\$1":"\$4-1"-"\$5,(\$15*1.0)/\$14 }' \
    | bedtools groupby -g 4 -o mean -c 5 \
    > ${feature_type}_methylation.tsv
    """
}

process combine_methylation {
    label 'low_cpu'
    publishDir "${params.outputDir}", mode: 'copy'

    input:
        path(methylation_files)

    output:
        path('combined_methylation.tsv'), emit: combined

    script:
    """
    echo 'File\tLocus\tFrac Methylated' > combined_methylation.tsv

    for f in ${methylation_files}; do
        filebase=\$(basename "\${f}" _methylation.tsv)
        lines=\$(wc -l < <(grep -ve '^\\s*\$' -e '^#' "\$f") | cut -f 1 -d ' ')
        if [ \$lines -gt 0 ]; then
            paste <(yes \${filebase} | head -n \$lines) <(grep -ve '^\\s*\$' -e '^#' "\$f") >> combined_methylation.tsv
        fi
    done
    """
}
