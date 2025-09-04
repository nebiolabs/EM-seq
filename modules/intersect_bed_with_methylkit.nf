process intersect_bed_with_methylkit {
    label 'process_single'
    tag "${library}"
    conda "bioconda::bedtools=2.31.1"
    publishDir "${params.outputDir}/methylKit_intersections", pattern: "*summary.tsv"

    input:
        tuple val(library), path(methylkit_bed)
        path(target_bed_prepared)
        val(genome_fa)
        val(genome_fai)

    output:
        tuple val(library), path('*_intersect.tsv'), emit: tsv
        tuple val(library), path('*_region_summary.tsv'), emit: region_summary
        tuple val("${task.process}"), val('bedtools'), eval('bedtools --version | sed \'s/^bedtools //\''), topic: versions

    script:
    """
    target_basename=\$(basename "${target_bed_prepared}" _slop_sorted.bed)

    # nonamecheck since control contigs don't conform
    # save the intersection output and also summarize mean methylation per region/context
    bedtools intersect -nonamecheck \\
        -a "${target_bed_prepared}" \\
        -b "${methylkit_bed}" \\
        -wa -wb -sorted -g "${genome_fai}" | \\
    tee ${library}_vs_\${target_basename}_intersect.tsv | \\
    awk -v FS='\t' -v OFS='\t' 'NR>1 {print \$2,\$3,\$4,\$7,\$6, \$5}' | sort -k4,4 -k6,6 | \\
    bedtools groupby -g 4,6 -o mean -c 5 > ${library}_vs_\${target_basename}_region_summary.tsv
    """
}
