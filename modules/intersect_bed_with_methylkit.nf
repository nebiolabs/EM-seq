process intersect_bed_with_methylkit {
    label 'process_single'
    tag "${library}"
    conda "bioconda::bedtools=2.31.1"

    input:
        tuple val(library), path(methylkit_bed)
        path(target_bed_prepared)
        val(genome_fa)
        val(genome_fai)

    output:
        tuple val(library), path('*_intersect.tsv'), emit: tsv
        tuple val("${task.process}"), val('bedtools'), eval('bedtools --version | sed \'s/^bedtools //\''), topic: versions

    script:
    """
    methylkit_basename=\$(basename "${methylkit_bed}" .bed)
    target_basename=\$(basename "${target_bed_prepared}" _slop_sorted.bed)

    # Uses bedtools intersect with -wa -wb to keep both annotations
    # nonamecheck since control contigs don't conform
    bedtools intersect -nonamecheck \\
        -a "${target_bed_prepared}" \\
        -b "${methylkit_bed}" \\
        -wa -wb -sorted -g "${genome_fai}" > \${methylkit_basename}_vs_\${target_basename}_intersect.tsv
    """
}
