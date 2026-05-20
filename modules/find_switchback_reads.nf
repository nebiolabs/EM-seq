process find_switchback_reads {
    label 'low_cpu'
    tag { library }
    publishDir "${params.outputDir}/stats/switchbacks"
    conda "bioconda::fgbio"

    input:
        tuple val(library), path(bam), path(bai)
        val(genome_fa)

    output:
        tuple val(library), path("${library}.switchbacks.summary.txt"), emit: for_agg
        tuple val(library), path("${library}.switchbacks.bam")
        tuple val("${task.process}"), val('fgbio'), eval('fgbio --version 2>&1 | head -n 1'), topic: versions

    script:
    """
    fgbio -Xmx${task.memory.toGiga()}g FindSwitchbackReads \\
        --input ${bam} \\
        --output ${library}.switchbacks.bam \\
        --metrics ${library}.switchbacks \\
        --ref ${genome_fa}
    """

    stub:
    """
    touch ${library}.switchbacks.bam
    printf 'sample\tlibrary\ttemplates\taligned_templates\tswitchback_templates\tfraction_switchbacks\tread_based_switchbacks\tmean_length\tmean_offset\ttandem_based_switchbacks\tmean_gap\n' > ${library}.switchbacks.summary.txt
    printf '${library}\t${library}\t0\t0\t0\t0.0\t0\t0.0\t0.0\t0\t0.0\n' >> ${library}.switchbacks.summary.txt
    """
}
