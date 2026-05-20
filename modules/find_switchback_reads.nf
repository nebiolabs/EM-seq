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
}
