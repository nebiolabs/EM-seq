process picard_metrics {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.22"
    publishDir "${params.outputDir}/stats/picard_alignment_metrics"

    input:
        tuple val(library), path(bam), path(bai)
        val(genome_fa)
        val(genome_fai)

    output:
        tuple val(library), path("${library}.alignment_summary_metrics.txt"), emit: for_agg
        tuple val("${task.process}"), val('picard'), eval('picard CollectAlignmentSummaryMetrics --version 2>&1 | cut -f 2 -d ":"'), topic: versions

    script:
    """
    picard -Xmx${task.memory.toGiga()}g CollectAlignmentSummaryMetrics \
        --VALIDATION_STRINGENCY SILENT -BS true -R ${genome_fa} \
        -I ${bam} -O ${library}.alignment_summary_metrics.txt
    """
}
