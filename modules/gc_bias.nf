
process gc_bias {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.22"
    publishDir "${params.outputDir}/stats/gc_bias"

    input:
        tuple val(library), path(bam), path(bai)
        val(genome_fa)
        val(genome_fai)
    output:
        tuple val(library), path("${library}.gc_metrics"), emit: for_agg
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions
        tuple val("${task.process}"), val('picard'), eval('picard CollectGcBiasMetrics --version 2>&1 | cut -f 2 -d ":"'), topic: versions

    script:
    def prefix = group ? "${library}.${group}" : "${library}"
    // The group rides in ACCUMULATION_LEVEL because GcBias.parse in ngs-aggregate_results already
    // reads that column and maps picard's own labels to 'all'. Anchoring on 'All Reads' leaves the
    // ACCUMULATION_LEVEL header line alone; \t is GNU sed, which is what the conda env provides.
    def relabel_group = group ? "sed -i 's/^All Reads\\t/${group}\\t/' ${prefix}.gc_metrics" : ''
    """
    """
}
