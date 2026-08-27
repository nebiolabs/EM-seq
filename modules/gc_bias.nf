// One curve per invocation. With an empty group this is the whole-reference curve; with a group it
// is that organism's curve, and picard_reference MUST then be the group's subset reference, never
// the composite
//
// Invoke it twice in one workflow under an alias (Nextflow forbids calling a process twice by the
// same name).
process gc_bias {
    tag { "${library}${group ? ':' + group : ''}" }
    label 'medium_cpu'
    conda "bioconda::picard=3.3.0"
    publishDir "${params.outputDir}/stats/gc_bias"

    input:
        tuple val(library), path(bam), path(bai), val(group), val(picard_reference)

    output:
        tuple val(library), path("${library}${group ? '.' + group : ''}.gc_metrics"), emit: for_agg
        tuple val("${task.process}"), val('picard'), eval('picard CollectGcBiasMetrics --version 2>&1 | cut -f 2 -d ":"'), topic: versions

    script:
    def prefix = group ? "${library}.${group}" : "${library}"
    // The group rides in ACCUMULATION_LEVEL because GcBias.parse in ngs-aggregate_results already
    // reads that column and maps picard's own labels to 'all'. Anchoring on 'All Reads' leaves the
    // ACCUMULATION_LEVEL header line alone; \t is GNU sed, which is what the conda env provides.
    def relabel_group = group ? "sed -i 's/^All Reads\\t/${group}\\t/' ${prefix}.gc_metrics" : ''
    """
    picard -Xmx${task.memory.toGiga()}g CollectGcBiasMetrics \\
        --IS_BISULFITE_SEQUENCED true --VALIDATION_STRINGENCY SILENT \\
        -I ${bam} -O ${prefix}.gc_metrics -S ${prefix}.gc_summary_metrics \\
        --CHART /dev/null -R ${picard_reference}
    ${relabel_group}
    """
}
