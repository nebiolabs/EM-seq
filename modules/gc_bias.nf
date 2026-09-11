// The whole-reference GC bias curve for one library.
//
// Per-organism curves for a composite genome come from gc_bias_by_contig_group, which streams a
// contig subset of the same BAM against that organism's own reference.
process gc_bias {
    tag { library }
    label 'medium_cpu'
    conda "bioconda::picard=3.3.0"
    publishDir "${params.outputDir}/stats/gc_bias"

    input:
        tuple val(library), path(bam), path(bai), val(picard_reference)

    output:
        tuple val(library), path("${library}.gc_metrics"), emit: for_agg
        tuple val("${task.process}"), val('picard'), eval('picard CollectGcBiasMetrics --version 2>&1 | cut -f 2 -d ":"'), topic: versions

    script:
    // The picard flags here must stay in step with modules/gc_bias_by_contig_group.nf.
    """
    picard -Xmx${task.memory.toGiga()}g CollectGcBiasMetrics \\
        --IS_BISULFITE_SEQUENCED true --VALIDATION_STRINGENCY SILENT \\
        -I ${bam} -O ${library}.gc_metrics -S ${library}.gc_summary_metrics \\
        --CHART /dev/null -R ${picard_reference}
    """
}
