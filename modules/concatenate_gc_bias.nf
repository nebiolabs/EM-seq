// picard writes a multi-line '##' preamble before its header row. Taking line 1 of the first file
// and lines 2+ of every file leaves each file's preamble and header inline, which is harmless:
// GcBias.parse in ngs-aggregate_results skips '#' lines and the repeated column header.
process concatenate_gc_bias {
    tag { library }
    label 'process_single'
    publishDir "${params.outputDir}/stats/gc_bias"

    input:
        tuple val(library), path(files)

    output:
        tuple val(library), path("${library}.gc_metrics.combined.tsv"), emit: for_agg

    script:
    """
    head -n 1 ${files[0]} > ${library}.gc_metrics.combined.tsv
    for file in ${files}; do
        tail -n +2 \$file >> ${library}.gc_metrics.combined.tsv
    done
    """
}
