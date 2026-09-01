process tasmanian {
    label 'medium_cpu'
    tag { library }
    publishDir "${params.outputDir}/stats/tasmanian"
    conda "bioconda::tasmanian-mismatch=2.0.3"

    errorStrategy { retry < 1 ? 'retry' : 'terminate' }
    maxRetries 1
    memory { retry > 0 ? '16 GB' : '8 GB' }

    input:
        tuple val(library), path(bam), path(bai)
        val(genome_fa)
        val(genome_fai)

    output:
        tuple val(library), path("${library}.tasmanian.csv"), emit: for_agg
        tuple val("${task.process}"), val('tasmanian-mismatch'), val('*should be* 2.0.3'), topic: versions

    script:
    """
    set +e
    set +o pipefail
    tasmanian-mismatch ${bam} ${genome_fa} \
        --position-mode read \
        --min-base-quality 20 \
        --min-map-quality 30 \
        -F 3840 \
        -o ${library}.tasmanian.csv
    """

}
