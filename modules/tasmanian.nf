process tasmanian {
    label 'medium_cpu'
    tag { library }
    publishDir "${params.outputDir}/stats/tasmanian"
    conda "bioconda::tasmanian-mismatch=2.1.0"

    errorStrategy { task.attempt <= 1 ? 'retry' : 'terminate' }
    maxRetries 1
    memory { task.attempt > 1 ? '16 GB' : '8 GB' }

    input:
        tuple val(library), path(bam), path(bai)
        val(genome_fa)
        val(genome_fai)

    output:
        tuple val(library), path("${library}.tasmanian.tsv"), emit: for_agg
        tuple val("${task.process}"), val('tasmanian-mismatch'), val('*should be* 2.1.0'), topic: versions

    script:
    """
    set -euo pipefail
    for mode in read insert; do
        tasmanian-mismatch ${bam} ${genome_fa} \
            --threads ${task.cpus} \
            --position-mode \${mode} \
            --methylation-mode \
            --min-base-quality 20 \
            --min-map-quality 30 \
            -F 3840 \
            -o ${library}.tasmanian.\${mode}.tsv
    done

    head -n1 ${library}.tasmanian.read.tsv | awk '{print \$0"\tpair_mode"}' > ${library}.tasmanian.tsv
    tail -n +2 ${library}.tasmanian.read.tsv | awk '{print \$0"\tread"}' >> ${library}.tasmanian.tsv
    tail -n +2 ${library}.tasmanian.insert.tsv | awk '{print \$0"\tinsert"}' >> ${library}.tasmanian.tsv
    """

}
