process tasmanian {
    label 'medium_cpu'
    tag { library }
    publishDir "${params.outputDir}/stats/tasmanian"
    conda "bioconda::samtools=1.22 bioconda::tasmanian-mismatch=1.0.9"

    errorStrategy { retry < 1 ? 'retry' : 'terminate' }
    maxRetries 1
    memory { retry > 0 ? '16 GB' : '8 GB' }

    input:
        tuple val(library), path(bam), path(bai)
        val(genome_fa)
        val(genome_fai)

    output:
        tuple val(library), path('*.csv'), emit: for_agg

    script:
    """
    set +e
    set +o pipefail
    samtools view -q 30 -F 3840 ${bam} | head -n 2000000 | run_tasmanian -r ${genome_fa} > ${library}.csv
    """

}
