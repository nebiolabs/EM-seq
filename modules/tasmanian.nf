process tasmanian {
    label 'medium_cpu'
    tag { library }
    publishDir "${params.outputDir}/stats/tasmanian"
    conda "bioconda::tasmanian-mismatch=2.0.5"

    input:
        tuple val(library), path(bam), path(bai)
        val(genome_fa)
        val(genome_fai)

    output:
        tuple val(library), path("${library}.tasmanian.csv"), emit: for_agg
        tuple val("${task.process}"), val('tasmanian-mismatch'), eval('tasmanian-mismatch --version | cut -f 2 -d " "'), topic: versions

    script:
    """
    tasmanian-mismatch --min-base-quality 20 --min-map-quality 30 \\
        -F 3840 --threads ${task.cpus} -o ${library}.tasmanian.csv ${bam} ${genome_fa}
    """
}
