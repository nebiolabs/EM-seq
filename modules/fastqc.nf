process fastqc {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::fastqc=0.11.8"
    publishDir "${params.outputDir}/stats/fastqc", mode: 'copy'

    input:
        tuple val(library), path(bam), path(bai)

    output:
        tuple val(library), path('*_fastqc.zip'), emit: for_agg
        tuple val(library), path('*_fastqc.html'), emit: html
        tuple val("${task.process}"), val('fastqc'), eval('fastqc --version | cut -f 2 -d " "'), topic: versions

    shell:
    """
    fastqc -f bam ${bam}
    """
}
