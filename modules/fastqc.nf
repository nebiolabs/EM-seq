process fastqc {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::fastqc=0.11.8"
    publishDir "${params.outputDir}/stats/fastqc"

    input:
        tuple val(library), path(bam), path(bai)

    output:
        tuple val(library), path('*_fastqc.zip'), emit: for_agg
        tuple val(library), path('*_fastqc.html'), emit: html

    shell:
    """
    fastqc -f bam ${bam}
    """
}
