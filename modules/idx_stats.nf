process idx_stats {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::samtools=1.22"
    publishDir "${params.outputDir}/stats/idxstats"

    input:
        tuple val(library), path(bam), path(bai)

    output:
        tuple val(library), path("*idxstat"), emit: for_agg

    script:
    """
    samtools idxstats ${bam} > ${library}.idxstat
    """
}