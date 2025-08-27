process flag_stats {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::samtools=1.22"
    publishDir "${params.outputDir}/stats/flagstats"

    input:
        tuple val(library), path(bam), path(bai)

    output:
        tuple val(library), path("*flagstat"), emit: for_agg

    script:
    """
    samtools flagstat -@${task.cpus} ${bam} > ${library}.flagstat
    """
}
