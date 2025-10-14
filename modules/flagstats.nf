process flagstats {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::samtools=1.22"
    publishDir "${params.outputDir}/stats/flagstats", mode: 'copy'

    input:
        tuple val(library), path(bam), path(bai)

    output:
        tuple val(library), path("${library}.flagstat"), emit: for_agg
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions

    script:
    """
    samtools flagstat -@${task.cpus} ${bam} > ${library}.flagstat
    """
}
