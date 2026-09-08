process mergeAndPicodup {
    label 'dedup'
    tag { library }
    publishDir "${params.outputDir}/markduped_bams", mode: 'copy', pattern: '*.md.{bam,bai}'
    publishDir "${params.outputDir}/stats/markdups", mode: 'copy', pattern: '*log'
    conda "bioconda::samtools=1.22"
    memory { params.max_memory ?: 16.GB }

    input:
        tuple val(library), path(bams)
    output:
        tuple val(library), path("${library}.md.bam"), path("${library}.md.bai"), emit: md_bams
        tuple val(library), path('*.markdups_log'), emit: log
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions
        tuple val("${task.process}"), val('picodup'), eval("${params.picodup_bin} --version | cut -f 2 -d ' '"), topic: versions

    script:
    def merge_cpus = Math.min(2, task.cpus)
    def dedup_cpus = Math.max(1, task.cpus - merge_cpus)
    def optical_distance_line = params.picodup_optical_distance != null
        ? "optical_distance=${params.picodup_optical_distance}"
        : 'optical_distance=$(echo ${inst_name} | awk \'{if ($1~/^M0|^NS|^NB/) {print 100} else {print 2500}}\')'
    """
    set +o pipefail
    inst_name=\$(samtools view ${bams[0]} | head -n1 | cut -d ":" -f1)
    set -o pipefail

    ${optical_distance_line}

    samtools merge --threads ${merge_cpus} -O bam -o - ${bams} \\
    | ${params.picodup_bin} \\
        --input /dev/stdin \\
        --output ${library}.md.bam \\
        --metrics ${library}.markdups_log \\
        --threads ${dedup_cpus} \\
        --optical-model ${params.picodup_optical_model} \\
        --optical-distance \${optical_distance} \\
        --tagging-policy All

    samtools index -@ ${task.cpus} ${library}.md.bam ${library}.md.bai
    """
}
