process extract_cytosine_report {
    label 'medium_cpu'
    tag "${library}"
    conda "bioconda::methyldackel=0.6.1 bioconda::samtools=1.21 conda-forge::pigz=2.8"

    input:
        tuple val(library), path(md_bam), path(md_bai)
        val(genome_fa)
        val(genome_fai)

    output:
        tuple val(library), path('*cytosine_report.txt'), emit: report
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions
        tuple val("${task.process}"), val('methyldackel'), eval('MethylDackel --version 2>&1 | cut -f 2 -d ":"'), topic: versions

    script:
    """
    MethylDackel extract --cytosine_report \\
        -q ${params.methyl_quality_threshold} \\
        --nOT 0,0,0,5 --nOB 0,0,5,0 \\
        -@ ${task.cpus} \\
        --CHH --CHG \\
        -o ${library} \\
        ${genome_fa} \\
        "${md_bam}"
    """
}
