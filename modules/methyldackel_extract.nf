// MethylDackel methylKit extraction. contexts is 'all' for CpG, CHG and CHH, or 'cpg' for CpG
// only. The depth-dependent analyses of pre-aligned BAMs (depth_dependent.nf) use 'cpg',
// because the CHG and CHH files are large and unused there. Both read the same reads, through
// params.methyl_quality_threshold and params.methyldackel_trim.
process methylDackel_extract {
    label 'medium_cpu'
    tag "${library}"
    publishDir "${params.outputDir}/methylDackelExtracts", mode: 'copy'
    conda "bioconda::methyldackel=0.6.1 bioconda::samtools=1.21 conda-forge::pigz=2.8"

    input:
        tuple val(library), path(md_bam), path(md_bai)
        val(genome_fa)
        val(genome_fai)
        val(contexts)

    output:
        tuple val(library), path("${library}_CpG.methylKit.gz"), emit: cpg
        tuple val(library), path("${library}_CHG.methylKit.gz"), path("${library}_CHH.methylKit.gz"), optional: true, emit: chg_chh
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions
        tuple val("${task.process}"), val('methyldackel'), eval('MethylDackel --version 2>&1 | cut -f 2 -d ":"'), topic: versions

    script:
    if (!(contexts in ['all', 'cpg'])) {
        error "methylDackel_extract: contexts must be 'all' or 'cpg', got '${contexts}'"
    }
    def context_args = contexts == 'all' ? '--CHH --CHG' : ''
    """
    MethylDackel extract --methylKit -q ${params.methyl_quality_threshold} ${params.methyldackel_trim} -@ ${task.cpus} \\
        ${context_args} -o ${library} ${genome_fa} "${md_bam}"
    pigz -p ${task.cpus} *.methylKit
    """
}
