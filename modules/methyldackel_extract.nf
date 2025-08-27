process methylDackel_extract {
    label 'high_cpu'
    tag "${library}"
    publishDir "${params.outputDir}/methylDackelExtracts", mode: 'copy'
    conda "bioconda::methyldackel=0.6.1 bioconda::samtools=1.21 conda-forge::pigz=2.8"

    input:
        tuple val(library), path(md_bam), path(md_bai)
        val(genome_fa)
        val(genome_fai)

    output:
        tuple val(library), path('*CHG.methylKit.gz'), path('*CHH.methylKit.gz'),path('*CpG.methylKit.gz')

    script:
    """
    MethylDackel extract --methylKit -q 20 --nOT 0,0,0,5 --nOB 0,0,5,0 -@ ${task.cpus} \
        --CHH --CHG -o ${library} ${genome_fa} "${md_bam}"
    pigz -p ${task.cpus} *.methylKit
    """
}
