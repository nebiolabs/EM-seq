#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_glob = ['*.{1,2}.fastq.gz', '*_R{1,2}_*.fastq.gz']
params.read_format = 'paired-end'
params.max_memory = 4.GB
params.outdir = './ubam'
params.platform = "Illumina"
params.platform_model = "?"
params.center_name = null

input_glob = params.input_glob
read_format = params.read_format

process FastqToBamPaired {
    conda "bioconda::samtools=1.22"
    publishDir "${params.outdir}", mode: 'copy'
    memory { params.max_memory ?: 300.GB }

    input:
        tuple val(library), path(read1), path(read2)

    output:
        path('*.bam')

    script:
    def cn_flag = params.center_name ? "-r \"CN:${params.center_name}\"" : ""
    """
    set +o pipefail

    barcode=\$(zcat ${read1} | head -n 40000 | awk 'NR % 4 == 1' | cut -d ":" -f 10 | sort | uniq -c | sort -rn | head -n 1 | awk '{print \$2}')

    set -o pipefail

    samtools import -i -1 ${read1} -2 ${read2} \\
        -r "PL:${params.platform}" \\
        -r "PM:${params.platform_model}" \\
        -r "ID:1" \\
        -r "LB:${library}" \\
        -r "SM:${library}" \\
        -r "BC:\$barcode" \\
        ${cn_flag} \\
        -o ${library}.bam
    """
}

process FastqToBamSingle {
    conda "bioconda::samtools=1.22"
    publishDir "${params.outdir}", mode: 'copy'
    memory { params.max_memory ?: 300.GB }

    input:
        tuple val(library), path(read1)

    output:
        path('*.bam')

    script:
    def cn_flag = params.center_name ? "-r \"CN:${params.center_name}\"" : ""
    """
    set +o pipefail

    barcode=\$(zcat ${read1} | head -n 40000 | awk 'NR % 4 == 1' | cut -d ":" -f 10 | sort | uniq -c | sort -rn | head -n 1 | awk '{print \$2}')

    set -o pipefail

    samtools import -i -1 ${read1} \\
        -r "PL:${params.platform}" \\
        -r "PM:${params.platform_model}" \\
        -r "ID:1" \\
        -r "LB:${library}" \\
        -r "SM:${library}" \\
        -r "BC:\$barcode" \\
        ${cn_flag} \\
        -o ${library}.bam
    """
}

workflow {

    // Display read group information
    log.info """
    Read Group Settings:
      Platform (PL):       ${params.platform}
      Platform Model (PM): ${params.platform_model}
      Center Name (CN):    ${params.center_name ?: 'not set'}
      Barcode (BC):        extracted from FASTQ header
      Sample ID (SM):      derived from FASTQ filename
      Library ID (LB):     derived from FASTQ filename

    To customize, use: --platform "..." --platform_model "..." --center_name "..."
    """.stripIndent()

    if (read_format == 'paired-end') {
        fastq_files = Channel.fromFilePairs(input_glob, flat: true)
        FastqToBamPaired(fastq_files)
    }
    else if (read_format == 'single-end') {
        fastq_files = Channel.fromPath(input_glob).map{it-> [it.baseName.split('.fastq')[0], it]}
        FastqToBamSingle(fastq_files)
    }
    else {
        error "Unknown read format -- accepted is 'paired-end' or 'single-end'"
    }
}
