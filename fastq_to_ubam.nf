#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_glob = params.input_glob ?: ['*.{1,2}.fastq.gz']
params.read_format = params.read_format ?: 'paired-end'
params.outdir = './ubam'
params.sequencing_center = params.sequencing_center ?: 'Unknown'
params.max_memory = null

process DetectBarcode {
    conda "bioconda::samtools=1.21"
    cache 'lenient'

    input:
        val(library)
        path(fastq)

    output:
        tuple val(library), stdout

    script:
    """
    set +o pipefail
    barcode=\$(zcat ${fastq} | sed -n '1~4p' | head -n 10000 | grep -oP '[GCATN+\\-]+\$' | sort | uniq -c | sort -rn | head -n 1 | awk '{print \$2}')
    set -o pipefail
    echo "\${barcode:-unknown}"
    """
}

process FastqToBamPaired {
    conda "bioconda::fgbio=2.3.0"
    publishDir "${params.outdir}", mode: 'copy'
    memory { params.max_memory ?: 300.GB }

    input:
        tuple val(library), path(read1), path(read2), val(barcode)

    output:
        path("${library}.bam")

    script:
    """
    fgbio FastqToBam --input ${read1} ${read2} --output ${library}.bam --sample ${library} --library ${library} --barcode ${barcode.trim()} --sequencing-center "${params.sequencing_center}"
    """
}

process FastqToBamSingle {
    conda "bioconda::fgbio=2.3.0"
    publishDir "${params.outdir}", mode: 'copy'
    memory { params.max_memory ?: 300.GB }

    input:
        tuple val(library), path(read1), val(barcode)

    output:
        path("${library}.bam")

    script:
    """
    fgbio FastqToBam --input ${read1} --output ${library}.bam --sample ${library} --library ${library} --barcode ${barcode.trim()} --sequencing-center "${params.sequencing_center}"
    """
}

workflow {

    if (params.read_format == 'paired-end') {
        fastq_files = Channel.fromFilePairs(params.input_glob, flat: true)
        barcodes = DetectBarcode(fastq_files.map{it -> it[0]}, fastq_files.map{it -> it[1]})
        FastqToBamPaired(fastq_files.join(barcodes, by: 0))
    }
    else if (params.read_format == 'single-end') {
        fastq_files = Channel.fromPath(params.input_glob).map{it-> [it.baseName.split('.fastq')[0], it]}
        barcodes = DetectBarcode(fastq_files.map{it -> it[0]}, fastq_files.map{it -> it[1]})
        FastqToBamSingle(fastq_files.join(barcodes, by: 0))
    }
    else {
        error "Unknown read format -- accepted is 'paired-end' or 'single-end'"
    }
}
