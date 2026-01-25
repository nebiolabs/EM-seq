#!/usr/bin/env nextflow
nextflow.enable.dsl=2

input_glob = params.input_glob ?: ['*.{1,2}.fastq.gz']
read_format = params.read_format ?: 'paired-end'
params.outdir = './ubam'

process FastqToBamPaired {
    conda "bioconda::samtools=1.23"
    publishDir "${params.outdir}", mode: 'copy'
    memory { params.max_memory ?: 300.GB }
    
    input:
        tuple val(library), path(read1), path(read2)
    
    output:
        path('*.bam')

    script:
    """
    # Extract most frequent barcode from first 10k reads
    # Split on : and take last field, then remove everything after any space
    barcode=\$(zcat ${read1} | sed -n '1~4p' | head -n 10000 | awk -F: '{print \$NF}' | awk '{print \$1}' | sort | uniq -c | sort -rn | head -n 1 | awk '{print \$2}')
    
    samtools import -i \
        -r ID:${library} \
        -r SM:${library} \
        -r LB:${library} \
        -r PL:ILLUMINA \
        -r CN:"New England Biolabs" \
        -r BC:\${barcode} \
        -1 ${read1} \
        -2 ${read2} \
        -o ${library}.bam
    """
}

process FastqToBamSingle {
    conda "bioconda::samtools=1.23"
    publishDir "${params.outdir}", mode: 'copy'
    memory { params.max_memory ?: 300.GB }
    
    input:
        tuple val(library), path(read1)
    
    output:
        path('*.bam')

    script:
    """
    # Extract most frequent barcode from first 10k reads
    # Split on : and take last field, then remove everything after any space
    barcode=\$(zcat ${read1} | sed -n '1~4p' | head -n 10000 | awk -F: '{print \$NF}' | awk '{print \$1}' | sort | uniq -c | sort -rn | head -n 1 | awk '{print \$2}')
    
    samtools import -i \
        -r ID:${library} \
        -r SM:${library} \
        -r LB:${library} \
        -r PL:ILLUMINA \
        -r CN:"New England Biolabs" \
        -r BC:\${barcode} \
        ${read1} \
        -o ${library}.bam
    """
}

workflow {

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
