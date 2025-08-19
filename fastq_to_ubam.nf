#!/usr/bin/env nextflow
nextflow.enable.dsl=2

input_glob = params.input_glob ?: ['*.{1,2}.fastq.gz']
read_format = params.read_format ?: 'paired-end'
params.outdir = './ubam'

process FastqToBamPaired {
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.21"
    publishDir "${params.outdir}", mode: 'copy'
    memory { params.max_memory ?: 300.GB }
    
    input:
        tuple val(library), path(read1), path(read2)
    
    output:
        path('*.bam')

    shell:
    '''
    set +o pipefail

    flowcell=$(zcat !{read1} | head -n 1 | cut -d ":" -f 3)
    barcode=$(zcat !{read1} | head -n 1 | cut -d ":" -f 10)

    set -o pipefail

    picard FastqToSam TMP_DIR=/state/partition1/sge_tmp F1=!{read1} F2=!{read2} OUTPUT=/dev/stdout SM=!{library} LB=!{library} CN="New England Biolabs" PU=Illumina QUIET=true COMPRESSION_LEVEL=0 \
    | samtools view -h - \
    | while read -r line; do if [[ "$line" =~ ^@RG ]]; then printf "%s\\t%s\\n" "$line" "BC:${barcode}"; else printf "%s\\n" "$line"; fi; done \
    | samtools view -hb -o !{library}.bam -
    '''
}

process FastqToBamSingle {
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.21"
    publishDir "${params.outdir}", mode: 'copy'
    memory { params.max_memory ?: 300.GB }
    
    input:
        path(read1)
    
    output:
        path('*.bam')

    shell:
    '''
    set +o pipefail

    library=$(basename !{read1} .fastq.gz)
    flowcell=$(zcat !{read1} | head -n 1 | cut -d ":" -f 3)
    barcode=$(zcat !{read1} | head -n 1 | cut -d ":" -f 10)

    picard FastqToSam F1=!{read1} OUTPUT=/dev/stdout SM=${library} LB=${library} CN="New England Biolabs" PU=Illumina QUIET=true COMPRESSION_LEVEL=0 \
    | samtools view -h - \
    | while read -r line; do if [[ "$line" =~ ^@RG ]]; then printf "%s\\t%s\\n" "$line" "BC:${barcode}"; else printf "%s\\n" "$line"; fi; done \
    | samtools view -hb -o ${library}.bam -
    '''
}

workflow {

    if (read_format == 'paired-end') {
        fastq_files = Channel.fromFilePairs(input_glob, flat: true)
        FastqToBamPaired(fastq_files)
    }
    else if (read_format == 'single-end') {
        fastq_files = Channel.fromPath(input_glob)
        FastqToBamSingle(fastq_files)
    }
    else {
        error "Unknown read format -- accepted is 'paired-end' or 'single-end'"
    }
}
