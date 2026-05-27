#!/usr/bin/env nextflow
nextflow.enable.dsl=2

input_glob = params.input_glob ?: ['*.{1,2}.fastq.gz']
read_format = params.read_format ?: 'paired-end'
params.outdir = './ubam'

// Shared shell function to extract and validate barcode from FASTQ header
// Samples first 10k reads and returns the most frequent valid barcode
def extractBarcodeFunction = '''
extract_barcode() {
    local fastq_file="$1"

    # Extract last colon-field from comment (after space) of first 10k read headers
    # Filter to valid barcodes (nucleotides with optional +), count occurrences, return most frequent
    barcode=$(zcat "$fastq_file" \
        | head -n 40000 \
        | awk 'NR % 4 == 1 {sub(/.*[[:space:]]/, ""); n=split($0,a,":"); print a[n]}' \
        | grep -E '^[ACGTN+-]+$' \
        | sort | uniq -c | sort -rn | head -1 | awk '{print $2}')

    # Fallback to unknown if no valid barcode found
    echo "${barcode:-unknown}"
}
'''

process FastqToBamPaired {
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.21"
    publishDir "${params.outdir}", mode: 'copy'
    memory { params.max_memory ?: 300.GB }

    input:
        tuple val(library), path(read1), path(read2)

    output:
        path('*.bam')

    script:
    """
    set +o pipefail
    ${extractBarcodeFunction}
    barcode=\$(extract_barcode "${read1}")
    set -o pipefail

    picard FastqToSam TMP_DIR=/state/partition1/sge_tmp F1=${read1} F2=${read2} OUTPUT=temp.bam SM=${library} LB=${library} CN="New England Biolabs" PU=Illumina QUIET=true

    samtools reheader -c "sed \\"s/RG/RG\\tBC:\$barcode/\\"" temp.bam > ${library}.bam
    rm temp.bam
    """
}

process FastqToBamSingle {
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.21"
    publishDir "${params.outdir}", mode: 'copy'
    memory { params.max_memory ?: 300.GB }

    input:
        tuple val(library), path(read1)

    output:
        path('*.bam')

    script:
    """
    set +o pipefail
    ${extractBarcodeFunction}
    barcode=\$(extract_barcode "${read1}")
    set -o pipefail

    picard FastqToSam F1=${read1} OUTPUT=temp.bam SM=${library} LB=${library} CN="New England Biolabs" PU=Illumina QUIET=true

    samtools reheader -c "sed \\"s/RG/RG\\tBC:\$barcode/\\"" temp.bam > ${library}.bam
    rm temp.bam
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
