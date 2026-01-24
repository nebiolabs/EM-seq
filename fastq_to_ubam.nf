#!/usr/bin/env nextflow
nextflow.enable.dsl=2

input_glob = params.input_glob ?: ['*.{1,2}.fastq.gz']
read_format = params.read_format ?: 'paired-end'
params.outdir = './ubam'

/**
 * Extracts and validates barcode from FASTQ file header.
 *
 * This function reads the first FASTQ record header and attempts to extract a valid barcode,
 * handling both standard EM-seq format and non-standard external dataset formats (e.g., SRR).
 *
 * Extraction strategy:
 * 1. Try to extract field 10 (colon-delimited) for standard EM-seq format
 * 2. If field 10 is empty or contains spaces, fall back to the last colon-delimited field
 * 3. Remove any trailing text after spaces from the candidate
 *
 * Validation rules:
 * - Barcode must match the pattern: ^[ACGTN+-]+$
 * - Only nucleotide bases (A, C, G, T, N) and barcode separators (+, -) are allowed
 * - No spaces or other characters are permitted
 *
 * @param fastqFile The FASTQ file path to extract barcode from
 * @return Shell script that sets the 'barcode' variable to either:
 *         - The extracted valid barcode (e.g., "GCTTCACAAT+TAGCTTTAAC")
 *         - "UNKNOWN" if no valid barcode is found or validation fails
 *
 * Examples:
 * - Standard EM-seq: "@AV100001:...:0063 1:N:0:GCTTCACAAT+TAGCTTTAAC" → "GCTTCACAAT+TAGCTTTAAC"
 * - Non-standard SRR: "@SRR20318439.1 A00536:248:HFHTKDSX3:1:1101:2736:1000 length=111" → "UNKNOWN"
 */
def extractBarcode(fastqFile) {
    """
    set +o pipefail

    barcode=\$(zcat ${fastqFile} | head -n 1 | awk -F: '{
        candidate = \$10
        if (candidate == "" || candidate ~ / /) {
            candidate = \$NF
        }
        gsub(/ .*/, "", candidate)
        if (candidate ~ /^[ACGTN+-]+\$/) {
            print candidate
        } else {
            print "UNKNOWN"
        }
    }')

    set -o pipefail
    """
}

process FastqToBamPaired {
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.21"
    publishDir "${params.outdir}", mode: 'copy'
    memory { params.max_memory ?: 300.GB }
    
    input:
        tuple val(library), path(read1), path(read2)
    
    output:
        path('*.bam')
        path('*.barcode.txt')

    script:
    """
    ${extractBarcode(read1)}
    
    echo "\$barcode" > ${library}.barcode.txt

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
        path('*.barcode.txt')

    script:
    """
    ${extractBarcode(read1)}
    
    echo "\$barcode" > ${library}.barcode.txt

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
