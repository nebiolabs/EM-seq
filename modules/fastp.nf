process fastp {
	tag "${library}"
    label 'low_cpu'
    conda "bioconda::samtools=1.21 bioconda::fastp=1.0.1"
    publishDir "${params.outputDir}/fastp"

    input:
        tuple val(library), path(bam)
    
    output:
        tuple val(library), path("*.trimmed.fastq"), emit: trimmed_fastq
        tuple val(library), path("${library}.fastp.json"), emit: fastp_json
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions
        tuple val("${task.process}"), val('fastp'), eval('fastp --version 2>&1 | cut -f 2 -d " "'), topic: versions
    
    script:
    def fastp_args = params.single_end ? "--out1 ${library}.1.trimmed.fastq" : "--interleaved_in --out1 ${library}.1.trimmed.fastq --out2 ${library}.2.trimmed.fastq"
    """
    set +o pipefail
    inst_name=\$(samtools view ${bam} | head -n 1 | cut -d ":" -f 1)
    set -o pipefail

    trim_polyg=\$(echo "\${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')
    echo \${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'
    
    samtools fastq -n ${bam} | \\
    fastp --stdin \\
                    -l 2 -Q \${trim_polyg} \\
                    --thread 1 \\
                    --overrepresentation_analysis \\
                    -j "${library}.fastp.json" \\
                    --split_by_lines ${params.fastq_split_lines} \\
                    ${fastp_args}
    """
}
