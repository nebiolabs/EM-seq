process fastp {
	tag "${library}"
    label 'medium_cpu'
    conda "bioconda::samtools=1.21 bioconda::fastp=1.0.1"
    publishDir "${params.outputDir}/fastp"

    input:
        tuple val(library), path(bam)
    
    output:
        tuple val(library), path("*.1.trimmed.fastq"), path("*.2.trimmed.fastq"), emit: trimmed_fastq
        tuple val(library), path("${library}.fastp.json"), emit: fastp_json

    script:
    """
    set +o pipefail
    inst_name=\$(samtools view ${bam} | head -n 1 | cut -d ":" -f 1)
    set -o pipefail

    trim_polyg=\$(echo "\${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')
    echo \${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'
    
    samtools fastq --threads ${task.cpus} -n ${bam} | \\
    fastp --interleaved_in --stdin \\
                    -l 2 -Q \${trim_polyg} \\
                    --thread ${task.cpus} \\
                    --overrepresentation_analysis \\
                    -j "${library}.fastp.json" \\
                    --split_by_lines ${params.fastq_split_lines} \\
                    --out1 ${library}.1.trimmed.fastq \\
                    --out2 ${library}.2.trimmed.fastq
    """
}
