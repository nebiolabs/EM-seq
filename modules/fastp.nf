process fastp {
	tag "${library}"
    label 'low_cpu'
    conda "bioconda::samtools=1.21 bioconda::fastp=1.3.6 bioconda::bamslice=0.2.1"
    publishDir "${params.outputDir}/fastp"

    input:
        tuple val(library), path(bam), val(start_offset), val(end_offset)

    output:
        tuple val(library), val("${library}_${start_offset}_${end_offset}"), path("${library}_${start_offset}_${end_offset}*.trimmed.fastq.gz"), emit: trimmed_fastq
        tuple val(library), val("${library}_${start_offset}_${end_offset}"), path("${library}_${start_offset}_${end_offset}.fastp.json"), emit: fastp_json
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions
        tuple val("${task.process}"), val('fastp'), eval('fastp --version 2>&1 | cut -f 2 -d " "'), topic: versions
        tuple val("${task.process}"), val('bamslice'), eval('bamslice --version | cut -f 2 -d " "'), topic: versions

    script:
    def chunk = "${library}_${start_offset}_${end_offset}"
    def fastp_args = params.single_end ? "--out1 ${chunk}.1.trimmed.fastq.gz" : "--interleaved_in --out1 ${chunk}.1.trimmed.fastq.gz --out2 ${chunk}.2.trimmed.fastq.gz"
    """
    set +o pipefail
    inst_name=\$(samtools view ${bam} | head -n 1 | cut -d ":" -f 1)
    set -o pipefail

    trim_polyg=\$(echo "\${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')
    echo \${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'

    bamslice --input ${bam} --start-offset ${start_offset} --end-offset ${end_offset} | \\
    fastp --stdin \\
                    -l 2 -Q \${trim_polyg} \\
                    --thread 1 \\
                    --overrepresentation_analysis \\
                    -j "${chunk}.fastp.json" \\
                    ${fastp_args}
    """
}
