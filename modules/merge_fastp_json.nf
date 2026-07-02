process mergeFastpJson {
	tag "${library}"
    label 'low_cpu'
    conda 'bioconda::bamslice=0.2.1'
    publishDir "${params.outputDir}/fastp"

    input:
        tuple val(library), val(_slices), path(json_files)

    output:
        tuple val(library), path("${library}.fastp.json"), emit: merged_json
        tuple val("${task.process}"), val('bamslice'), eval('bamslice --version | cut -f 2 -d " "'), topic: versions

    script:
    """
    fastp-merge ${json_files.join(' ')} -o ${library}.fastp.json
    """
}
