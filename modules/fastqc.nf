process fastqc {
    label 'single_threaded_qc'
    tag { library }
    conda "bioconda::falco=1.0.0 conda-forge::zip=3.0"
    publishDir "${params.outputDir}/stats/fastqc"

    input:
        tuple val(library), path(bam), path(bai)

    output:
        tuple val(library), path('*_fastqc.zip'), emit: for_agg
        tuple val(library), path('*_fastqc.html'), emit: html
        tuple val("${task.process}"), val('falco'), eval('falco --version | cut -f 2 -d " "'), topic: versions

    script:
    // falco's own output is flat (fastqc_data.txt, fastqc_report.html, summary.txt;
    // no per-sample name or archive). Repackage into a <prefix>_fastqc.zip containing
    // a <prefix>_fastqc/ directory — the layout MultiQC's fastqc module and
    // aggregate_results.nf's `unzip *fastqc.zip` + "fastqc/fastqc_data.txt"
    // substitution both expect. <prefix> must be the BAM filename minus ".bam"
    // (bam.baseName), not the library tag, to match what FastQC itself would name it.
    """
    falco -f bam ${bam}
    mkdir ${bam.baseName}_fastqc
    mv fastqc_data.txt summary.txt ${bam.baseName}_fastqc/
    zip -rq ${bam.baseName}_fastqc.zip ${bam.baseName}_fastqc
    mv fastqc_report.html ${bam.baseName}_fastqc.html
    """
}
