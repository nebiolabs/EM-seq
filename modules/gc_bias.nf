
process gc_bias {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.22"
    publishDir "${params.outputDir}/stats/gc_bias", mode: 'copy'

    input:
        tuple val(library), path(bam), path(bai)
        val(genome_fa)
        val(genome_fai)
    output:
        tuple val(library), path("${library}.gc_metrics"), emit: for_agg
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions
        tuple val("${task.process}"), val('picard'), eval('picard CollectGcBiasMetrics --version 2>&1 | cut -f 2 -d ":"'), topic: versions

    script:
    """
    samtools view -H ${bam} | grep "^@SQ" \
    | grep -v "plasmid_puc19\\|phage_lambda\\|phage_Xp12\\|phage_T4\\|EBV\\|chrM" \
    | awk -F":|\\t" '{print \$3"\\t"0"\\t"\$5}' > include_regions.bed

    samtools view -h -L include_regions.bed ${bam} | \
    picard -Xmx${task.memory.toGiga()}g CollectGcBiasMetrics \
        --IS_BISULFITE_SEQUENCED true --VALIDATION_STRINGENCY SILENT \
        -I /dev/stdin -O ${library}.gc_metrics -S ${library}.gc_summary_metrics \
        --CHART /dev/null -R ${genome_fa}
    """
}
