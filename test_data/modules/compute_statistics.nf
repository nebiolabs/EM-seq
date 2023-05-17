
process gc_bias {
    cpus 1
    tag { library }
    conda "picard"
    publishDir "${library}/stats/gc_bias"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(library), path('*gc_metrics'), emit: for_agg

    shell:
    '''
    picard -Xmx4g CollectGcBiasMetrics IS_BISULFITE_SEQUENCED=true VALIDATION_STRINGENCY=SILENT I=!{bam} O=!{library}.gc_metrics S=!{library}.gc_summary_metrics CHART=!{library}.gc.pdf R=!{params.genome}
    '''
}

process idx_stats {
    label 'cpus_8'
    tag { flowcell }
    conda "samtools"
    publishDir "${library}/stats/idxstats"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(library), path(idxstats), emit: for_agg

    shell
    '''
    samtools idxstats !{bam} > !{library}.idxstat
    '''
}

process flag_stats {
    label 'cpus_8'
    tag { flowcell }
    conda "samtools"
    publishDir "${library}/stats/idxstats"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(library), path(flagstat), emit: for_agg

    shell
    '''
    samtools flagstat -@!{task.cpus} !{bam} > !{library}.flagstat
    '''
}

process fast_qc {
    cpus 1
    tag { flowcell }
    conda "picard"
    publishDir "${library}/stats/fastqc"

    input:
        tuple val(library), path(bam), path(bai)

    output:
        tuple val(library), path('*_fastqc.zip'), emit: for_agg 

    shell:
    '''
    fastqc -f bam !{bam}
    '''
}

process insert_size_metrics {
    cpus 1
    tag { flowcell }
    conda "picard"
    publishDir "${library}/stats/insert_size"

    input:
        tuple val(library), path(bam), path(bai)

    output:
        tuple val(library), path('*_metrics'), emit: for_agg

    shell:
    '''
    picard -Xmx16g CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT  I=!{bam} O=${library}.insertsize_metrics MINIMUM_PCT=0 HISTOGRAM_FILE=/dev/null
    '''
}

process picard_metrics {
    cpus 1
    tag { flowcell }
    conda "picard"
    publishDir "${library}/stats/picard_alignment_metrics"

    input:
        tuple val(library), path(bam), path(bai)

    output:
        tuple val(library), path('*alignment_summary_metrics.txt'), emit: for_agg

    shell:
    '''
    picard -Xmx16g CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=SILENT BS=true R=!{params.genome} I=!{bam} O=!{library}.alignment_summary_metrics.txt
    '''
}

process tasmanian {
    label 'cpus_8'
    tag { flowcell }
    publishDir "${library}/stats/tasmanian"
    conda "samtools tasmanian-mismatch"

    input:
        tuple val(library), path(bam), path(bai)

    output:
        tuple val(library), path('*.csv'), emit: for_agg

    shell:
    '''
    set +e
    set +o pipefail

    samtools view -q 30 -F 3840 !{bam} | head -n 2000000 | run_tasmanian -r !{params.genome} > !{library}.csv
    '''

}