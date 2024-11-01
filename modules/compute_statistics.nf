
process gc_bias {
    cpus 1
    tag { library }
    conda "bioconda::picard=2.20.7 bioconda::samtools=1.9"
    publishDir "${params.outputDir}/stats/gc_bias"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(params.email), val(library), path('*gc_metrics'), emit: for_agg

    shell:
    '''
    samtools view -H !{bam} | grep "^@SQ" | grep -v "plasmid_puc19\\|phage_lambda\\|phage_Xp12\\|phage_T4\\|EBV\\|chrM" | awk -F":|\\t" '{print $3"\\t"0"\\t"$5}' > include_regions.bed
    samtools view -h -L include_regions.bed !{bam} | \
    picard -Xmx4g CollectGcBiasMetrics IS_BISULFITE_SEQUENCED=true VALIDATION_STRINGENCY=SILENT I=/dev/stdin O=!{library}.gc_metrics S=!{library}.gc_summary_metrics CHART=!{library}.gc.pdf R=!{params.genome}
    '''
}

process idx_stats {
    label 'cpus_8'
    tag { library }
    conda "bioconda::samtools=1.9"
    publishDir "${params.outputDir}/stats/idxstats"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(params.email), val(library), path("*idxstat"), emit: for_agg

    shell:
    '''
    samtools idxstats !{bam} > !{library}.idxstat
    '''
}

process flag_stats {
    label 'cpus_8'
    tag { library }
    conda "bioconda::samtools=1.9"
    publishDir "${params.outputDir}/stats/flagstats"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(params.email), val(library), path("*flagstat"), emit: for_agg

    shell:
    '''
    samtools flagstat -@!{task.cpus} !{bam} > !{library}.flagstat
    '''
}

process fastqc {
    cpus 1
    tag { library }
    conda "bioconda::fastqc=0.11.8"
    publishDir "${params.outputDir}/stats/fastqc"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(params.email), val(library), path('*_fastqc.zip'), emit: for_agg

    shell:
    '''
    fastqc -f bam !{bam}
    '''
}

process insert_size_metrics {
    cpus 1
    tag { library }
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.21"
    publishDir "${params.outputDir}/stats/insert_size"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(params.email), val(library), path('*_metrics'), emit: for_agg

    shell:
    '''

    good_mapq=$(mktemp -u /tmp/good_mapq.XXXXXX) # -u is unique
    bad_mapq=$(mktemp -u /tmp/bad_mapq.XXXXXX)
    mkfifo "$good_mapq"
    mkfifo "$bad_mapq"
    trap "rm -f $good_mapq $bad_mapq" EXIT # cleanup upon exit

    samtools view -hq20 -U "$bad_mapq" !{bam} > "$good_mapq" &
    samtools_pid=$!

    picard -Xmx16g CollectInsertSizeMetrics --INCLUDE_DUPLICATES --VALIDATION_STRINGENCY SILENT -I "$good_mapq" -O good_mapq.out.txt --MINIMUM_PCT 0 --Histogram_FILE /dev/null &
    picard -Xmx16g CollectInsertSizeMetrics --INCLUDE_DUPLICATES --VALIDATION_STRINGENCY SILENT -I "$bad_mapq" -O bad_mapq.out.txt --MINIMUM_PCT 0 --Histogram_FILE /dev/null &

    # Wait for samtools to finish, then close named pipes
    wait $samtools_pid

    # do we need these pipes to stay open ?
    #exec 3<>"$good_mapq"
    #exec 3<>"$bad_mapq"

    # Wait for remaining processes
    wait

    # how about this?
    #extract the comments from the "good" mapq file
    grep '^#' good_mapq.out.txt > !{library}_insertsize_metrics
    #get the header line from the "good" mapq file and add a column for the mapq category
    grep -i "rf.count" good_mapq.out.txt  | sed 's/$/\tcategory/' >> testing_insertsize_metrics
    #add a column for the mapq category to the data 
    grep '^[0-9][^A-Z]*$' good_mapq.out.txt | awk '{print $0"\\t>=20"}' >> !{library}_insertsize_metrics
    grep '^[0-9][^A-Z]*$' bad_mapq.out.txt | awk '{print $0"\\t<20"}' >> !{library}_insertsize_metrics
    '''
}

process picard_metrics {
    cpus 1
    tag { library }
    conda "bioconda::picard=2.20.7"
    publishDir "${params.outputDir}/stats/picard_alignment_metrics"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(params.email), val(library), path('*alignment_summary_metrics.txt'), emit: for_agg

    shell:
    '''
    picard -Xmx16g CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=SILENT BS=true R=!{params.genome} I=!{bam} O=!{library}.alignment_summary_metrics.txt
    '''
}

process tasmanian {
    cpus 8
    tag { library }
    publishDir "${params.outputDir}/stats/tasmanian"
    conda "bioconda::samtools=1.9 bioconda::tasmanian-mismatch=1.0.7"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(params.email), val(library), path('*.csv'), emit: for_agg

    shell:
    '''
    set +e
    set +o pipefail

    samtools view -q 30 -F 3840 !{bam} | head -n 2000000 | run_tasmanian -r !{params.genome} > !{library}.csv
    '''

}
