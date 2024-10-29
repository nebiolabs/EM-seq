
process gc_bias {
    cpus 1
    tag { library }
    conda "bioconda::picard=2.20.7 bioconda::samtools=1.9"
    publishDir "${library}/stats/gc_bias"

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
    publishDir "${library}/stats/idxstats"

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
    publishDir "${library}/stats/flagstats"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(params.email), val(library), path("*flagstat"), emit: for_agg

    shell:
    '''
    samtools flagstat -@!{task.cpus} !{bam} > !{library}.flagstat
    '''
}

process fast_qc {
    cpus 1
    tag { library }
    conda "bioconda::fastqc=0.11.8"
    publishDir "${library}/stats/fastqc"

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
    publishDir "${library}/stats/insert_size"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(params.email), val(library), path('*_metrics'), emit: for_agg

    shell:
    '''
    # Above params.min_mapq value (default=20)
    echo "insert_size   2bparsedByRailsModel" > !{library}_insertsize_metrics

    mkfifo good_mapq
    mkfifo bad_mapq
    samtools view -hq20 -U bad_mapq !{bam} > good_mapq &

    picard -Xmx16g CollectInsertSizeMetrics --INCLUDE_DUPLICATES --VALIDATION_STRINGENCY SILENT -I good_mapq -O temp.out.txt --MINIMUM_PCT 0 --Histogram_FILE /dev/null &
    picard -Xmx16g CollectInsertSizeMetrics --INCLUDE_DUPLICATES --VALIDATION_STRINGENCY SILENT -I bad_mapq -O temp2.out.txt --MINIMUM_PCT 0 --Histogram_FILE /dev/null &

    while [[ ! -f "temp.out.txt" || ! -f "temp2.out.txt" ]]; do sleep 2; done

    rm -f good_mapq bad_mapq

    cat temp.out.txt | awk -v mapq=">=!{params.min_mapq}" '
        BEGIN{n="_"; cols=""}
        {
            if (length($0)==0) {} 
            else if (n!="_"){
                for (i=0; i<4-NF;i++) {cols=cols"\t0"};
                print $0""cols"\\t"mapq; cols=""} 
            else {print $0}; 
            if ($1~/^insert_size/) {n=mapq;} 
        }' > !{library}_insertsize_metrics
        
    # Below params.min_mapq value 
    cat temp2.out.txt | awk -v mapq="<!{params.min_mapq}" 'BEGIN{n="_";cols=""}{ if (length($0)==0) {} else if (n!="_"){
        for (i=0; i<4-NF;i++) {cols=cols"\t0"}; print $0""cols"\\t"mapq; cols=""}; if ($1~/^insert_size/) {n=mapq;}}' \
    >> !{library}_insertsize_metrics

    # We could have 2 (fr) or 3 columns (fr AND rf) columns. BTW, we need 4 columns for the model to accept this data. We fill the rest with zeros 
    '''
}

process picard_metrics {
    cpus 1
    tag { library }
    conda "bioconda::picard=2.20.7"
    publishDir "${library}/stats/picard_alignment_metrics"

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
    label 'cpus_8'
    tag { library }
    publishDir "${library}/stats/tasmanian"
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
