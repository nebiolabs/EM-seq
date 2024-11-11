
process gc_bias {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.21"
    publishDir "${params.outputDir}/stats/gc_bias"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(params.email), val(library), path('*gc_metrics'), emit: for_agg

    shell:
    '''
    samtools view -H !{bam} | grep "^@SQ" \
    | grep -v "plasmid_puc19\\|phage_lambda\\|phage_Xp12\\|phage_T4\\|EBV\\|chrM" \
    | awk -F":|\\t" '{print $3"\\t"0"\\t"$5}' > include_regions.bed

    samtools view -h -L include_regions.bed !{bam} | \
    picard -Xmx!{task.memory.toGiga()}g CollectGcBiasMetrics \
        --IS_BISULFITE_SEQUENCED true --VALIDATION_STRINGENCY SILENT \
        -I /dev/stdin -O !{library}.gc_metrics -S !{library}.gc_summary_metrics \
        --CHART !{library}.gc.pdf -R !{params.genome}
    '''
}

process idx_stats {
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
    label 'medium_cpu'
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
    label 'medium_cpu'
    tag { library }
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.21"
    publishDir "${params.outputDir}/stats/insert_size"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(params.email), val(library), path('*_metrics'), emit: for_agg
        tuple val(params.email), val(library), path('*good_mapq.insert_size_metrics.txt'), emit: high_mapq_insert_size_metrics

    shell:
    '''

    good_mapq=$(mktemp -u good_mapq.XXXXXX)
    bad_mapq=$(mktemp -u bad_mapq.XXXXXX)
    mkfifo "$good_mapq"
    mkfifo "$bad_mapq"
    trap "rm -f $good_mapq $bad_mapq" EXIT # cleanup upon exit

    samtools view -h -q 20 -U "$bad_mapq" !{bam} > "$good_mapq" &
    samtools_pid=$!
    picard -Xmx!{task.memory.toGiga()}g CollectInsertSizeMetrics \

        --INCLUDE_DUPLICATES --VALIDATION_STRINGENCY SILENT -I "$good_mapq" -O good_mapq.out.txt \
        --MINIMUM_PCT 0 -H /dev/null &
    picard_good_mapq_pid=$!

    picard -Xmx!{task.memory.toGiga()}g CollectInsertSizeMetrics \
        --INCLUDE_DUPLICATES --VALIDATION_STRINGENCY SILENT -I "$bad_mapq" -O bad_mapq.out.txt \
        --MINIMUM_PCT 0 -H /dev/null &
    picard_bad_mapq_pid=$!

    # Wait for programs to finish, named pipes will be closed by the trap
    wait $samtools_pid
    wait $picard_good_mapq_pid
    wait $picard_bad_mapq_pid


    # extract the leading lines from the "good" mapq file
    grep -B 1000 '^insert_size' good_mapq.out.txt > !{library}_insertsize_metrics
    # add a column for the mapq category to the hisogram header
    sed -i 's/count$/count\\tcategory/' !{library}_insertsize_metrics
    grep '^[0-9][^A-Z]*$' good_mapq.out.txt | awk '{print $0"\\t>=20"}' >> !{library}_insertsize_metrics
    grep '^[0-9][^A-Z]*$' bad_mapq.out.txt | awk '{print $0"\\t<20"}' >> !{library}_insertsize_metrics

    # make sure we have 5 columns not 4 (if tandem_counts = 0 )
    cat !{library}_insertsize_metrics | awk '{
        if (header_seen==1) {
            if (add_tandem==1) {
                print $1"\t"$2"\t"$3"\t"0"\t"$4
            }
            else {print $0}
        } 
        else {print $0}
        if ($1~/^insert_size/) {
            header_seen=1; 
            if (NF==4) {add_tandem=1;}
        }
    }' > tmp && mv tmp !{library}_insertsize_metrics
    '''
}

process picard_metrics {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.21"
    publishDir "${params.outputDir}/stats/picard_alignment_metrics"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(params.email), val(library), path('*alignment_summary_metrics.txt'), emit: for_agg

    shell:
    '''
    picard -Xmx!{task.memory.toGiga()}g CollectAlignmentSummaryMetrics \
        --VALIDATION_STRINGENCY SILENT -BS true -R !{params.genome} \
        -I !{bam} -O !{library}.alignment_summary_metrics.txt
    '''
}

process tasmanian {
    label 'medium_cpu'
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
