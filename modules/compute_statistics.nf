
process gc_bias {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.21"
    publishDir "${params.outputDir}/stats/gc_bias"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)
        path(genome_path)
    output:
        tuple val(params.email), val(library), path('*gc_metrics'), emit: for_agg

    script:
    """
    genome=\$(ls *.bwameth.c2t.bwt | sed 's/.bwameth.c2t.bwt//')
    samtools view -H ${bam} | grep "^@SQ" \
    | grep -v "plasmid_puc19\\|phage_lambda\\|phage_Xp12\\|phage_T4\\|EBV\\|chrM" \
    | awk -F":|\\t" '{print \$3"\\t"0"\\t"\$5}' > include_regions.bed

    samtools view -h -L include_regions.bed ${bam} | \
    picard -Xmx${task.memory.toGiga()}g CollectGcBiasMetrics \
        --IS_BISULFITE_SEQUENCED true --VALIDATION_STRINGENCY SILENT \
        -I /dev/stdin -O ${library}.gc_metrics -S ${library}.gc_summary_metrics \
        --CHART ${library}.gc.pdf -R \${genome}
    """
}

process idx_stats {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::samtools=1.9"
    publishDir "${params.outputDir}/stats/idxstats"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(params.email), val(library), path("*idxstat"), emit: for_agg

    script:
    """
    samtools idxstats ${bam} > ${library}.idxstat
    """
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

    script:
    """
    samtools flagstat -@${task.cpus} ${bam} > ${library}.flagstat
    """
}

process fastqc {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::fastqc=0.11.8"
    publishDir "${params.outputDir}/stats/fastqc"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(params.email), val(library), path('*_fastqc.zip'), emit: for_agg

    shell:
    """
    fastqc -f bam ${bam}
    """
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

    script:
    """

    good_mapq=\$(mktemp -u good_mapq.XXXXXX)
    bad_mapq=\$(mktemp -u bad_mapq.XXXXXX)
    mkfifo "\$good_mapq"
    mkfifo "\$bad_mapq"
    trap "rm -f \$good_mapq \$bad_mapq" EXIT # cleanup upon exit

    samtools view -h -q 20 -U "\$bad_mapq" ${bam} > "\$good_mapq" &
    samtools_pid=\$!
    picard -Xmx${task.memory.toGiga()}g CollectInsertSizeMetrics \
    --INCLUDE_DUPLICATES --VALIDATION_STRINGENCY SILENT -I "\$good_mapq" -O good_mapq.out.txt \
        --MINIMUM_PCT 0 -H /dev/null &
    picard_good_mapq_pid=\$!

    picard -Xmx${task.memory.toGiga()}g CollectInsertSizeMetrics \
        --INCLUDE_DUPLICATES --VALIDATION_STRINGENCY SILENT -I "\$bad_mapq" -O bad_mapq.out.txt \
        --MINIMUM_PCT 0 -H /dev/null &
    picard_bad_mapq_pid=\$!

    # Wait for programs to finish, named pipes will be closed by the trap
    wait \$samtools_pid
    wait \$picard_good_mapq_pid
    wait \$picard_bad_mapq_pid

    # extract the leading lines from the "good" mapq file
    grep -B 1000 '^insert_size' good_mapq.out.txt | grep -v "insert_size" > ${library}_insertsize_metrics
    echo -e "insert_size\tAll_Reads.fr_count\tAll_Reads.rf_count\tAll_Reads.tandem_count\tcategory" >> ${library}_insertsize_metrics

    grep -h -A1000 '^insert_size' good_mapq.out.txt bad_mapq.out.txt | awk 'BEGIN{flag=0} {
        if (! \$2) {if (\$2 != 0) {next}}
        if (\$1~/^insert_size/) {
            if (flag==0) { category = ">=20"}
            else {category = "<20"}
            flag++;

            # set header with indices of arr to keep order in good and bad mapq files.
            for (i=1; i<=5; i++){
                no_cols[i] = 0
                if      (\$i ~ /insert_size/) {isize=i}
                else if (\$i ~ /rf_count/)    {rf=i}
                else if (\$i ~ /fr_count/)    {fr=i}
                else if (\$i ~ /tandem/)      {tandem=i}
                else                         {no_cols[i]++}
            }
            # columns that are not present still need to be printed (with 0 value)
            for (i in no_cols) {
                if (no_cols[i] > 0) {
                    if      (! isize )  { isize = i}
                    else if (! rf )     { rf = i}
                    else if (! fr )     { fr = i}
                    else if (! tandem ) { tandem = i}
                }
            }
        }
        else {
            # get values
            for (n=1; n<=4;n++) { 
                arr[n]=0
                if (\$n) { arr[n] = \$n }
            }

            # sort values on header index and print
            print arr[isize]"\\t"arr[fr]"\\t"arr[rf]"\\t"arr[tandem]"\\t"category
        }
    }' >> ${library}_insertsize_metrics

    # for multiqc channel
    mv good_mapq.out.txt ${library}.good_mapq.insert_size_metrics.txt
    """
}

process picard_metrics {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.21"
    publishDir "${params.outputDir}/stats/picard_alignment_metrics"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)
        path(genome_path)

    output:
        tuple val(params.email), val(library), path('*alignment_summary_metrics.txt'), emit: for_agg

    script:
    """
    genome=\$(ls *.fa 2>/dev/null || ls *.fasta 2>/dev/null)
    picard -Xmx${task.memory.toGiga()}g CollectAlignmentSummaryMetrics \
        --VALIDATION_STRINGENCY SILENT -BS true -R \${genome} \
        -I ${bam} -O ${library}.alignment_summary_metrics.txt
    """
}

process tasmanian {
    label 'medium_cpu'
    tag { library }
    publishDir "${params.outputDir}/stats/tasmanian"
    conda "bioconda::samtools=1.9 bioconda::tasmanian-mismatch=1.0.7"

    errorStrategy { retry < 1 ? 'retry' : 'terminate' }
    maxRetries 1
    memory { retry > 0 ? '8 GB' : '4 GB' }

    input:
        tuple val(library), path(bam), path(bai), val(barcodes)
        path(genome_path)

    output:
        tuple val(params.email), val(library), path('*.csv'), emit: for_agg

    script:
    """
    set +e
    set +o pipefail
    genome=\$(ls *.bwameth.c2t.bwt | sed 's/.bwameth.c2t.bwt//')
    samtools view -q 30 -F 3840 ${bam} | head -n 2000000 | run_tasmanian -r \${genome} > ${library}.csv
    """

}
