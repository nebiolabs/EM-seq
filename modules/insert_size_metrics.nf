process insert_size_metrics {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.22"
    publishDir "${params.outputDir}/stats/insert_size"

    input:
        tuple val(library), path(bam), path(bai)

    output:
        tuple val(library), path('*_metrics'), emit: for_agg
        tuple val(library), path('*good_mapq.insert_size_metrics.txt'), emit: high_mapq

    script:
    """

    good_mapq_pipe="${library}.good_mapq"
    bad_mapq_pipe="${library}.bad_mapq"
    trap "rm -f \$good_mapq_pipe \$bad_mapq_pipe" EXIT

    mkfifo "\$good_mapq_pipe"
    mkfifo "\$bad_mapq_pipe"

    samtools view -h -q 20 -U "\$bad_mapq_pipe" ${bam} > "\$good_mapq_pipe" &
    samtools_pid=\$!

    picard -Xmx${task.memory.toGiga()}g CollectInsertSizeMetrics \
        --INCLUDE_DUPLICATES --VALIDATION_STRINGENCY SILENT \
        -I "\$good_mapq_pipe" -O ${library}.good_mapq.insert_size_metrics.txt \
        --MINIMUM_PCT 0 -H /dev/null &
    picard_good_mapq_pid=\$!

    picard -Xmx${task.memory.toGiga()}g CollectInsertSizeMetrics \
        --INCLUDE_DUPLICATES --VALIDATION_STRINGENCY SILENT \
        -I "\$bad_mapq_pipe" -O ${library}.bad_mapq.insert_size_metrics.txt \
        --MINIMUM_PCT 0 -H /dev/null &
    picard_bad_mapq_pid=\$!

    wait \$samtools_pid
    wait \$picard_good_mapq_pid
    wait \$picard_bad_mapq_pid

    # Combine and annotate insert size metrics for aggregation
    grep -B 1000 '^insert_size' ${library}.good_mapq.insert_size_metrics.txt | grep -v "insert_size" > ${library}.insertsize_metrics
    echo -e "insert_size\tAll_Reads.fr_count\tAll_Reads.rf_count\tAll_Reads.tandem_count\tcategory" >> ${library}.insertsize_metrics

    for file in ${library}.good_mapq.insert_size_metrics.txt ${library}.bad_mapq.insert_size_metrics.txt; do
        if [[ "\$file" == *good_mapq* ]]; then
            category=">=20"
        else
            category="<20"
        fi
        grep -A1000 '^insert_size' "\$file" | awk -v cat="\$category" 'NR>1 {print \$1"\t"\$2"\t"\$3"\t"\$4"\t"cat}' >> ${library}.insertsize_metrics
    done

    """
}
