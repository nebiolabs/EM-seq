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

    grep -h -A1000 '^insert_size' ${library}.good_mapq.insert_size_metrics.txt ${library}.bad_mapq.insert_size_metrics.txt | awk 'BEGIN{flag=0} {
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
    }' >> ${library}.insertsize_metrics

    """
}
