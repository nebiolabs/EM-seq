
process gc_bias {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.22"

    input:
        tuple val(library), path(bam), path(bai)
        tuple path(genome_fa), path(genome_fai)
    output:
        tuple val(library), path('*gc_metrics'), emit: for_agg

    script:
    """
    samtools view -H ${bam} | grep "^@SQ" \
    | grep -v "plasmid_puc19\\|phage_lambda\\|phage_Xp12\\|phage_T4\\|EBV\\|chrM" \
    | awk -F":|\\t" '{print \$3"\\t"0"\\t"\$5}' > include_regions.bed

    samtools view -h -L include_regions.bed ${bam} | \
    picard -Xmx${task.memory.toGiga()}g CollectGcBiasMetrics \
        --IS_BISULFITE_SEQUENCED true --VALIDATION_STRINGENCY SILENT \
        -I /dev/stdin -O ${library}.gc_metrics -S ${library}.gc_summary_metrics \
        -R ${genome_fa} --CHART /dev/null
    """
}

process idx_stats {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::samtools=1.22"

    input:
        tuple val(library), path(bam), path(bai)

    output:
        tuple val(library), path("*idxstat"), emit: for_agg

    script:
    """
    samtools idxstats ${bam} > ${library}.idxstat
    """
}

process flag_stats {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::samtools=1.22"

    input:
        tuple val(library), path(bam), path(bai)

    output:
        tuple val(library), path("*flagstat"), emit: for_agg

    script:
    """
    samtools flagstat -@${task.cpus} ${bam} > ${library}.flagstat
    """
}

process fastqc {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::fastqc=0.11.8"

    input:
        tuple val(library), path(bam), path(bai)

    output:
        tuple val(library), path('*_fastqc.zip'), emit: for_agg

    shell:
    """
    fastqc -f bam ${bam}
    """
}

process insert_size_metrics {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.22"

    input:
        tuple val(library), path(bam), path(bai)

    output:
        tuple val(library), path('*_metrics'), emit: for_agg
        tuple val(library), path('*.good_mapq.insert_size_metrics.txt'), emit: high_mapq_insert_size_metrics

    script:
    """

    good_mapq="${library}.tmp.good_mapq.\$\$"
    bad_mapq="${library}.tmp.bad_mapq.\$\$"
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

    # Explicitly remove fifos before Nextflow collects outputs
    rm -f "\$good_mapq" "\$bad_mapq"
    """
}

process picard_metrics {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.22"

    input:
        tuple val(library), path(bam), path(bai)
        tuple path(genome_fa), path(genome_fai)

    output:
        tuple val(library), path('*alignment_summary_metrics.txt'), emit: for_agg

    script:
    """
    picard -Xmx${task.memory.toGiga()}g CollectAlignmentSummaryMetrics \
        --VALIDATION_STRINGENCY SILENT -BS true -R ${genome_fa} \
        -I ${bam} -O ${library}.alignment_summary_metrics.txt
    """
}

process tasmanian {
    label 'medium_cpu'
    tag { library }
    conda "bioconda::samtools=1.22 bioconda::tasmanian-mismatch=1.0.9"
    if (params.publishOutput) {
        publishDir params.publishOutput ? "{params.outputDir}/${task.process}_logs" : null, mode: 'symlink'
    }

    errorStrategy { retry < 1 ? 'retry' : 'terminate' }
    maxRetries 1
    memory { retry > 0 ? '16 GB' : '8 GB' }

    input:
        tuple val(library), path(bam), path(bai)
        tuple path(genome_fa), path(genome_fai)

    output:
        tuple val(library), path('*.csv'), emit: for_agg

    script:
    """
    set +e
    set +o pipefail
    samtools view -q 30 -F 3840 ${bam} | head -n 2000000 | run_tasmanian -r ${genome_fa} > ${library}.csv
    """

}
