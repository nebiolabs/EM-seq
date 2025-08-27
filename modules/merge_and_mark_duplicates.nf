process mergeAndMarkDuplicates {
    label 'high_cpu'
    tag { library }
    publishDir "${params.outputDir}/markduped_bams", mode: 'copy', pattern: '*.md.{bam,bai}'
    publishDir "${params.outputDir}/stats/markdups", mode: 'copy', pattern: '*log'
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.22"

    input:
        tuple val(library), path(bams), path(bais)
    output:
        tuple val(library), path('*.md.bam'), path('*.md.bai'), emit: md_bams
        path('*.markdups_log'), emit: log_files
        tuple val(library), path('*.markdups_log'), emit: log

    script:
    """
    set +o pipefail
    inst_name=\$(samtools view ${bams[0]} | head -n1 | cut -d ":" -f1);
    set -o pipefail

    optical_distance=\$(echo \${inst_name} | awk '{if (\$1~/^M0|^NS|^NB/) {print 100} else {print 2500}}')

    samtools merge ${library}.bam ${bams}

    picard -Xmx${task.memory.toGiga()}g MarkDuplicates \
        --TAGGING_POLICY All \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE \${optical_distance} \
        --TMP_DIR ${params.tmp_dir} \
        --CREATE_INDEX true \
        --MAX_RECORDS_IN_RAM 5000000 \
        --BARCODE_TAG "RX" \
        --ASSUME_SORT_ORDER coordinate \
        --VALIDATION_STRINGENCY SILENT \
        --ADD_PG_TAG_TO_READS false \
        -I ${library}.bam \
        -O ${library}.md.bam \
        -M ${library}.markdups_log
    """
}
