process mergeAndMarkDuplicates {
    label 'cpus_8'
    tag { library }
    publishDir "${params.outputDir}/markduped_bams", mode: 'copy', pattern: '*.md.{bam,bai}'
    publishDir "${params.outputDir}/stats/markdups", mode: 'copy', pattern: '*log'
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.22"
    memory {
    try {
        def totalFileSize = bams.collect { it.size() }.sum() / (1024 * 1024 * 1024)
        def calculatedMemory = totalFileSize < 2 ? '32.GB' : totalFileSize < 8 ? '64.GB' : '128.GB'
        return params.max_memory ?: calculatedMemory
    }
    catch (Exception _e) {
        return params.max_memory ?: '64.GB'
    }
}

    input:
        tuple val(library), path(bams), path(bais)
    output:
        tuple val(library), path('*.md.bam'), path('*.md.bai'), emit: md_bams
        tuple val(library), path('*.markdups_log'), emit: log
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions
        tuple val("${task.process}"), val('picard'), eval('picard MarkDuplicates --version 2>&1 | cut -f 2 -d ":"'), topic: versions

    script:
    """
    set +o pipefail
    inst_name=\$(samtools view ${bams[0]} | head -n1 | cut -d ":" -f1);
    set -o pipefail

    optical_distance=\$(echo \${inst_name} | awk '{if (\$1~/^M0|^NS|^NB/) {print 100} else {print 2500}}')

    samtools merge --threads ${task.cpus} ${library}.bam ${bams}

    picard -Xmx${task.memory.toGiga().intdiv(2)}g MarkDuplicates \
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
