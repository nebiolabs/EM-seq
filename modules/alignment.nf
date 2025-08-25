process fastp {
	tag "${library}"
    label 'process_single'
    conda "bioconda::samtools=1.21 bioconda::fastp=0.23.4"
    publishDir "${params.outputDir}/stats/fastp", mode: 'copy', pattern: "*json"

    input:
        tuple val(library), path(bam)
    
    output:
        tuple val(library), path("*.1.trimmed.fastq"), path("*.2.trimmed.fastq"), emit: trimmed_fastq
        tuple val(library), path("${library}.fastp.json"), emit: fastp_json

    script:
    """
    set +o pipefail
    inst_name=\$(samtools view ${bam} | head -n 1 | cut -d ":" -f 1)
    set -o pipefail

    trim_polyg=\$(echo "\${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')
    echo \${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'
    
    samtools fastq -n ${bam} > ${library}.fastq 
    fastp --interleaved_in --in1 ${library}.fastq \
                    -l 2 -Q \${trim_polyg} \
                    --overrepresentation_analysis \
                    -j "${library}.fastp.json" \
                    --split ${params.fastq_split_count} \
                    --out1 ${library}.1.trimmed.fastq \
                    --out2 ${library}.2.trimmed.fastq
    """
}

process alignReads {
    label 'high_cpu'
    tag { library }
    conda "conda-forge::python=3.10 bioconda::bwameth=0.2.7 bioconda::fastp=0.26 bioconda::mark-nonconverted-reads=1.2 bioconda::samtools=1.22"
    publishDir "${params.outputDir}/bwameth_align", mode: 'symlink'
    memory {
        try { 
            def fileSize = fastq.size() / (1024 * 1024 * 1024)
            if (fileSize < 1.8) return '64 GB'
            else if (fileSize < 6.5) return '128 GB'
            else return '256 GB'
        }
        catch (Exception _e) {
            return '128 GB'  // Default memory if size cannot be determined
        }
    }

    input:
        tuple val(library), val(chunk_name), path(read1), path(read2)
        val(genome_fa)

    output:
        tuple val(library), path("*.aln.bam"), path("*.aln.bam.bai"), emit: bam_files
        tuple val(library), path("*.nonconverted.tsv"), emit: nonconverted_counts

    script:
    """
    echo "Input file size: ${read1.size() / (1024 * 1024 * 1024)} GB"
    echo "Memory allocated: ${task.memory}" 

    bwameth.py -p -t ${Math.max(1,(task.cpus*7).intdiv(8))} --reference ${genome_fa} ${read1} ${read2} 2> "${library}.log.bwamem" \
    | mark-nonconverted-reads.py --reference ${genome_fa} 2> "${library}.nonconverted.tsv" \
    | samtools view -u /dev/stdin \
    | samtools sort -T ${params.tmp_dir}/samtools_sort_tmp -@ ${Math.max(1,task.cpus.intdiv(8))} \
       -m ${(task.memory.toGiga()*5).intdiv(8)}G --write-index \
       -o "${chunk_name}.aln.bam##idx##${chunk_name}.aln.bam.bai" /dev/stdin


    """
}

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
