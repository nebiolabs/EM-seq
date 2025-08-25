process alignReads {
    label 'high_cpu'
    tag { library }
    conda "conda-forge::python=3.10 bioconda::bwameth=0.2.7 bioconda::fastp=0.26 bioconda::mark-nonconverted-reads=1.2 bioconda::samtools=1.22"
    publishDir "${params.outputDir}/bwameth_align", mode: 'symlink'
    memory {
        try { 
            def fileSize = bam.size() / (1024 * 1024 * 1024)
            if (fileSize < 1.8) return '64 GB'
            else if (fileSize < 6.5) return '128 GB'
            else return '256 GB'
        }
        catch (Exception _e) {
            return '128 GB'  // Default memory if size cannot be determined
        }
    }

    input:
        tuple val(library), path(bam)
        val(genome_fa)

    output:
        tuple val(library), path("*.fastp.json"), emit: fastp_reports
        tuple val(library), path("*.nonconverted.tsv"), emit: nonconverted_counts
        tuple val(library), path("*.aln.bam"), path("*.aln.bam.bai"), emit: bam_files

    script:
    """
    echo "Input file size: ${bam.size() / (1024 * 1024 * 1024)} GB"
    echo "Memory allocated: ${task.memory}" 

    # Define helper functions
    flowcell_from_bam(){
        set +o pipefail
        # check this is NOT mgi:
        fc=\$(samtools view \$1 | head -n1 | cut -f1)
        if echo \$fc | grep -q ":"; then
            echo "\$fc" | cut -d":" -f3
        elif echo \$fc | grep -q "L"; then
            echo "\$fc" | cut -d "L" -f1
        fi
        set -o pipefail
    }

    # Determine barcodes and read group line
    get_rg_line() {
        set +o pipefail
        local file=\$1
        barcodes=\$(samtools view -H \$file | grep @RG | awk '{for (i=1;i<=NF;i++) {if (\$i~/BC:/) {print substr(\$i,4,length(\$i))} } }' | head -n1)
        rg_line=\$(samtools view -H \$file | grep "^@RG" | sed 's/\\t/\\\\t/g' | head -n1)
        set -o pipefail
        # Export variables to parent scope
        export barcodes
        export rg_line
    }

    get_rg_line ${bam}
    stream_reads="samtools view -u -h ${bam}"
    flowcell=\$(flowcell_from_bam ${bam})

    if [ "${params.flowcell}" == "undefined" ]; then
        flowcell="\${flowcell}"
    else
        flowcell="${params.flowcell}"
    fi

    set +o pipefail
    inst_name=\$(samtools view ${bam} | head -n 1 | cut -d ":" -f 1)
    set -o pipefail

    trim_polyg=\$(echo "\${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')
    echo \${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'

    samtools view -u -h ${bam} \
    | samtools fastq -n  /dev/stdin \
    | fastp --stdin --stdout -l 2 -Q \${trim_polyg} --interleaved_in --overrepresentation_analysis -j "${library}.fastp.json" 2> fastp.stderr \
    | bwameth.py -p -t ${Math.max(1,(task.cpus*7).intdiv(8))} --read-group "\${rg_line}" --reference ${genome_fa} /dev/stdin 2> "${library}.log.bwamem" \
    | mark-nonconverted-reads.py --reference ${genome_fa} 2> "${library}.nonconverted.tsv" \
    | samtools view -u /dev/stdin \
    | samtools sort -T ${params.tmp_dir}/samtools_sort_tmp -@ ${Math.max(1,task.cpus.intdiv(8))} \
       -m ${(task.memory.toGiga()*5).intdiv(8)}G --write-index \
       -o "${library}.aln.bam##idx##${library}.aln.bam.bai" /dev/stdin


    """
}

process mergeAndMarkDuplicates {
    label 'high_cpu'
    tag { library }
    publishDir "${params.outputDir}/markduped_bams", mode: 'copy', pattern: '*.md.{bam,bai}'
    publishDir "${params.outputDir}/stats/markdups", mode: 'copy', pattern: '*log'
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.22"

    input:
        tuple val(library), path(bam), path(bai)
    output:
        tuple val(library), path('*.md.bam'), path('*.md.bai'), emit: md_bams
        path('*.markdups_log'), emit: log_files
        tuple val(library), path('*.markdups_log'), emit: log

    script:
    """
    set +o pipefail
    inst_name=\$(samtools view ${bam} | head -n1 | cut -d ":" -f1);
    set -o pipefail

    optical_distance=\$(echo \${inst_name} | awk '{if (\$1~/^M0|^NS|^NB/) {print 100} else {print 2500}}')

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
        -I ${bam} \
        -O ${library}.md.bam \
        -M ${library}.markdups_log
    """
}
