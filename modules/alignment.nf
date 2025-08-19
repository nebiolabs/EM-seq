process enough_reads {
    tag {library}
    conda "bioconda::samtools=1.22"

    input:
        tuple val(library), 
              path(bam)

        output:
            tuple val(library), path(bam), path("*passes_or_fails.txt") 

        script:
        """
        full_path=\$(realpath ${bam})

        if [ \$(stat -c%s \${full_path}) -lt 100 ]; then
            passes_or_fails="fail"
        else
            passes_or_fails="pass"
        fi

        echo -e "$library\\t\${passes_or_fails}" > ${library}_passes_or_fails.txt
        """
}

process send_email {
    
    input:
        file libraries

    script:
    """
    touch tmp
    for f in ${libraries}
    do
        cat \$f | awk '{print \$1"<br>"}' >> tmp
    done
    libs=\$(cat tmp)

    sendmail -t <<EOF
    To: ${params.email}
    Subject: File Read Check
    Content-Type: text/html

    <html>
      <body>
        <p>The following libraries:<br> <strong>\${libs}</strong> do not have enough reads. <br> Continuing with other libraries. </p>
      </body>
    </html>
    EOF
    """
}


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
        tuple val(library),
              path(bam)
        val(genome)

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

    get_frac_reads() {
        local file=\$1
        if [ "${params.max_input_reads}" == "all_reads" ]; then
            frac_reads=1
        else
            n_reads=\$(samtools view -c -F 2304 \$file)
            if [ \$n_reads -le ${params.max_input_reads} ]; then
                frac_reads=1
            else
                frac_reads=\$(echo \$n_reads | awk '{print ${params.max_input_reads}/\$1}')
            fi 
        fi
        export frac_reads
}

    reheader_sam() {
      local input_file="\$1"

      cat "\$input_file" | \
      awk 'BEGIN{
        header_ids = "@HD @SQ @RG @CO" # exclude @pg
        split(header_ids, headers_arr, " ")
        flag=0;
      } {
        if (\$1~/^@/) {gsub(/\\\\t/,"\\t",\$0); id = substr(\$1,1,3); arr[id] = arr[id]"\\n"\$0}
        else {
          if (flag==0) {
        for (id in headers_arr){
          printf "%s", arr[headers_arr[id]]
        }
        flag=1;
        print ""
          }
          print \$0
        }
      }' | tail -n +2
    }
    
    get_rg_line ${bam}
    get_frac_reads ${bam}
    stream_reads="samtools view -u -h ${bam}"
    flowcell=\$(flowcell_from_bam ${bam})

    if [ "${params.flowcell}" == "undefined" ]; then
        flowcell="\${flowcell}"
    else
        flowcell="${params.flowcell}"
    fi

    if [ "\${frac_reads}" != "1" ]; then
        downsample_seed_frac=\$(awk -v seed=${params.downsample_seed} -v frac=\${frac_reads} 'BEGIN { printf "%.4f", seed + frac }')
        stream_reads="\${stream_reads} | samtools view -u -s \${downsample_seed_frac}"
    fi

    set +o pipefail
    inst_name=\$(samtools view ${bam} | head -n 1 | cut -d ":" -f 1)
    set -o pipefail

    trim_polyg=\$(echo "\${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')
    echo \${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'

    bam2fastq="| samtools collate -f -r 100000 -u /dev/stdin -O | samtools fastq -n  /dev/stdin"
    # -n in samtools because bwameth needs space not "/" in the header (/1 /2)

    eval \${stream_reads} \${bam2fastq} \
    | fastp --stdin --stdout -l 2 -Q \${trim_polyg} --interleaved_in --overrepresentation_analysis -j "${library}.fastp.json" 2> fastp.stderr \
    | bwameth.py -p -t ${Math.max(1,(task.cpus*7).intdiv(8))} --read-group "\${rg_line}" --reference ${genome} /dev/stdin 2> "${library}.log.bwamem" | reheader_sam /dev/stdin \
    | mark-nonconverted-reads.py --reference ${genome} 2> "${library}.nonconverted.tsv" \
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
