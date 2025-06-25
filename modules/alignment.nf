process alignReads {
    label 'high_cpu'
    tag { library }
    conda "conda-forge::python=3.10 bioconda::bwameth=0.2.7 bioconda::fastp=0.26 bioconda::mark-nonconverted-reads=1.2 bioconda::samtools=1.22 bioconda::seqtk=1.4 bioconda::gatk4=4.6.2.0 conda-forge::gawk=5.3.1"
    publishDir "${params.outputDir}/bwameth_align"
    memory {
        try { 
            def fileSize = input_file1.size() / (1024 * 1024 * 1024)
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
        path(input_file1),
        path(input_file2),
        val(fileType)
        path(genome)

    output:
        tuple val(library), path("*.fastp.json"), emit: fastp_reports
        tuple val(library), path("*.nonconverted.tsv"), emit: nonconverted_counts
        tuple val(library), path("*.aln.bam"), path("*.aln.bam.bai"), emit: aligned_bams
        tuple val(library), path("*.metadata.bam"), emit: metadata_bams

    script:
    """
    echo "Input file size: ${input_file1.size() / (1024 * 1024 * 1024)} GB"
    echo "Memory allocated: ${task.memory}" 

    # Determine the genome index
    genome=\$(ls *.bwameth.c2t.bwt | sed 's/.bwameth.c2t.bwt//')

    # Define helper functions
    flowcell_from_fastq() {
        set +o pipefail
        zcat -f "\$1" | head -n1 | cut -d ":" -f3
        set -o pipefail
    }

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

    # assumes that the barcodes are in the last part of the read name 
    # extracts the most frequent barcode in the first 10k reads
    barcodes_from_fastq() {
        zcat -f "\$1" | awk 'NR%4==1' | head -n10000 | \
        grep -o '[GCATN+-]\{6,\}' | \
        sort | uniq -c | sort -nr | head -n1 | \
        awk '{gsub(/N+/, "-", \$2); print \$2}'
    }

    # Determine barcodes and read group line
    get_rg_line() {
        set +o pipefail
        local file=\$1
        local type=\$2
        if [ "\$type" == "bam" ]; then
            rg_line=\$(samtools view -H \$file | grep "^@RG" | sed 's/\\t/\\\\t/g' | head -n1)
        else
            barcodes=(\$(barcodes_from_fastq \$file))
            rg_line="@RG\\tID:\${barcodes}\\tSM:${library}\\tBC:\${barcodes}"
        fi
        set -o pipefail
    }

    get_frac_reads() {
        local file=\$1
        local type=\$2
        if [ "${params.max_input_reads}" == "all_reads" ]; then
            frac_reads=1
        else
            if [ "\$type" == "bam" ]; then
                n_reads=\$(samtools view -c -F 2304 \$file)
            elif [ "\$file" == "fastq_paired_end" || "\$file" == "fastq_single_end" ]; then
                n_reads=\$(zcat -f \$file | grep -c "^+\$")
            else
                echo "Error: Unsupported file type \$type" >&2
                exit 1
            fi
            if [ \$n_reads -le ${params.max_input_reads} ]; then
                frac_reads=1
            else
                frac_reads=\$(echo \$n_reads | awk '{print ${params.max_input_reads}/\$1}')
            fi 
        fi
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


   case ${fileType} in 
        "fastq_paired_end")
            get_rg_line ${input_file1} "fastq"
            get_frac_reads ${input_file1} "fastq"
            stream_reads="samtools import -u -1 ${input_file1} -2 ${input_file2}"
            flowcell=\$(flowcell_from_fastq ${input_file1})
            ;;
        "bam")
            get_rg_line ${input_file1} "bam"
            get_frac_reads ${input_file1} "bam"
            stream_reads="samtools view -u -h ${input_file1}"
            flowcell=\$(flowcell_from_bam ${input_file1})
            ;;
        "fastq_single_end")
            get_rg_line ${input_file1} "fastq"
            get_frac_reads ${input_file1} "fastq"
            stream_reads="samtools import -u -s ${input_file1}"
            flowcell=\$(flowcell_from_fastq ${input_file1})
            ;;
    esac
    if [ "${params.flowcell}" == "undefined" ]; then
        flowcell="\${flowcell}"
    else
        flowcell="${params.flowcell}"
    fi

    if [ "\${frac_reads}" != "1" ]; then
        downsample_seed_frac=\$(awk -v seed=${params.downsample_seed} -v frac=\${frac_reads} 'BEGIN { printf "%.4f", seed + frac }')
        stream_reads="\${stream_reads} | samtools view -u -s \${downsample_seed_frac}"
    fi

    base_outputname="${library}_\${barcodes}_\${flowcell}"
   
    set +o pipefail
    inst_name=\$(samtools view ${input_file1} | head -n 1 | cut -d ":" -f 1)
    set -o pipefail

    trim_polyg=\$(echo "\${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')
    echo \${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'
    
    # hard clips all but the first base to create minimal metadata bam
    metadata_tee="tee >(gatk ClipReads -I /dev/stdin -O \"\${base_outputname}.metadata.bam\" --clip-representation HARDCLIP_BASES -CT 2-10000)"
    
    bam2fastq="| \${metadata_tee} | samtools collate -f -r 100000 -u /dev/stdin | samtools fastq -n  /dev/stdin"
    # -n in samtools because bwameth needs space not "/" in the header (/1 /2)

 
    eval \${stream_reads} \${bam2fastq} \
    | fastp --stdin --stdout -l 2 -Q \${trim_polyg} --interleaved_in --overrepresentation_analysis -j "\${base_outputname}.fastp.json" 2> fastp.stderr \
    | bwameth.py -p -t ${Math.max(1,(task.cpus*7).intdiv(8))} --read-group "\${rg_line}" --reference \${genome} /dev/stdin 2> "\${base_outputname}.log.bwamem" | reheader_sam /dev/stdin \
    | mark-nonconverted-reads.py --reference \${genome} 2> "\${base_outputname}.nonconverted.tsv" \
    | samtools view -u /dev/stdin \
    | samtools sort -T ${params.tmp_dir}/samtools_sort_tmp -@ ${Math.max(1,task.cpus.intdiv(8))} \
       -m ${(task.memory.toGiga()*5).intdiv(8)}G --write-index \
       -o "\${base_outputname}.aln.bam##idx##\${base_outputname}.aln.bam.bai" /dev/stdin 

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
        tuple val(library), path('*.markdups_log'), emit: for_agg

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
        -O ${library}_${barcodes}.md.bam \
        -M ${library}.markdups_log
    """
}

process genome_index {
    /* Combined process for genome indexing.
     * Creates both bwameth aligner index and samtools faidx index.
     * Attempts to link existing indices. If there is no index
     * we download it from the provided URL.
     * If no index and no URL, User will have to debug.
     */

    label 'low_cpu'
    tag { genome_basename }
    conda "bioconda::samtools=1.22 bioconda::bwameth=0.2.7"
    
    output:
    path "bwameth_index/*.{fa,fai,amb,ann,bwt,pac,sa,c2t}", emit: aligner_files
    tuple path("genome_index/*.fa"), path("genome_index/*.fai"), emit: genome_index

    script:
    genome_basename = file(params.path_to_genome_fasta).baseName
    """
    # Create output directories
    mkdir -p bwameth_index genome_index
    
    real_genome_file="\$(basename ${params.path_to_genome_fasta})"
    
    # Handle genome file (download if URL or link if local)
    if [ ! -f "bwameth_index/\${real_genome_file}" ]; then
        # if the reference .fa file is a url, not a local path
        if [[ "${params.path_to_genome_fasta}" =~ ^https?:// ]]; then
            echo "Trying to download the reference"
            if ! curl -f -o "bwameth_index/\${real_genome_file}" ${params.path_to_genome_fasta}; then
                echo "Error: Failed to download \${params.path_to_genome_fasta}" >&2
                exit 1
            fi
        else
            # Link local files
            ln -sf "\$(dirname ${params.path_to_genome_fasta})/\${real_genome_file}"* bwameth_index/
        fi
    fi
    
    # Create bwameth index
    cd bwameth_index
    if [ ! -f "\${real_genome_file}.bwameth.c2t.bwt" ]; then
        echo "Creating bwameth index for \${real_genome_file}"
        bwameth.py index \${real_genome_file}
    else
        echo "Bwameth index files already exist for \${real_genome_file}"
    fi
    cd ..
    
    # Create samtools faidx index
    cd genome_index
    # Copy the fasta file to genome_index directory
    cp "../bwameth_index/\${real_genome_file}" "\${real_genome_file}"
    
    # Check if index exists adjacent to the original FASTA file (using the full path parameter)
    if [ -f "${params.path_to_genome_fasta}.fai" ]; then
        echo "Found existing genome index: ${params.path_to_genome_fasta}.fai"
        ln -sf "${params.path_to_genome_fasta}.fai" "\${real_genome_file}.fai"
    else
        echo "No existing genome index found. Creating new index for \${real_genome_file}"
        samtools faidx "\${real_genome_file}"
    fi
    
    # Verify the index file exists and is not empty
    if [ ! -s "\${real_genome_file}.fai" ]; then
        echo "Error: Failed to create or link genome index file"
        exit 1
    fi
    
    echo "Genome indices ready:"
    echo "  - Bwameth index: bwameth_index/\${real_genome_file}.bwameth.c2t.*"
    echo "  - Samtools index: genome_index/\${real_genome_file}.fai"
    """
}

process touchFile {   
    input:
        val filename

    script:
    """
    touch ${filename}
    """
}
