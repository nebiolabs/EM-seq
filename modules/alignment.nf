
process alignReads {
    label 'high_cpu'
    tag { library }
    conda "conda-forge::python=3.10 bioconda::bwameth=0.2.7 bioconda::fastp=0.23.4 bioconda::mark-nonconverted-reads=1.2 bioconda::sambamba=1.0 bioconda::samtools=1.19 bioconda::seqtk=1.4"
    publishDir "${params.outputDir}/bwameth_align"
    input:
        tuple path(input_file1),
              path(input_file2),
              val(genome),
              val(fileType)
    output:
        tuple val(params.email), val(library), env(barcodes), path("*.nonconverted.tsv"), path("*.fastp.json"), emit: for_agg
        path "*.aln.bam", emit: aligned_bams
        tuple val(library), path("*.nonconverted.tsv"), emit: nonconverted_counts
        tuple val(library), path("*.aln.bam"), path("*.aln.bam.bai"), env(barcodes), emit: bam_files

    shell:

    library = input_file1.baseName.replaceFirst(/.fastq|.fastq.gz|.bam/,"").replaceFirst(/_R1$|_1$|.1$/,"")
    '''
    get_nreads_from_fastq() {
        zcat -f $1 | grep -c "^+$" \
        | awk '{
            frac=!{params.max_input_reads}/$1; 
            if (frac>=1) {frac=0.999}; 
            split(frac, numParts, "."); print numParts[2]
            }'
    }

    flowcell_from_fastq() {
        set +o pipefail
        zcat -f $1 | head -n1 | cut -d ":" -f3
        set -o pipefail
    }

    flowcell_from_bam(){
        set +o pipefail
        samtools view $1 | head -n1 | cut -d":" -f3
        set -o pipefail
    }

    barcodes_from_fastq () {

        set +o pipefail

        zcat -f $1 \
        | head -n10000 \
        | awk '{
            if (NR%4==1) {
                split($0, parts, ":"); 
                arr[ parts[ length(parts) ] ]++
            }} END { for (i in arr) {print arr[i]"\\t"i} }' \
        | sort -k1nr | head -n1 | cut -f2 
        # | tr -c "[ACGTN]" "\\t"

        set -o pipefail

    }    
    
    get_barcodes_and_rg_line() {
        local file=$1
        local type=$2
        if [ "$type" == "bam" ]; then
            barcodes=$(samtools view -H $file | grep @RG | awk '{for (i=1;i<=NF;i++) {if ($i~/BC:/) {print substr($i,4,length($i))} } }' | head -n1)
            rg_line=$(samtools view -H $file | grep "^@RG" | sed 's/\\t/\\\\t/g')
        else
            barcodes=($(barcodes_from_fastq $file))
            rg_line="@RG\\tID:${barcodes}\\tSM:!{library}\\tBC:${barcodes}"
        fi
    }

    get_frac_reads() {
        local file=$1
        local type=$2
        if [ "!{params.max_input_reads}" == "all_reads" ]; then
            frac_reads=1
        else
            if [ "$type" == "bam" ]; then
                n_reads=$(samtools view -c -F 2304 $file)
            else
                n_reads=$(get_nreads_from_fastq $file)
            fi
            if [ $n_reads -le !{params.max_input_reads} ]; then
                frac_reads=1
            else
                frac_reads=$(echo $n_reads | awk '{!{params.max_input_reads}/$1}')
            fi 
        fi
    }

    case !{fileType} in 
        "fastq_paired_end")
            get_barcodes_and_rg_line !{input_file1} "fastq"
            get_frac_reads !{input_file1} "fastq"
            stream_reads="samtools import -u -1 !{input_file1} -2 !{input_file2}"
            flowcell=$(flowcell_from_fastq !{input_file1})
            ;;
        "bam")
            get_barcodes_and_rg_line !{input_file1} "bam"
            get_frac_reads !{input_file1} "bam"
            stream_reads="samtools view -u -h !{input_file1}"
            flowcell=$(flowcell_from_bam !{input_file1})
            ;;
        "fastq_single_end")
            get_barcodes_and_rg_line !{input_file1} "fastq"
            get_frac_reads !{input_file1} "fastq"
            stream_reads="samtools import -u -s !{input_file1}"
            flowcell=$(flowcell_from_fastq !{input_file1})
            ;;
    esac
    if [ "!{params.flowcell}" == "undefined" ]; then
        flowcell="${flowcell}"
    else
        flowcell="!{params.flowcell}"
    fi

    if [ ${frac_reads} -lt 1 ]; then
        downsample_seed_frac=$(awk -v seed=!{params.downsample_seed} -v frac=${frac_reads} 'BEGIN { printf "%.4f", seed + frac }')
        stream_reads="${stream_reads} | samtools view -u -s ${downsample_seed_frac}"
    fi

    base_outputname="!{library}_${barcodes}_${flowcell}"
   
    inst_name=$(samtools view !{input_file1} | head -n 1 | cut -d ":" -f 1)
    trim_polyg=$(echo "${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')
    echo ${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'
    bam2fastq="| samtools collate -f -r 100000 -u /dev/stdin -O | samtools fastq -n  /dev/stdin"
    # -n in samtools because bwameth needs space not "/" in the header (/1 /2)
     
    eval ${stream_reads} ${bam2fastq} \
    | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis -j "${base_outputname}.fastp.json" 2> fastp.stderr \
    | bwameth.py -p -t !{Math.max(1,(task.cpus*7).intdiv(8))} --read-group "${rg_line}" --reference !{params.genome} /dev/stdin 2> "${base_outputname}.log.bwamem" \
    | mark-nonconverted-reads.py --reference !{params.genome} 2> "${base_outputname}.nonconverted.tsv" \
    | samtools view -u /dev/stdin \
    | sambamba sort -l 3 --tmpdir=!{params.tmp_dir} -t !{Math.max(1,task.cpus.intdiv(8))} -m !{(task.memory.toGiga()*3).intdiv(4)}GB -o "${base_outputname}.aln.bam" /dev/stdin
    '''
}

process mergeAndMarkDuplicates {
    label 'high_cpu'
    tag { library }
    publishDir "${params.outputDir}/markduped_bams", mode: 'copy', pattern: '*.md.{bam,bai}'
    conda "bioconda::picard=3.1 bioconda::samtools=1.19"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes) 

    output:
        tuple val(library), path('*.md.bam'), path('*.md.bai'), val(barcodes), emit: md_bams
        tuple val( params.email ), val(library), path('*.md.bam'), path('*.md.bai'), emit: for_agg
        path('*.markdups_log'), emit: log_files

    shell:
    '''
    set +o pipefail
    inst_name=$(samtools view !{bam} | head -n1 | cut -d ":" -f1);
    set -o pipefail
    
    optical_distance=$(echo ${inst_name} | awk '{if ($1~/^M0|^NS|^NB/) {print 100} else {print 2500}}')

    picard -Xmx!{task.memory.toGiga()}g MarkDuplicates \
        --TAGGING_POLICY All \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE ${optical_distance} \
        --TMP_DIR !{params.tmp_dir} \
        --CREATE_INDEX true \
        --MAX_RECORDS_IN_RAM 5000000 \
        --BARCODE_TAG "RX" \
        --ASSUME_SORT_ORDER coordinate \
        --VALIDATION_STRINGENCY SILENT \
        -I !{bam} \
        -O !{library}_!{barcodes}.md.bam \
        -M !{library}.markdups_log
    '''
}


process bwa_index {
    // publishDir on wherever the genome is located.
    // Hence genome param should be full path
    label 'low_cpu'
    tag { genome }
    conda "bioconda::samtools=1.19 bioconda::bwa=0.7.18"

    shell:
    '''
    samtools faidx !{params.genome}
    bwa index !{params.genome}
    '''
}