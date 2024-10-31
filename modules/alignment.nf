
process alignReads {
    cpus 8
    tag { library }
    errorStrategy { task.exitStatus == 140 ? 'ignore' : 'terminate' }
    conda "conda-forge::python=3.10 bioconda::bwameth=0.2.7 bioconda::fastp=0.23.4 bioconda::mark-nonconverted-reads=1.2 bioconda::sambamba=1.0 bioconda::samtools=1.19 bioconda::seqtk=1.4"
    publishDir "${library}/bwameth_align"


    input:
        tuple path(input_file1),
              path(input_file2),
              val(genome),
              val(fileType)
    output:
        tuple env(flowcell), val(params.email), val(library), env(barcodes), path("*.nonconverted.tsv"), path("*_fastp.json"), emit: for_agg
        path "*.aln.bam", emit: aligned_bams
        tuple val(library), path("*.nonconverted.tsv"), emit: nonconverted_counts
        tuple val(library), path("*.aln.bam"), path("*.aln.bam.bai"), env(barcodes), emit: bam_files
        env(bam_barcode), emit: seen_barcode

    shell:

    library = input_file1.baseName.replaceFirst(/.fastq|.fastq.gz|.bam/,"").replaceFirst(/_R1$|_1$|.1$/,"")
    '''
    set -eo pipefail   

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

    get_read_length_from_fastq() {
        set +o pipefail
        zcat -f $1 | head -n2 | tail -n1 | wc -c
        set -o pipefail
    }
    get_read_length_from_bam() {
        set +o pipefail
        samtools view $1 | head -n 100 | awk 'BEGIN{m=0}{if (length($10)>m) {m=length($10)}}END{print m}'
        set -o pipefail
    }


    barcodes_from_fastq () {

    set +o pipefail

    zcat -f $1 \
    | awk '{
        if (NR%4==1) {
            split($0, parts, ":"); 
            arr[ parts[ length(parts) ] ]++
        }} END { for (i in arr) {print arr[i]"\\t"i} }' \
    | sort -k1nr | head -n1 | cut -f2 
    # | tr -c "[ACGTN]" "\\t"

    # set -o pipefail

    }    

    case !{fileType} in 
        "fastq_paired_end")
            barcodes=($(barcodes_from_fastq !{input_file1}))
            n_reads=$(get_nreads_from_fastq !{input_file1})
            downsampling="samtools import -u -1 !{input_file1} -2 !{input_file2} -O bam -@!{task.cpus} | samtools view -h -s!{params.downsample_seed}.${n_reads} "
            flowcell=$(flowcell_from_fastq !{input_file1})
            read_length=$(get_read_length_from_fastq !{input_file1})
            ;;
        "bam")
            barcodes=$(samtools view -H !{input_file1} | grep @RG | awk '{for (i=1;i<=NF;i++) {if ($i~/BC:/) {print substr($i,4,length($i))} } }' | head -n1)
            frac_reads=$(samtools view -c !{input_file1} | awk '{frac=!{params.max_input_reads}/$1; if (frac>=1) {print 1} else {split(frac, parts, "."); print parts[2]}}')
            if [ ${frac_reads} -lt 1 ]; then
                ds_suffix="-s !{params.downsample_seed}.${frac_reads}"
            else
                ds_suffix=""
            fi
            downsampling="samtools view -h !{input_file1} ${ds_suffix}"
            flowcell=$(flowcell_from_bam !{input_file1})
            read_length=$(get_read_length_from_bam !{input_file1})
            ;;
        "fastq_single_end")
            barcodes=($(barcodes_from_fastq !{input_file1}))
            n_reads=$(get_nreads_from_fastq !{input_file1})
            downsampling="samtools import -u -s !{input_file1} -@!{task.cpus} | samtools view -h -s!{params.downsample_seed}.${n_reads} "
            flowcell=$(flowcell_from_fastq !{input_file1})
            read_length=$(get_read_length_from_fastq !{input_file1})
            ;;
    esac

    flowcell=$( (echo !{params.flowcell} | grep -q "undefined") && echo "${flowcell}" || echo "!{params.flowcell}")
    downsample=$( (echo !{params.downsample} | grep -q "all_reads") && echo "${downsampling}" || echo " ")

    bwa_mem_log_filename="!{library}_${barcodes}_${flowcell}.log.bwamem"
    bam_filename="!{library}_${barcodes}_${flowcell}.aln.bam"
    rg_id="@RG\\tID:${barcodes}\\tSM:!{library}"
    inst_name=$(echo $barcodes | sed 's/^@//')
    trim_polyg=$(echo "${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')
    echo ${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'
    bam2fastq="| samtools collate -u -@!{task.cpus} /dev/stdin -O | samtools fastq -n -@ !{task.cpus} /dev/stdin"
    # -n in samtools because bwameth needs space not "/" in the header (/1 /2)
 
    set +o pipefail

    eval ${downsampling} ${bam2fastq}  \
    | paste - - - - | tr "\\t" "\\n" \
    | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis -j !{library}_fastp.json 2> fastp.stderr \
    | awk -v rl=${read_length} '{if (NR%4==2 || NR%4==0) {print substr($0,1,rl)} else print $0 }' \
    | bwameth.py -p -t !{task.cpus} --read-group "${rg_id}" --reference !{params.genome} /dev/stdin 2> ${bwa_mem_log_filename} \
    | mark-nonconverted-reads.py --reference !{params.genome} 2> "!{library}_${barcodes}_${flowcell}.nonconverted.tsv" \
    | samtools view -hu /dev/stdin \
    | sambamba sort -l 3 --tmpdir=!{params.tmp_dir} -t !{task.cpus} -m !{task.cpus*8}GB -o ${bam_filename} /dev/stdin
    bam_barcode=$(samtools view ${bam_filename} | head -n1 | cut -f3 -d ":")
    '''
}
process mergeAndMarkDuplicates {
    label 'cpus_8'
    cpus 8
    errorStrategy 'retry'
    tag { library }
    publishDir "${library}/markduped_bams", mode: 'copy', pattern: '*.md.{bam,bai}'
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

    fastq_barcode=$(samtools view !{bam} | head -n1 | cut -d ":" -f1);
    optical_distance=$(echo ${fastq_barcode} | awk '{if ($1~/^M0|^NS|^NB/) {print 100} else {print 2500}}')
    picard -Xmx20g MarkDuplicates TAGGING_POLICY=All OPTICAL_DUPLICATE_PIXEL_DISTANCE=${optical_distance} TMP_DIR=!{params.tmp_dir} CREATE_INDEX=true MAX_RECORDS_IN_RAM=5000000 BARCODE_TAG="RX" ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT I=!{bam} O=!{library}_!{barcodes}.md.bam M=!{library}.markdups_log
    '''
}
