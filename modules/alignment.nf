
process alignReads {
    //label 'cpus_8'
    cpus 8
    tag { flowcell }
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
        tuple env(flowcell), val(library), path("*.aln.bam"), path("*.aln.bam.bai"), env(barcodes), emit: bam_files
        env(bam_barcode), emit: seen_barcode

    /* 2 caveats in the following shell script:
    1. Could not include  sed 's/.*BC:Z:\([ACTGN-]*\).*@/\1/' (@ symbol to avoid stop commenting this) 
       hence grep -o replaces it.
    2. bam2fastq_or_fqmerge="samtools fastq -nT BC -@ !{task.cpus} !{input_file}" generates
       a successfull fastq with barcodes in header. HOWEVER, bwameth doesn't know how to parse it.
    3. We first get rid of secondary/supp alignments, in case the bam was already aligned to a genome.
       Then we collate it for the same reason, as mapped files can be coordinate sorted.
    */ 
    shell:

    library = input_file1.baseName.replaceFirst(/.fastq|.fastq.gz|.bam/,"").replaceFirst(/_R1$|_1$|.1$/,"")
    // def bamFile = (fileType == "bam") ? "${library}.bam" : null
    '''
    # set -eo pipefail   

    get_nreads_from_fastq() {
        zcat -f $1 | grep -c "^+$" \
        | awk '{
            frac=!{params.max_input_reads}/$1; 
            if (frac>=1) {frac=0.999}; 
            split(frac, numParts, "."); print numParts[2]
            }'
    }

    flowcell_from_fastq() {
        zcat -f $1 | head -n1 | cut -d ":" -f3
    }
    flowcell_from_bam(){
        samtools view $1 | head -n1 | cut -d":" -f3
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

    set -o pipefail

    }    

    case !{fileType} in 
        "fastq_paired_end")
            barcodes=($(barcodes_from_fastq !{input_file1}))
            n_reads=$(get_nreads_from_fastq !{input_file1})
            downsampling="samtools import -u -1 !{input_file1} -2 !{input_file2} -O bam -@!{task.cpus} | samtools view -h -s!{params.downsample_seed}.${n_reads} "
            flowcell=$(flowcell_from_fastq !{input_file1})
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
            ;;
        "fastq_single_end")
            barcodes=($(barcodes_from_fastq !{input_file1}))
            n_reads=$(get_nreads_from_fastq !{input_file1})
            downsampling="samtools import -u -s !{input_file1} -O bam -@!{task.cpus} | samtools view -h -s!{params.downsample_seed}.${n_reads} "
            flowcell=$(flowcell_from_fastq !{input_file1})
            ;;
    esac

    flowcell=$( (echo !{params.flowcell} | grep -q "undefined") && echo "${flowcell}" || echo "!{params.flowcell}")

    bwa_mem_log_filename="!{library}_${barcodes}_!{flowcell}_!{lane}_!{tile}.log.bwamem"
    bam_filename="!{library}_${barcodes}_!{flowcell}_!{lane}_!{tile}.aln.bam"
    rg_id="@RG\\tID:${barcodes}\\tSM:!{library}"
    inst_name=$(echo $barcodes | sed 's/^@//')
    trim_polyg=$(echo "${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')
    echo ${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'
    bam2fastq="| samtools collate -@!{task.cpus} /dev/stdin -O | samtools fastq -n -@ !{task.cpus} /dev/stdin"
    # -n in samtools because bwameth needs space not "/" in the header (/1 /2)
 
    set +o pipefail

    eval ${downsampling} ${bam2fastq}  \
    | tee >(paste - - - - | sed -n '1~2!p' | tr "\\t" "\\n"  | gzip > !{library}_!{lane}_!{tile}.fq.gz) \
    | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis -j !{library}_fastp.json 2> fastp.stderr \
    | awk '{if (NR%4==2 || NR%4==0) {print substr($0,1,!{params.read_length})} else print $0 }' \
    | bwameth.py -p -t !{task.cpus} --read-group "${rg_id}" --reference !{params.genome} /dev/stdin 2> ${bwa_mem_log_filename} \
    | mark-nonconverted-reads.py --reference !{params.genome} 2> "!{library}_${barcodes}_!{params.flowcell}_!{lane}_!{tile}.nonconverted.tsv" \
    | samtools view -hu /dev/stdin \
    | sambamba sort -l 3 --tmpdir=!{params.tmp_dir} -t !{task.cpus} -m !{task.cpus*8}GB -o ${bam_filename} /dev/stdin
    bam_barcode=$(samtools view ${bam_filename} | head -n1 | cut -f3 -d ":")
    # zcat !{library}_*.fq.gz | gzip > !{library}.fq.gz  #for metadata_fastq_channel
    '''
}

// If distinguishing PCR or sequencing duplicates doesn't make any difference
// We don't have to use picard and hence we don't neeed optical distance.
// Also, to write the least we need to the disk, I piped samblaster in the previous step.
process mergeAndMarkDuplicates {
    label 'cpus_8'
    cpus 8
    errorStrategy 'retry'
    tag { library }
    publishDir "${library}/markduped_bams", mode: 'copy', pattern: '*.{md.bam}*'
    conda "bioconda::picard=3.1 bioconda::samtools=1.19"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes) 

    output:
        tuple val(library), path('*.md.bam'), path('*.md.bai'), val(barcodes), emit: md_bams
        tuple val( params.email ), val(library), path('*.md.bam'), path('*.md.bai'), emit: for_agg
        path('*.markdups_log'), emit: log_files
        //tuple val(library), path()

    shell:
    '''
    set +o pipefail

    fastq_barcode=$(samtools view !{bam} | head -n1 | cut -d ":" -f1);
    optical_distance=$(echo ${fastq_barcode} | awk '{if ($1~/^M0|^NS|^NB/) {print 100} else {print 2500}}')
    picard -Xmx20g MarkDuplicates TAGGING_POLICY=All OPTICAL_DUPLICATE_PIXEL_DISTANCE=${optical_distance} TMP_DIR=!{params.tmp_dir} CREATE_INDEX=true MAX_RECORDS_IN_RAM=5000000 BARCODE_TAG="RX" ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT I=!{bam} O=!{library}_!{barcodes}.md.bam M=!{library}.markdups_log
    '''
}

process find_soft_clips {
    /* I can not see this variable in the output of picard::CollectAlignmentSummaryMetrics.
     * see alignment_metrics.rb or the picard tool for more information. Perhaps, older
     * versions of the picard tool will contain this information. 
     */ 
    label 'cpus_8'
    tag { library }
    publishDir "${params.flowcell}/${library}/softclips", more: 'copy', pattern: '*.{cigar_stats.tsv}*'
    conda "bioconda::samtools=1.9"

    input: 
        tuple val(library), path(bam), path(bai), val(barcodes)

    output:
        tuple val(library), val(barcodes), emit: cigar_stats

    // bam file should be *md.bam
    shell:
    '''
        sam=$(echo !{bam} | sed 's/bam/cigar_stats.tsv/')
        samtools view !{bam} | awk '{}'
    '''
}

//process downsample {
//    label 'cpus_8'
//    tag { library }
//    publishDir "${params.flowcell}/${library}/downsampled"
//    conda "agbiome:bbtools"
//
//    input:
//        tuple val(flowcell), 
//              path(input_file),
//              val(lane),
//              val(tile),
//              val(genome)
//
//    output:
//        tuple val(flowcell), 
//              path(downsampled_reads),
//              val(lane),
//              val(tile),
//              val(genome), emit: downsampled_reads_and_metadata
//
//    shell:
//    '''
//        if $( echo !{input_file} | grep -q ".bam$"); then
//            output_file=$(echo !{input_file} | sed 's/.bam/.dspl.bam/')
//            reformat.sh in=!{input_file} out=!{input_file} reads=4000000 sampleseed=42
//        else
//            read1=!{input_file}
//            read2=$(ls -l ${read1} | awk '{print $NF}' | 's/1.fastq/2.fastq/')
//                        
//        fi
//    '''
//
//}
