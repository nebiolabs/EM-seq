
process alignReads {
    //label 'cpus_8'
    cpus 8
    tag { flowcell }
    conda "python=3.10 bwameth seqtk sambamba fastp mark-nonconverted-reads samtools"
    publishDir "${params.flowcell}/${library}/bwameth_align"


    input:
        tuple val(flowcell), 
              path(input_file),
              val(lane),
              val(tile),
              val(genome)

    output:
        path "*.aln.bam", emit: aligned_bams
        tuple val(library), path("*.nonconverted.tsv"), emit: nonconverted_counts
        path "*_fastp.json", emit: fastp_log_files
        tuple val(library), path("*.aln.bam"), path("*.aln.bam.bai"), env(barcodes), emit: bam_files
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

    library = input_file.baseName.replaceFirst(/.fastq|.bam/,"")

    '''
    # set -eo pipefail

    barcodes=$(samtools view !{input_file} | head -n 10000 | grep -o BC:Z:[ACTGN-]* | awk '{tot++; arr[$1]++}END{for (i in arr) { print i"\t"arr[i]/tot*100} }' | sort -k2nr | head -n1 | awk '{print substr($1,6,length($1))}')
    
    shared_operations() {
        bwa_mem_log_filename="!{library}_${fastq_barcode}!{flowcell}_!{lane}_!{tile}.log.bwamem"
        bam_filename="!{library}_${fastq_barcode}_${barcodes}_!{flowcell}_!{lane}_!{tile}.aln.bam"
        rg_id="@RG\\tID:${fastq_barcode}\\tSM:!{library}"
        inst_name=$(echo $fastq_barcode | sed 's/^@//')
        trim_polyg=$(echo "${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')
        echo ${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'
    }
    
    downsampling=$([[ !{params.max_input_reads} != -1 ]] && echo  " | seqtk sample -s!{params.downsample_seed} /dev/stdin !{params.max_input_reads}" || echo "")
    # add reformat.sh from bbmap

    if $( echo !{input_file} | grep -q ".bam$")
    then
        fastq_barcode=$(samtools view -F 2304 !{input_file} | head -n1 | cut -d ":" -f1);
        shared_operations;
        bam2fastq_or_fqmerge=" samtools collate -@ 2 !{input_file} -O | samtools fastq -n -@ 2 /dev/stdin"
        # -n in samtools because bwameth needs space not "/" in the header (/1 /2)

    else  
        read1=!{input_file}
        read2=$(readlink ${read1} | awk '{print $NF}' | sed -E 's/1.fastq/2.fastq/')
        if [ ${#read2} -eq 0 ]; then 
            echo "fastq files HAVE TO END WITH 1.fastq and 2.fastq" && exit
        fi
        fastq_barcode=$(zcat -f ${read1} | head -n 1 | cut -d ":" -f1)
        shared_operations;
        bam2fastq_or_fqmerge="seqtk mergepe <(zcat -f ${read1}) <(zcat -f ${read2})"
    fi

    eval ${bam2fastq_or_fqmerge} ${downsampling} \
    | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis -j !{library}_fastp.json 2> fastp.stderr \
    | awk '{if (NR%4==2 || NR%4==0) {print substr($0,1,!{params.read_length})} else print $0 }' \
    | bwameth.py -p -t !{task.cpus} --read-group "${rg_id}" --reference !{params.genome} /dev/stdin 2> ${bwa_mem_log_filename} \
    | mark-nonconverted-reads.py --reference !{params.genome} 2> "!{library}_${fastq_barcode}_!{params.flowcell}_!{lane}_!{tile}.nonconverted.tsv" \
    | samtools view -hb /dev/stdin \
    | sambamba sort -l 3 --tmpdir=!{params.tmp_dir} -t !{task.cpus} -m !{task.cpus*8}GB -o ${bam_filename} /dev/stdin
    bam_barcode=$(samtools view ${bam_filename} | head -n1 | cut -f3 -d ":")
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
    publishDir "${params.flowcell}/${library}/markduped_bams", mode: 'copy', pattern: '*.{md.bam}*'
    conda "picard samtools"

    input:
        tuple val(library), path(bam), path(bai), val(barcodes) 

    output:
        tuple val(library), path('*.md.bam'), path('*.md.bai'), val(barcodes), emit: md_bams
        path('*.markdups_log'), emit: log_files
        //tuple val(library), path()

    shell:
    '''
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
    conda "bioconda::samtools"

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
