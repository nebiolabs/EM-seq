library = params.library
tmp_dir = params.tmp_dir


process formatInput_trim_bwamethAlign {
    label 'cpus_8'
    tag { [flowcell, library] }
    conda "bwameth seqtk sambamba fastp mark-nonconverted-reads samtools"
    publishDir "${library}/bwameth_align"


    input:
        tuple val(flowcell), 
              val(library),
              path(input_file),
              val(lane),
              val(tile),
              val(genome)

    output:
        path "*.aln.bam", emit: aligned_bams
       // path "*.nonconverted.tsv", emit: nonconverted_counts
        path "*_fastp.json", emit: fastp_log_files
        tuple val(library), path("*.aln.bam"), emit: tuple_lib_bam
        env(barcodes), emit: barcodes

    /* 2 caveats in the following shell script:
    1. Could not include  sed 's/.*BC:Z:\([ACTGN-]*\).*@/\1/' (@ symbol to avoid stop commenting this) 
       hence grep -o replaces it.
    2. bam2fastq_or_fqmerge="samtools fastq -nT BC -@ !{task.cpus} !{input_file}" generates
       a successfull fastq with barcodes in header. HOWEVER, bwameth doesn't know how to parse it.
    */ 
    shell:
    '''
#    set -eo pipefail

    # barcodes=$(samtools view !{input_file} | head -n 10000 | sed 's/.*BC:Z:\\([ACTGN-]*\\).*/\1/' | awk '{tot++; arr[$1]++}END{for (i in arr) { print i"\t"arr[i]/tot*100} }' | sort -k2nr | head -n1 | awk '{print $1}')
    barcodes=$(samtools view !{input_file} | head -n 10000 | grep -o BC:Z:[ACTGN-]* | awk '{tot++; arr[$1]++}END{for (i in arr) { print i"\t"arr[i]/tot*100} }' | sort -k2nr | head -n1 | awk '{print substr($1,6,length($1))}')

    shared_operations() {
        rg_id="@RG\\tID:${fastq_barcode}\\tSM:!{library}"
        bwa_mem_log_filename="!{library}_${fastq_barcode}!{flowcell}_!{lane}_!{tile}.log.bwamem"
        bam_filename="!{library}_${fastq_barcode}_${barcodes}_!{flowcell}_!{lane}_!{tile}.aln.bam"
        inst_name=$(echo $fastq_barcode | sed 's/^@//')   
        trim_polyg=$(echo "${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')
        echo ${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'
    }

#    if $( echo !{input_file} | grep -q "bam$")
#    then
        fastq_barcode=$(samtools view !{input_file} | head -n1 | cut -d ":" -f1);
        shared_operations;
        bam2fastq_or_fqmerge="samtools fastq -n -@ !{task.cpus} !{input_file}"
        # -n in samtools because bwameth needs space not "/" in the header (/1 /2)

#    else  
#        read1=!{input_file}
#        read2=$(ls -l ${read1} | awk '{print $NF}' | 's/1.fastq/2.fastq/')
#        fastq_barcode=$(zcat -f ${read1} | head -n 1 | cut -d ":" -f1)
#        shared_operations;
#        bam2fastq_or_fqmerge="seqtk mergepe <(zcat -f ${read1}) <(zcat -f ${read2})"
#    fi


    ${bam2fastq_or_fqmerge} \
    | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis -j !{library}_fastp.json 2> fastp.stderr \
    | bwameth.py -p -t !{task.cpus} --read-group "${rg_id}" --reference !{params.genome} /dev/stdin 2> ${bwa_mem_log_filename} \
    | sambamba view -t 2 -S -f bam -o ${bam_filename} /dev/stdin 2> sambamba.stderr;
    '''

// THERE IS A ISSUE WITH pysam.calignmentfile.AlignedSegment, where `reference_name` is not recognized`. Unfortunately, 
// in newer versions, there is no such pysam.calignmentfile.AlignedSegment but pysam.AlignedSegment. I will look into this later.
//    ${bam2fastq_or_fqmerge} \
//    | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis -j !{library}_fastp.json 2> fastp.stderr \
//    | bwameth.py -p -t !{task.cpus} --read-group "${rg_id}" --reference !{params.genome} /dev/stdin 2> ${bwa_mem_log_filename} \
//    | mark-nonconverted-reads.py --reference !{params.genome} 2> "!{library}_${fastq_barcode}_!{flowcell}_!{lane}_!{tile}.nonconverted.tsv" \
//    | sambamba view -t 2 -S -f bam -o ${bam_filename} /dev/stdin 2> sambamba.stderr;



}


process bam2fastq {
    conda 'samtools'
    label 'cpus_2'
    tag 'converting unaligned bam onto fastq'
    publishDir "${library}/bam2fastq"

    input:
        path bam_file

    output:
        tuple stdout, path('fq1'), path('fq2')

    shell:
    '''
    samtools fastq !{bam_file} -1 fq1 -2 fq2
    echo "!{bam_file}" | awk -F"/" '{print $NF}' | sed 's/.bam//' | tr -d '\\n'   
    '''
}


process mergeAndMarkDuplicates {
    label 'cpus_8'
    errorStrategy 'retry'
    tag { library }
    publishDir "${library}/markduped_bams", mode: 'copy', pattern: '*.{md.bam}*'
    conda "samtools=1.9 samblaster=0.1.24 sambamba=0.7.0"

    input:
        tuple val(library), file(libraryBam) // from aligned_files.groupTuple()

    output:
        tuple val(library), path('*.md.bam'), path('*.md.bam.bai'), emit: md_bams
        path('*.samblaster'), emit: samblaster_logs

    shell:
    '''
    samtools cat  -b <( find . -name '*.aln.bam' ) \
    | samtools view -h /dev/stdin \
    | samblaster 2> !{library}.log.samblaster \
    | sambamba view -t 2 -l 0 -S -f bam /dev/stdin \
    | sambamba sort --tmpdir=!{params.tmp_dir} -t !{task.cpus} -m 20GB -o !{library}.md.bam /dev/stdin
    '''
}


    // we need to send the same file into multiple channels for downstream 
    // analysis whether they came from aligned bams or fastqs
//    md_bams.into {  md_files_for_mbias; md_files_for_extract; md_files_for_fastqc; 
//                    md_files_for_samstats; md_files_for_picard; 
//                    md_files_for_goleft; md_files_for_picard_gc; md_files_for_samflagstats; 
//                    md_files_for_aggregate; md_files_for_human_reads;
//                 }
