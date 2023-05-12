tmp_dir = params.tmp_dir

process alignReads {
    label 'cpus_8'
    tag { flowcell }
    conda "python=3.10 bwameth seqtk sambamba fastp mark-nonconverted-reads samtools samblaster"
    publishDir "${library}/bwameth_align"


    input:
        tuple val(flowcell), 
              path(input_file),
              val(lane),
              val(tile),
              val(genome)

    output:
        path "*.aln.bam", emit: aligned_bams
        path "*.nonconverted.tsv", emit: nonconverted_counts
        path "*_fastp.json", emit: fastp_log_files
        tuple env(library), path("*.aln.bam"), path("*aln.bam.bai"), env(barcodes), emit: bam_files

    /* 2 caveats in the following shell script:
    1. Could not include  sed 's/.*BC:Z:\([ACTGN-]*\).*@/\1/' (@ symbol to avoid stop commenting this) 
       hence grep -o replaces it.
    2. bam2fastq_or_fqmerge="samtools fastq -nT BC -@ !{task.cpus} !{input_file}" generates
       a successfull fastq with barcodes in header. HOWEVER, bwameth doesn't know how to parse it.
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

    if $( echo !{input_file} | grep -q ".bam$")
    then
        fastq_barcode=$(samtools view !{input_file} | head -n1 | cut -d ":" -f1);
        shared_operations;
        #bam2fastq_or_fqmerge="samtools fastq -n -@ !{task.cpus} !{input_file}"
        bam2fastq_or_fqmerge="samtools fastq -@ !{task.cpus} !{input_file} | paste - - - - | sort | uniq | tr \"\\t\" \"\\n\" | sed 's/\\(@.*\\)\\/./\\1/'"
        # -n in samtools because bwameth needs space not "/" in the header (/1 /2)

    else  
        read1=!{input_file}
        read2=$(ls -l ${read1} | awk '{print $NF}' | 's/1.fastq/2.fastq/')
        if [ ${#read2} -eq 0 ]; then 
            echo "fastq files HAVE TO END WITH 1.fastq and 2.fastq" && exit
        fi
        fastq_barcode=$(zcat -f ${read1} | head -n 1 | cut -d ":" -f1)
        shared_operations;
        bam2fastq_or_fqmerge="seqtk mergepe <(zcat -f ${read1}) <(zcat -f ${read2})"
    fi


    ${bam2fastq_or_fqmerge} \
    | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis -j !{library}_fastp.json 2> fastp.stderr \
    | bwameth.py -p -t !{task.cpus} --read-group "${rg_id}" --reference !{params.genome} /dev/stdin 2> ${bwa_mem_log_filename} \
    | mark-nonconverted-reads.py --reference !{params.genome} 2> "!{library}_${fastq_barcode}_!{params.flowcell}_!{lane}_!{tile}.nonconverted.tsv" \
    | samtools sort -m20G -@4 | samtools view -h \
    | samblaster 2> !{library}.log.samblaster | samtools view -bh !{library}_${barcodes}_!{flowcell}.aln.bam && \
    samtools index -@4 -o !{library}_${barcodes}_!{flowcell}.aln.bam

    #| sambamba sort -m !{task.memory} -t 2 --tmpdir=!{params.tmp_dir} /dev/stdin \
    '''
}

// If distinguishing PCR or sequencing duplicates doesn't make any difference
// We don't have to use picard and hence we don't neeed optical distance.
// Also, to write the least we need to the disk, I piped samblaster in the previous step.
process mergeAndMarkDuplicates {
    label 'cpus_8'
    errorStrategy 'retry'
    tag { library }
    publishDir "${library}/markduped_bams", mode: 'copy', pattern: '*.{md.bam}*'
    conda "samtools=1.9 samblaster=0.1.24 sambamba=0.7.0"

    input:
        tuple val(library), file(libraryBam), val(barcodes) 

    output:
        tuple val(library), path('*.md.bam'), path('*.md.bam.bai'), val(barcodes), emit: md_bams
        path('*.samblaster'), emit: samblaster_logs

    shell:
    '''
    optical_distance=$(echo !{barcodes} | awk '{if ($1~/^M0|^NS|^NB/) {print 100} else {print 2500}}')
    #samtools cat  -b <( find . -name '*.aln.bam' ) \
    | samtools view -h /dev/stdin \
    | samblaster 2> !{library}.log.samblaster \
    | sambamba view -t 2 -l 0 -S -f bam /dev/stdin \
    | sambamba sort --tmpdir=!{params.tmp_dir} -t !{task.cpus} -m 20GB -o !{library}_!{barcodes}.md.bam /dev/stdin
    '''
}




