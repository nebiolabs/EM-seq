

process formatInput_trim_bwamethAlign {
    label 'cpus_8'
    tag { [flowcell, library] }
    conda "bwameth seqtk sambamba fastp mark-nonconverted-reads"
    publishDir 'output/bwameth_align'


    input:
        tuple val(flowcell), 
              val(library),
              //path(input_file),
              val(input_file),
              val(lane),
              val(tile) 

    // `input_file` = *1.fastq or *.bam (not 2.fastq)

    output:
        path "*.aln.bam", emit: aligned_bams
        path "*.nonconverted.tsv", emit: nonconverted_counts
        // path "*_fastp.json", emit: fastp_log_files

    shell:
    '''
    set -eo pipefail

    if $(grep -q "2.fastq" !{input_file}); 
    then 
        echo "nothing wrong with this process. Skipping read2.fastq... "
        exit 0
    fi

    shared_operations() {
        rg_id="@RG\\tID:${fastq_barcode}\\tSM:!{library}"
        bwa_mem_log_filename="!{library}_${fastq_barcode}!{flowcell}_!{lane}_!{tile}.log.bwamem"
        bam_filename="!{library}_${fastq_barcode}_!{flowcell}_!{lane}_!{tile}.aln.bam"
        inst_name=$(echo $fastq_barcode | sed 's/^@//')   
        trim_polyg=$(echo "${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')
        echo ${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'
    }

    if $( echo !{input_file} | grep -q "bam$")
    then
        fastq_barcode=$(samtools view !{input_file} | head -n1 | cut -d ":" -f1);
        shared_operations;
        bam2fastq_or_fqmerge="samtools fastq -@ !{task.cpus}"
    else  
        read1=!{input_file}
        read2=$(echo !{input_file} | sed 's/1.fastq/2.fastq/')
        fastq_barcode=$(zcat -f !{input_file} | head -n 1 | cut -d ":" -f1)
        shared_operations;
        bam2fastq_or_fqmerge="seqtk mergepe <(zcat -f ${read1}) <(zcat -f ${read2})"
    fi


    ${bam2fastq_or_fqmerge} \
    | fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis -j !{library}_fastp.json 2> fastp.stderr \
    | bwameth.py -p -t !{task.cpus} --read-group ${rg_id} --reference !{params.genome} /dev/stdin 2> ${bwa_mem_log_filename} \
    | mark-nonconverted-reads.py 2> "!{library}_${fastq_barcode}_!{flowcell}_!{lane}_!{tile}.nonconverted.tsv" \
    | sambamba view -t 2 -S -f bam -o ${bam_filename} /dev/stdin 2> sambamba.stderr;
    '''
}


process bam2fastq {
    conda 'samtools'
    label 'cpus_2'
    tag 'converting unaligned bam onto fastq'
    publishDir 'output/bam2fastq'

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