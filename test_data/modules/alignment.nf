

process bwameth_align {
    label 'cpus_8'
    tag { [flowcell, library] }
    conda "bwameth seqtk sambamba fastp mark-nonconverted-reads"
    publishDir 'output/bwameth_align'

    input:
        tuple val(flowcell), 
              val(library),
              path(read1),
              path(read2),
              val(lane),
              val(tile) 

    output:
        path "*.aln.bam", emit: aligned_bams
        path "*.nonconverted.tsv", emit: nonconverted_counts
        // path "*_fastp.json", emit: fastp_log_files

    shell:
    '''
    fastq_barcode=$(zcat -f !{read1} | head -n 1 | sed -r 's/.*://')
    rg_id="@RG\\tID:${fastq_barcode}\\tSM:!{library}"
    bwa_mem_log_filename="!{library}_${fastq_barcode}!{flowcell}_!{lane}_!{tile}.log.bwamem"
    bam_filename="!{library}_${fastq_barcode}_!{flowcell}_!{lane}_!{tile}.aln.bam"
    inst_name=$(zcat -f !{read1} | head -n 1 | cut -f 1 -d ':' | sed 's/^@//')
    
    trim_polyg=$(echo "${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH/) {print "--trim_poly_g"} else {print ""}}')

    echo ${trim_polyg} | awk '{ if (length(\$1)>0) { print "2-color instrument: poly-g trim mode on" } }'

    echo "${rg_id} ${inst_name} ${bam_filename} ${bwa_mem_log_filename}"

    seqtk mergepe <(zcat -f !{read1}) <(zcat -f !{read2}) \
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