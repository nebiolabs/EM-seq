

process bwameth_align {
    label '8_cpus'
    tag { [flowcell, library] }
    conda "bwameth seqtk sambamba fastp mark-nonconverted-reads"
    publishDir 'log_records/bwameth_align'

    input:
        tuple val(flowcell), 
              val(library),
              path(read1),
              path(read2),
              val(barcode),
              val(lane),
              val(tile) 

    output:
        path "*.aln.bam", emit: aligned_bams
        path "*.nonconverted.tsv", emit: nonconverted_counts
        path "*_fastp.json", emit: fastp_log_files

    shell:
    '''
    bash metadata_extraction_and_alignment.sh \
        !{genome} \
        !{flowcell} \
        !{insert_read1} \
        !{insert_read2} \
        !{library} \
        !{task.cpus} \
        !{tile} \
        !{lane} \
    '''
}


process bam2fastq {
    conda 'samtools'
    label '2_cpus'
    tag 'converting unaligned bam onto fastq'

    input:
        path bam_file

    output:
        tuple path('fq1'), path('fq2'), stdout

    shell:
    '''
    samtools fastq !{bam_file} -1 fq1 -2 fq2
    echo !{bam_file} | awk -F"/" '{print $NF}' | sed 's/.bam//'
    '''
}