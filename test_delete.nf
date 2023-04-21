nextflow.enable.dsl=2

params.input_files = params.input_files ?: '*.{1,2}.fastq*'


process bam2fastq {
    label 'use_8_cpus'
    conda 'samtools'
    publishDir '!{bam_file}/fastq_files'

    input:
        path bam_file

    output:
        tuple stdout, path('fq1'), path('fq2')
    
    shell:
    '''
    samtools fastq !{bam_file} -1 fq1 -2 fq2
    echo !{bam_file} | awk -F"/" '{print $NF}' | sed 's/.bam//'
    '''    
}


process aver {
    input:
        tuple val(library), path(read1), path(read2), val(barcode), val(lane), val(tile)
        // tuple val(library), path(fastq1), path(fastq2)

    output:
        stdout

    shell:
    '''
    echo !{read1} and !{library}
    '''    
}


workflow {

    files = Channel.fromPath( params.input_files )
    
    bam2fastq(files)
    .map{ lib, read1, read2 -> [ 
                       library:lib,
                       insert_read1:read1, 
                       insert_read2:read2, 
                       barcode:'N', 
                       lane:'all', 
                       tile:'all' ]
    }.set{fq_set_channel}
   
   aver(fq_set_channel).view() 
}
