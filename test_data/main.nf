nextflow.enable.dsl=2

include { bwameth_align } from './modules/alignment'


# --------------- #
# INPUT ARGUMENTS #
# --------------- #

flowcell           = params.flowcell
genome             = params.genome
params.tmp_dir     = '/tmp'
outputPath         = params.outdir     = 'output'
fastq_mode         = params.fastq_mode = 'run_fastqs'
params.input_files = params.input_files ?: '*.{1,2}.fastq*'

println "Processing " + flowcell + "... => " + outputPath

# --------------- #
#  FASTQ or BAM   #
# --------------- #

if (params.input_files =~ /fastq/) {

    fastq_pairs = Channel.fromFilePairs(fastq_glob)
    .map{ lib, reads -> [flowcell: flowcell, library:lib, insert_read1:reads[0], 
                       insert_read2:reads[1], barcode:'N', lane:'all', tile:'all' 
                       ]
    }.set{fq_set_channel}
       
}
else if (params.input_files =~ /bam/) {

    bam_files = Channel.fromPath(params.input_files)
        .map{ lib, read1, read2 -> [flowcell: flowcell, library:lib, insert_read1:read1, 
                                    insert_read2:read2, barcode:'N', lane:'all', tile:'all' 
                                ]
    }.set{fq_set_channel}
    
}
else {
    println 'WRONG EXTENSION FILE. CHOOSE BAM OR FASTQ'
    System.exit(1)
}


workflow{
    take:
        fq_set_channel

    main:
        fullWorkflowFromFastq()
}