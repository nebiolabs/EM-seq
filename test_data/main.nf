nextflow.enable.dsl=2

include { bwameth_align; bam2fastq } from './modules/alignment'


/* --------------- *
 * INPUT ARGUMENTS *
 * --------------- */

flowcell           = params.flowcell
genome             = params.genome
params.tmp_dir     = '/tmp'
outputPath         = params.outdir     = 'output'
fastq_mode         = params.fastq_mode = 'run_fastqs'
params.input_files = params.input_files ?: '*.{1,2}.fastq*'

println "Processing " + flowcell + "... => " + outputPath
println "Cmd line: $workflow.commandLine"


/* --------------- *
 *  FASTQ or BAM   *
 * --------------- */

workflow harmonize_input_type {
    take:
        input_files

    main:
        if (input_files =~ /fastq/) {

            Channel.fromFilePairs(input_files)
            .map{ lib, reads -> [flowcell: flowcell, library:lib, insert_read1:reads[0], 
                            insert_read2:reads[1], lane:'all', tile:'all' 
                            ]
            }.set{fq_set_channel}
            
        }
        else if (input_files =~ /bam/) {

            bam_files = Channel.fromPath(input_files)
            bam2fastq(bam_files)
            .map{ lib, read1, read2 -> [flowcell: flowcell, library:lib, insert_read1:read1, 
                                            insert_read2:read2, lane:'all', tile:'all' 
                                        ]
            }.set{fq_set_channel}
            
        }
        else {
            println 'WRONG EXTENSION FILE. CHOOSE BAM OR FASTQ'
            System.exit(1)
        }

    emit:
        fq_set_channel

}

workflow {
    fq_set_channel = harmonize_input_type(params.input_files)

    bwameth_align(fq_set_channel)

}
//  workflow{
//      take:
//          fq_set_channel
//  
//      main:
//          //fullWorkflowFromFastq()
//          bwameth_align(fq_set_channel).view()
//  }