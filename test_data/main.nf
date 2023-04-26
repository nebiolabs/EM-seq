nextflow.enable.dsl=2

// include { PATH_TO_TILES_KNOWN } from './modules/path_to_tiles_provided'
include { formatInput_trim_bwamethAlign; bam2fastq } from './modules/alignment'


/* ------------------------ *
 * INCLUDE IN CONFIG FILE!! *
 * ------------------------ */
default_dest_path = '/mnt/galaxy/tmp/users'
// tmp_dir           =  params['other_paths'].tmp_dir
path_to_ngs_agg   = '/mnt/bioinfo/prg/ngs-aggregate_results/current'


/* --------------- *
 * INPUT ARGUMENTS *
 * --------------- */
params.flowcell    = 'N'
params.genome      = 'N'
params.tmp_dir     = '/tmp'
params.input_glob =  '*.{1,2}.fastq*'

outputDir = 'output_for_now' // params.outdir ?: new File([default_dest_path, "email",flowcell].join(File.separator))


println "Processing " + params.flowcell + "... => " + outputDir
println "Cmd line: $workflow.commandLine"


// This channel will search bam OR fastq. Since we do the matching in bash
// (because if it's bam, we don't want to write to disk but make the fastq
// files on the fly), we exclude the read2. We will include it inside the 
// bash script during alignment.
Channel
    .fromPath(params.input_glob)
    .map{ it -> it.toString()}
    .filter(  ~/.*1.bam/ )
    .ifEmpty {error "${params.input_glob} could not find any required file."}
    .set{input_listoffiles_Channel}

workflow {
    main:
        // formatInput_trim_bwamethAlign(input_listoffiles_Channel)
        formatInput_trim_bwamethAlign( [params.flowcell, "some_library", input_listoffiles_Channel, 'some_lane', 'some_tile'] )


        // bwameth_align(fq_set_channel)

}




// lines 93/94 in https://github.com/nebiolabs/seq-shepherd/blob/d255cef0c0f7dc0d2568908808967019180f3f9f/post_run_scripts/em-seq.nf ?
// How is cpus requesting more than 4 cpus working?