nextflow.enable.dsl=2

// include { PATH_TO_TILES_KNOWN } from './modules/path_to_tiles_provided'
include { formatInput_trim_bwamethAlign; mergeAndMarkDuplicates } from './modules/alignment'
include { methylDackel_mbias; methylDackel_extract} from './modules/methylation'
include { aggregate_emseq } from './modules/aggregation.nf'


/* ------------------------ *
 * INCLUDE IN CONFIG FILE!! *
 * ------------------------ */
params.default_dest_path = '/mnt/galaxy/tmp/users'
params.tmp_dir           =  'tmp' // params['other_paths'].tmp_dir
params.path_to_ngs_agg   = '/mnt/bioinfo/prg/ngs-aggregate_results/current'


/* --------------- *
 * INPUT ARGUMENTS *
 * --------------- */
params.email       = 'testing@emseq.neb.com'
params.flowcell    = 'N'
params.genome      = 'undefined'
params.input_glob  =  '*.{1,2}.fastq*'
params.library     = 'undefined'
params.lane        = 'all'
params.tile        = 'all'
params.project     = 'test_project'
params.sample      = 'test_sample'
params.barcode     = ''  

outputDir = 'output_for_now' // params.outdir ?: new File([default_dest_path, "email",flowcell].join(File.separator))


println "Processing " + params.flowcell + "... => " + outputDir
println "Cmd line: $workflow.commandLine"


// This channel will search bam OR fastq. Since we do the matching in bash
// (because if it's bam, we don't want to write to disk but make the fastq
// files on the fly), we exclude the read2. We will include it inside the 
// bash script during alignment.

// modify glob later to 2 channels. Bam and Fastq and concat them!
Channel
    .fromPath(params.input_glob)
    .filter{ it.name =~ /(\.bam|1\.fastq)$/ }  
    .ifEmpty {error "${params.input_glob} could not find any required file."}
    .map{ filename -> 
       [flowcell: params.flowcell, 
        library: params.library,
        input_file: filename,
        lane: params.lane,
        tile: params.tile,
        genome: params.genome ]}
    .set{input_alignment}

 workflow {
    main:
        bwaMeth = formatInput_trim_bwamethAlign( input_alignment )
        markDup = mergeAndMarkDuplicates( bwaMeth.tuple_lib_bam )
        extract = methylDackel_extract( markDup.md_bams )
        mbias   = methylDackel_mbias( markDup.md_bams )

       aggregate_emseq( mbias.mbias_for_aggregate )
}



//.of( params.email, 
//params.library, 
//params.project, 
//params.sample, 
//params.barcode, 
//params.lane, 
//'genome_name_here', // val(genome_name),  
//params.genome, // val(genome_fasta), 
//path(markDup.mardDup.md_bams[1]), 
//path(markDup.mardDup.md_bams[2]), 
///* path(fastqc), 
//path(flagstats), 
//path(idxstats), 
//path(noncontrol_gc), */
//path(bwaMeth.nonconverted_counts), 
//// path(noncontrol_stats),  ???
//path(mbias.mbias_output_tsv)


// lines 93/94 in https://github.com/nebiolabs/seq-shepherd/blob/d255cef0c0f7dc0d2568908808967019180f3f9f/post_run_scripts/em-seq.nf ?
// How is cpus requesting more than 4 cpus working?
