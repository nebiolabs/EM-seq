nextflow.enable.dsl=2
include { inferReadLength } from './modules/input_processing'
/* ------------------------ *
 * INCLUDE IN CONFIG FILE!! * //TODO: why !! these are NEB internal parameters?
 * ------------------------ */
params.default_dest_path = '/mnt/galaxy/tmp/users'
params.tmp_dir           =  '/tmp' // params['other_paths'].tmp_dir
params.path_to_ngs_agg   = '/mnt/bioinfo/prg/ngs-aggregate_results/current'

/* --------------- *
 * INPUT ARGUMENTS *
 * --------------- */
params.email       = 'me@example.com'
params.flowcell    = '?'
params.genome      = 'unknown'
params.input_glob  =  '*.{1,2}.fastq*'

Channel.fromPath(params.input_glob, checkIfExists: true) // this just checks that the files exist, actual channel is set up later

params.read_length = inferReadLength(params.input_glob)

params.email = 'whoami'.execute().text.trim() + '@' + 'hostname -d'.execute().text.trim()
println "Using email: ${params.email}"

//TODO: try to infer the flowcell using the illumina read naming convention, prompt the user if not possible
params.flowcell    = params.flowcell ?: channel.prompt("Please enter the flowcell: ")
params.genome      = ''
if (!file(params.genome).exists()) {
    error "Genome path '${params.genome}' does not exist"
}
params.lane        = 'all'
params.tile        = 'all'
params.project     = 'test_project'
params.sample      = 'test_sample'
params.barcode     = '' //TODO: infer from bam/fastq, cannot be a parameter to the whole workflow
params.workflow    = 'EM-seq' //TODO: embed workflow version here
params.outputDir   = 'em-seq_output'
params.min_mapq    = 20 // for methylation assessment.
params.max_input_reads = -1 // number of reads to use from each library default(-1) to skip downsampling
params.downsample_seed = 42

// include { PATH_TO_TILES_KNOWN } from './modules/path_to_tiles_provided' //TODO: what is this?
include { alignReads; mergeAndMarkDuplicates }                                                          from './modules/alignment'
include { methylDackel_mbias; methylDackel_extract }                                                    from './modules/methylation'
include { gc_bias; idx_stats; flag_stats; fastqc; insert_size_metrics; picard_metrics; tasmanian }      from './modules/compute_statistics'
//TODO: use collect_multiple_metrics instead of picard_metrics run separately were possible
include { aggregate_emseq }                                                                             from './modules/aggregation'


println "Processing " + params.flowcell + "... => " + params.outputDir
println "Cmd line: $workflow.commandLine"


// This channel will search bam OR fastq. Since we do the matching in bash
// (because if it's bam, we don't want to write to disk but make the fastq
// files on the fly), we exclude the read2. We will include it inside the 
// bash script during alignment. So glob for fastq SHOULD ONLY INCLUDE read1
// and NOT read2. E.g. *[\._-]1.fastq

// modify glob later to 2 channels. Bam and Fastq and concat them!
inputChannel = Channel
    .fromPath(params.input_glob)
    .filter{ it.name =~ /(\.bam|\.1\.fastq(?:\.gz)?)$/ }
    .ifEmpty {error "${params.input_glob} could not find any required file."}
    .map{ filename -> 
       [flowcell: params.flowcell,
        input_file: filename,
        lane: params.lane,
        tile: params.tile,
        genome: params.genome ]}


 workflow {
    main:
        // process files 
        alignedReads = alignReads( inputChannel )
        markDup      = mergeAndMarkDuplicates( alignedReads.bam_files )
        extract      = methylDackel_extract( markDup.md_bams )
        mbias        = methylDackel_mbias( markDup.md_bams )

        // collect statistics
        gcbias       = gc_bias( markDup.md_bams )
        idxstats     = idx_stats( markDup.md_bams )
        flagstats    = flag_stats( markDup.md_bams )
        fastqc       = fastqc( markDup.md_bams ) // All reads go in here. Good and Bad mapq.
        insertsize   = insert_size_metrics( markDup.md_bams ) 
        metrics      = picard_metrics( markDup.md_bams )
        mismatches   = tasmanian ( markDup.md_bams )

        // Channel for aggregation 
        markDup.md_bams
        .join( alignedReads.nonconverted_counts )
        .join( gcbias.for_agg )
        .join( idxstats.for_agg )
        .join( flagstats.for_agg )
        .join( fastqc.for_agg )
        .join( insertsize.for_agg )
        .join( mismatches.for_agg )
        .join( mbias.for_agg)
        .join( metrics.for_agg)
        .set{ aggregation_Channel }
        
        // aggregate_Channel.view()
         aggregate_emseq( aggregation_Channel ) 
}
