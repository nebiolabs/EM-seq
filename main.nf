nextflow.enable.dsl=2

/* ------------------------ *
 * INCLUDE IN CONFIG FILE!! *
 * ------------------------ */
params.default_dest_path = '/mnt/galaxy/tmp/users'
params.tmp_dir           =  'tmp' // params['other_paths'].tmp_dir
params.path_to_ngs_agg   = '/mnt/bioinfo/prg/ngs-aggregate_results/current'
params.input_files_type  = 'This needs to be specified as bam or fastq'
params.input_file_delim  = ',' // the provided input file with filenames

/* --------------- *
 * INPUT ARGUMENTS *
 * --------------- */
params.read_length = 76
params.email       = 'testing@emseq.neb.com'
params.flowcell    = 'N'
params.genome      = 'undefined'
params.input_glob  =  '*.{1,2}.fastq*'
params.lane        = 'all'
params.tile        = 'all'
params.project     = 'test_project'
params.sample      = 'test_sample'
params.barcode     = ''
params.workflow    = 'Automated EM-seq'
params.develop_mode  = false // When set to true, workflow will not exit early. 
outputDir = 'output_for_now' // params.outdir ?: new File([default_dest_path, "email",flowcell].join(File.separator))
params.min_mapq = 20 // for methylation assessment.
params.max_input_reads = -1 // default is not downsampling 
params.downsample_seed = 42

// include { PATH_TO_TILES_KNOWN } from './modules/path_to_tiles_provided'
include { alignReads; mergeAndMarkDuplicates }                                                          from './modules/alignment'
include { methylDackel_mbias; methylDackel_extract }                                                    from './modules/methylation'
include { gc_bias; idx_stats; flag_stats; fast_qc; insert_size_metrics; picard_metrics; tasmanian }     from './modules/compute_statistics'
include { aggregate_emseq }                                                                             from './modules/aggregation'


println "Processing " + params.flowcell + "... => " + outputDir
println "Cmd line: $workflow.commandLine"


// Reading file names and barcodes from a file makes it compatible to
// both internal and external use. At NEB we will need to pass the 
// values from the sample sheet into the file that is used in this wf.
// tile and lane should have been merged priot to executing this pipeline

inputChannel = Channel
    .fromPath(params.input_file)
    .splitCsv(sep: params.input_file_delim)

    // input_file content -> filename.bam, fake_file, barcode1-barcode2
    // fake file is to match the input set cardinality to alignment, as there is no 
    // conditional input in nextflow (https://groups.google.com/g/nextflow/c/_ygESaTlCXg/m/zBEmF1glJgAJ)
    if (params.input_files_type == 'bam')
        inputChannel = inputChannel.map { row -> tuple( file(row[0]), file('not_a_file'), row[1] )}

    // input_file content -> read1.fastq, read2.fastq, barcode1-barcode2
    else if (params.input_files_type == 'fastq' || params.input_files_type == 'fastq.gz')
        inputChannel = InputChannel.map { row -> tuple( file(row[0]), file(row[1]), row[2] )}

    // include lane, tile, flowcell and genome data into channel in first process. 
    /*.map{  -> 
       [flowcell: params.flowcell,
        input_file: filename,
        lane: params.lane,
        tile: params.tile,
        genome: params.genome ]} */


// This channel will search bam OR fastq. Since we do the matching in bash
// (because if it's bam, we don't want to write to disk but make the fastq
// files on the fly), we exclude the read2. We will include it inside the 
// bash script during alignment. So glob for fastq SHOULD ONLY INCLUDE read1
// and NOT read2. E.g. *[\._-]1.fastq

// modify glob later to 2 channels. Bam and Fastq and concat them!
Channel
    .fromPath(params.input_glob)
    .filter{ it.name =~ /(\.bam|1\.fastq)$/ }
    .ifEmpty {error "${params.input_glob} could not find any required file."}
    .map{ filename -> 
       [flowcell: params.flowcell,
        input_file: filename,
        lane: params.lane,
        tile: params.tile,
        genome: params.genome ]}
    .set{ inputChannel }


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
        fastqc       = fast_qc( markDup.md_bams ) // All reads go in here. Good and Bad mapq.
        insertsize   = insert_size_metrics( markDup.md_bams ) 
        metrics      = picard_metrics( markDup.md_bams )
        mismatches   = tasmanian ( markDup.md_bams )


        // Channel for aggregation
        alignedReads.for_agg.groupTuple(by: [0, 1])
         .join( markDup.for_agg.groupTuple(by: [0,1]), by: [0,1] )
         .join( gcbias.for_agg.groupTuple(by: [0,1]), by: [0,1]  )
         .join( idxstats.for_agg.groupTuple(by: [0,1]), by: [0,1]  )
         .join( flagstats.for_agg.groupTuple(by: [0,1]), by: [0,1]  )
         .join( fastqc.for_agg.groupTuple(by: [0,1]), by: [0,1]  )
         .join( insertsize.for_agg.groupTuple(by: [0,1]), by: [0,1]  )
         .join( mismatches.for_agg.groupTuple(by: [0,1]), by: [0,1]  )
         .join( mbias.for_agg.groupTuple(by: [0,1]), by: [0,1] )
         .join( metrics.for_agg.groupTuple(by: [0,1]), by: [0,1] )
         .set{ aggregation_Channel }


        // aggregation_Channel.view()
        aggregate_emseq( aggregation_Channel ) 
        
}