nextflow.enable.dsl=2

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
params.lane        = 'all'
params.tile        = 'all'
params.project     = 'test_project'
params.sample      = 'test_sample'
params.barcode     = ''
params.develop_mode  = false // When set to true, workflow will not exit early. 
outputDir = 'output_for_now' // params.outdir ?: new File([default_dest_path, "email",flowcell].join(File.separator))

// include { PATH_TO_TILES_KNOWN } from './modules/path_to_tiles_provided'
include { alignReads; mergeAndMarkDuplicates }                                                          from './modules/alignment'
include { methylDackel_mbias; methylDackel_extract }                                                    from './modules/methylation'
include { gc_bias; idx_stats; flag_stats; fast_qc; insert_size_metrics; picard_metrics; tasmanian }     from './modules/compute_statistics.nf'
include { aggregate_emseq }                                                                             from './modules/aggregation.nf'


println "Processing " + params.flowcell + "... => " + outputDir
println "Cmd line: $workflow.commandLine"


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
        fastqc       = fast_qc( markDup.md_bams )
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