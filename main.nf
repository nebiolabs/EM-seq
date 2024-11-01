nextflow.enable.dsl=2

/* ------------------------ *
 * INCLUDE IN CONFIG FILE!! *
 * ------------------------ */
params.default_dest_path = '/mnt/galaxy/tmp/users'
params.tmp_dir           =  '/tmp'
params.path_to_ngs_agg   = '/mnt/bioinfo/prg/ngs-aggregate_results/current'

/* --------------- *
 * INPUT ARGUMENTS *
 * --------------- */
params.email       = 'undefined'
params.flowcell    = 'undefined'
params.genome      = 'undefined'
params.input_glob  =  '*.{1,2}.fastq*'
params.project     = 'project_undefined'
params.workflow    = 'Automated EM-seq'
params.outputDir = "em-seq_output" // params.outdir ?: new File([default_dest_path, "email",flowcell].join(File.separator))
params.min_mapq = 20 // for methylation assessment.
params.max_input_reads = "all_reads" // default is not downsampling 
params.downsample_seed = 42

// include { PATH_TO_TILES_KNOWN } from './modules/path_to_tiles_provided'
include { alignReads; mergeAndMarkDuplicates }                                                          from './modules/alignment'
include { methylDackel_mbias; methylDackel_extract }                                                    from './modules/methylation'
include { gc_bias; idx_stats; flag_stats; fastqc; insert_size_metrics; picard_metrics; tasmanian }      from './modules/compute_statistics'
include { aggregate_emseq }                                                                             from './modules/aggregation'


println "Processing " + params.flowcell + "... => " + params.outputDir
println "Cmd line: $workflow.commandLine"


// detect bam or fastq (or fastq.gz)
def detectFileType(file) {
    file_str = file.toString()
    if (file_str.endsWith('.bam')) {
        return 'bam'
    } else if (file_str.endsWith('.fastq.gz') || file_str.endsWith('.fastq')) {
        // read2 exists for paired-end FASTQ?
        def read2File = file_str.replace('_R1.fastq', '_R2.fastq').replace('_1.fastq', '_2.fastq').replace("_R1_","_R2_").replace(".R1.",".R2.")
        if (new File(read2File).exists()) {
            return 'fastq_paired_end'
        } else {
            return 'fastq_single_end'
        }
    } else {
        println "Unknown file type for $file_str. If fastq, check if _R1_ or .R1. patterns are used."
        return 'unknown'
    }
}

Channel
    .fromPath(params.input_glob)
    .map { input_file ->
        def fileType = detectFileType(input_file)
        def read1File = input_file
        def read2File = params.genome // fake it to have the same number of elements in the tuple.
        if (fileType == 'fastq_paired_end') {
            read2File = input_file.toString().replace('_R1.', '_R2.').replace('_1.fastq', '_2.fastq').replace("_R1_","_R2_").replace(".R1.",".R2.").replace('_1.', '_2.')
        }
        def genome = params.genome
        return [read1File, read2File, genome, fileType]
    }.set{ inputChannel}



 workflow {
    main:
        // inputChannel.view()
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
        mismatches   = tasmanian( markDup.md_bams )

        // Channel for programs that summarize all results
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
