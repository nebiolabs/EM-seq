nextflow.enable.dsl=2

/* --------------- *
 * INPUT ARGUMENTS *
 * --------------- */
params.email                     = 'undefined'
params.flowcell                  = 'undefined'
params.path_to_genome_fasta      = 'undefined'   // path to the genome FASTA file, e.g. /path/to/genome.fa 
params.input_glob                = '*_R1.fastq*' // either the .bam or fastq read 1
params.project                   = 'project_undefined'
params.workflow                  = 'EM-seq'
params.outputDir                 = "em-seq_output" 
params.tmp_dir                   = '/tmp'
params.min_mapq                  = 20 // for methylation assessment.
params.max_input_reads           = "all_reads" // default is not downsampling , set to a number to downsample e.g. 1000000 is 500k read pairs
params.downsample_seed           = 42
params.enable_neb_agg            = 'True'
params.target_bed                = 'undefined' // BED file to intersect with methylKit output

include { alignReads; mergeAndMarkDuplicates; bwa_index; enough_reads; send_email; touchFile; samtools_faidx }     from './modules/alignment'
include { methylDackel_mbias; methylDackel_extract; convert_methylkit_to_bed; prepare_target_bed; intersect_beds; process_intersections; concatenate_intersections } from './modules/methylation'
include { gc_bias; idx_stats; flag_stats; fastqc; insert_size_metrics; picard_metrics; tasmanian }      from './modules/compute_statistics'
include { aggregate_emseq; multiqc }                                                                    from './modules/aggregation'


// identify and replace common R1/R2 naming patterns and return the read2 file name
// e.g. _R1.fastq -> _R2.fastq, _1.fastq -> _2.fastq, .R1. -> .R2., etc.
def replaceReadNumber(inputString) {
    return inputString.replace('_R1.', '_R2.')
                      .replace('_1.fastq', '_2.fastq')
                      .replace('_R1_', '_R2_')
                      .replace('.R1.', '.R2.')
                      .replace('_1.', '_2.')
                      .replace('.1.fastq', '.2.fastq')
}

// detect bam or fastq (or fastq.gz)
def detectFileType(file) {
    def file_str = file.toString()
    if (file_str.endsWith('.bam')) {
        return 'bam'
    } else if (file_str.endsWith('.fastq.gz') || file_str.endsWith('.fastq')) {
        // read2 exists for paired-end FASTQ?
        def read2File = replaceReadNumber(file_str)
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


workflow create_placeholder { touchFile("${workflow.workDir}/placeholder.r2.fastq") }

workflow {
    create_placeholder()
    main:
        // placeholder for R2 file, can't be a random file as that breaks nextflow's caching features
        // create the FILE here so it actually exists (touch)
        placeholder_r2 = file("${workflow.workDir}/placeholder.r2.fastq")
        
        // if reference is not indexed, index it.
        if (!file(params.path_to_genome_fasta).exists()) {
            println "Workflow failed: Genome file does not exist."
            System.exit(1)  // Exit with a custom status code
        }

        genome_index_ch = bwa_index()
        genome_fai_ch = samtools_faidx(Channel.fromPath(params.path_to_genome_fasta))

        reads = Channel
          .fromPath(params.input_glob)
          .filter { !(it.name ==~ /^.*2\.fastq.*$/) }
          .map { input_file ->
            def fileType = detectFileType(input_file)
            def read1File = input_file
            def read2File = placeholder_r2.toString()
            if (fileType == 'fastq_paired_end') {
                read2File = replaceReadNumber(input_file.toString())
           }
            if (read1File.toString() == read2File) {
                log.error("Error: Detected paired-end file with read1: ${read1File} but no read2. What is different in the file name?")
                throw new IllegalStateException("Invalid paired-end file configuration")
            }
	        def library = read1File.baseName.replaceFirst(/.fastq|.fastq.gz|.bam/,"").replaceFirst(/_1|\.1|.R1/,"")
            return [params.email, library, read1File, read2File, fileType]
          }
          //.join(genome_index_ch)
        

        reads.view()


        println "Processing " + params.flowcell + "... => " + params.outputDir
        println "Cmd line: $workflow.commandLine"


        // files with few reads will be filtered out and user will get an email.
        checking_reads = enough_reads(reads)
        passed_reads = checking_reads
                       .filter { tuple -> tuple[5].text.contains('pass') }
                       .map { tuple -> tuple[0..4] }
        failed_reads = checking_reads.filter { tuple -> tuple[5].text.contains('fail') }.map { tuple -> tuple[5].text }
        send_email( failed_reads.collect() )
        

        // align and mark duplicates
        alignedReads = alignReads( passed_reads, genome_index_ch )
        markDup      = mergeAndMarkDuplicates( alignedReads.bam_files )
        extract      = methylDackel_extract( markDup.md_bams, genome_index_ch )
        mbias        = methylDackel_mbias( markDup.md_bams, genome_index_ch )

        // intersect methylKit files with target BED file if provided
        if (params.target_bed != 'undefined') {
            target_bed_ch = Channel.fromPath(params.target_bed)
            
            methylkit_beds = convert_methylkit_to_bed( extract.extract_output.combine(genome_fai_ch) )
            
            prepared_bed = prepare_target_bed( target_bed_ch, genome_fai_ch )
            
            intersections = intersect_beds( methylkit_beds.methylkit_bed, prepared_bed.prepared_bed, genome_fai_ch )
            
            intersection_results = process_intersections( intersections.intersections )
            
            combined_results = concatenate_intersections( 
                intersection_results.intersection_results.collect(),
                intersection_results.intersection_summary.collect()
            )
        }

        // collect statistics
        gcbias       = gc_bias( markDup.md_bams, genome_index_ch )
        idxstats     = idx_stats( markDup.md_bams )
        flagstats    = flag_stats( markDup.md_bams )
        fastqc       = fastqc( markDup.md_bams )
        insertsize   = insert_size_metrics( markDup.md_bams ) 
        metrics      = picard_metrics( markDup.md_bams, genome_index_ch )
        mismatches   = tasmanian( markDup.md_bams, genome_index_ch )

        // Channels and processes that summarize all results

        // channel for internal summaries
        grouped_email_library = reads
	      .join( alignedReads.for_agg.groupTuple(by: [0, 1]), by: [0,1])
            .join( markDup.for_agg.groupTuple(by: [0,1]), by: [0,1] )
            .join( gcbias.for_agg.groupTuple(by: [0,1]), by: [0,1] )
            .join( idxstats.for_agg.groupTuple(by: [0,1]), by: [0,1] )
            .join( flagstats.for_agg.groupTuple(by: [0,1]), by: [0,1] )
            .join( fastqc.for_agg.groupTuple(by: [0,1]), by: [0,1] )
            .join( insertsize.for_agg.groupTuple(by: [0,1]), by: [0,1] )
            .join( mismatches.for_agg.groupTuple(by: [0,1]), by: [0,1] )
            .join( mbias.for_agg.groupTuple(by: [0,1]), by: [0,1] )
            .join( metrics.for_agg.groupTuple(by: [0,1]), by: [0,1] )

        if (params.enable_neb_agg.toString().toUpperCase() == "TRUE") {
            aggregate_emseq( grouped_email_library ) 
        }
       
        // channel for multiqc analysis
        all_results = grouped_email_library
         .join(insertsize.high_mapq_insert_size_metrics.groupTuple(by: [0, 1]), by: [0, 1])
         .map { items -> [items[0], items[7..-1]] }
         .groupTuple()
         .flatten()
         .toList()
         .map { items -> [items[0], items[7..-1]] }

        multiqc( all_results )
}
