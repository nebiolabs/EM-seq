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
params.enable_neb_agg            = 'False'
params.target_bed                = 'undefined' // BED file to intersect with methylKit output
params.testing_mode              = false 

include { alignReads; mergeAndMarkDuplicates; genome_index; send_email; touchFile }                     from './modules/alignment'
include { methylDackel_mbias; methylDackel_extract; convert_methylkit_to_bed }                          from './modules/methylation'
include { prepare_target_bed; intersect_bed_with_methylkit;
          group_bed_intersections; concatenate_intersections }                                          from './modules/bed_processing'
include { gc_bias; idx_stats; flag_stats; fastqc; insert_size_metrics; picard_metrics; tasmanian }      from './modules/compute_statistics'
include { aggregate_emseq; multiqc }                                                                    from './modules/aggregation'
include { test_flagstats; test_alignment_metrics }                                                      from './modules/tests'    

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


def checkFileSize (path) {
    return path.toFile().length() >= 200   // Minimum size in bytes for a read file to be considered valid
}


workflow {
    main:
        placeholder_r2 = touchFile( "placeholder.r2.fastq" )

        // if reference is not indexed, index it.
        if (!file(params.path_to_genome_fasta).exists()) {
            println "Workflow failed: Genome file does not exist."
            System.exit(1)  // Exit with a custom status code
        }

        genome_indices = genome_index()
        bwa_index_ch = genome_indices.aligner_files
        genome_ch = genome_indices.genome_index

        reads = Channel
          .fromPath(params.input_glob)
          .filter { !(it.name ==~ /^.*2\.fastq.*$/) }
          .map { input_file ->
            def fileType = detectFileType(input_file)
            def read1File = input_file
            def read2File = placeholder_r2.val  // val blocks until the value is available
            if (fileType == 'fastq_paired_end') {
                read2File = replaceReadNumber(input_file.toString())
           }
            if (read1File.toString() == read2File.toString()) {
                log.error("Error: Detected paired-end file with read1: ${read1File} but no read2. What is different in the file name?")
                throw new IllegalStateException("Invalid paired-end file configuration")
            }
	        def library = read1File.baseName.replaceFirst(/.fastq|.fastq.gz|.bam/,"").replaceFirst(/_1|\.1|.R1/,"")
            return [library, read1File, read2File, fileType]
          }

        println "Processing " + params.flowcell + "... => " + params.outputDir
        println "Cmd line: $workflow.commandLine"

        passed_reads = reads.filter { library, read1File, read2File, fileType -> checkFileSize(read1File) }
        failed_reads = reads.filter { library, read1File, read2File, fileType -> !checkFileSize(read1File) }
        failed_library_names = failed_reads.map { library, read1File, read2File, fileType -> library }

        // Send email if there are failed libraries
        failed_library_names.collect().subscribe { names ->
            if (names.size() > 0) {
                def joined_names = names.join('<br>')
                sendMail {
                    to params.email
                    subject 'File Read Check'
                    body """
                    <html>
                      <body>
                        <p>The following libraries:<br> <strong>${joined_names}</strong> do not have enough reads. <br> Continuing with remaining libraries. </p>
                      </body>
                    </html>
                    """
                }
            }
        }

        alignedReads = alignReads( passed_reads, bwa_index_ch )
        markDup      = mergeAndMarkDuplicates( alignedReads.aligned_bams )

        extract      = methylDackel_extract( markDup.md_bams, genome_ch )
        mbias        = methylDackel_mbias( markDup.md_bams, genome_ch )

        // intersect methylKit files with target BED file if provided //
        if (params.target_bed != 'undefined') {
            target_bed_ch = Channel.fromPath(params.target_bed)

            methylkit_beds = convert_methylkit_to_bed( extract.extract_output.combine(genome_ch) )
            prepared_bed = prepare_target_bed( target_bed_ch, genome_ch )
            intersections = intersect_bed_with_methylkit(
                methylkit_beds.methylkit_bed,
                prepared_bed.prepared_bed.first(),
                genome_ch
            )

            intersection_results = group_bed_intersections( intersections.intersections )

            combined_results = concatenate_intersections(
                intersection_results.intersection_results.collect(),
                intersection_results.intersection_summary.collect()
            )
        }

        // collect statistics
        gcbias       = gc_bias( markDup.md_bams, genome_ch )
        idxstats     = idx_stats( markDup.md_bams )
        flagstats    = flag_stats( markDup.md_bams )
        fastqc       = fastqc( markDup.md_bams )
        insertsize   = insert_size_metrics( markDup.md_bams )
        metrics      = picard_metrics( markDup.md_bams, genome_ch )
        mismatches   = tasmanian( markDup.md_bams, genome_ch )


        // channel for internal summaries
        grouped_library_results = markDup.md_bams
            .join( alignedReads.metadata )
	        .join( alignedReads.fastp_reports )
            .join( alignedReads.nonconverted_counts )
            .join( markDup.for_agg )
            .join( gcbias.for_agg )
            .join( idxstats.for_agg )
            .join( flagstats.for_agg )
            .join( fastqc.for_agg )
            .join( mismatches.for_agg )
            .join( mbias.for_agg )
            .join( metrics.for_agg )

        if (params.enable_neb_agg.toString().toUpperCase() == "TRUE") {
            aggregate_emseq( grouped_library_results
                                .join( insertsize.for_agg )
                                .join( passed_reads.map { library, read1File, _read2File, _fileType -> [library, read1File] })
            )
        }

        // channel for multiqc analysis -
        all_results = grouped_library_results
            .join(insertsize.high_mapq_insert_size_metrics)
            .map { it[1,2] + it[6..-1].flatten() } //multiqc needs all the files (without the library name)

        multiqc( all_results )

        // ONLY for testing.
        if (params.run_tests.toString().toUpperCase() == "TRUE") {
            test_flagstats( flagstats.for_agg  )
            test_alignment_metrics( metrics.for_agg )
        }
        else
        {
            println("${params.run_tests.toString().toUpperCase()}")
        }
}
