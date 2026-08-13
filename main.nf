// Registered before any validation so ALL errors are caught
include { registerEmailNotifications } from './lib/notifications.nf'
registerEmailNotifications()

include { createVersionsFile }                                from './lib/versions.nf'
include { format_ngs_agg_opts }                               from './modules/aggregate_results'
include { fastp }                                             from './modules/fastp'
include { mergeFastpJson }                                    from './modules/merge_fastp_json'
include { alignReads }                                        from './modules/align_reads'
include { mergeAndMarkDuplicates }                            from './modules/merge_and_mark_duplicates'
include { methylDackel_mbias }                                from './modules/methyldackel_mbias'
include { methylDackel_extract }                              from './modules/methyldackel_extract'
include { extract_cytosine_report }                           from './modules/extract_cytosine_report'
include { combine_nonconverted_counts }                       from './modules/combine_nonconverted_counts'
include { convert_methylkit_to_bed }                          from './modules/convert_methylkit_to_bed'
include { prepare_target_bed }                                from './modules/prepare_target_bed'
include { intersect_bed_with_methylkit }                      from './modules/intersect_bed_with_methylkit'
include { group_bed_intersections }                           from './modules/group_bed_intersections'
include { concatenate_files as concat_intersections;
          concatenate_files as concat_positional_summaries;
          concatenate_files as concat_region_summaries;}      from './modules/concatenate_files'
include { gc_bias }                                           from './modules/gc_bias'
include { idx_stats }                                         from './modules/idx_stats'
include { flagstats }                                        from './modules/flagstats'
include { fastqc }                                            from './modules/fastqc'
include { insert_size_metrics }                               from './modules/insert_size_metrics'
include { picard_metrics }                                    from './modules/picard_metrics'
include { tasmanian }                                         from './modules/tasmanian'
include { find_switchback_reads }                             from './modules/find_switchback_reads'
include { aggregate_results }                                 from './modules/aggregate_results.nf'
include { multiqc }                                           from './modules/multiqc.nf'

def checkFileSize (path) {
    return path.toFile().length() >= 200   // Minimum size in bytes for a read file to be considered valid
}

workflow {
    main:

        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
          exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
        }

         bams = channel.fromPath(params.ubam_dir + '/*.bam', checkIfExists: true)
               .map{it -> tuple(it.baseName, it)}

         def reference_list = params.genomes[params.genome]
         genome_fa = channel.value(reference_list.genome_fa)
         genome_fai = channel.value(reference_list.genome_fai)

        // Validate reference indices exist before running
        def ref = reference_list.bwa_index
        if (!file("${ref}.bwameth.c2t.amb").exists()) {
            exit 1, "bwameth index not found for ${ref}. Run: bwameth.py index ${ref}"
        }
        if (!file("${reference_list.genome_fa}.fai").exists()) {
            exit 1, "Fasta index (.fai) not found for ${reference_list.genome_fa}. Run: samtools faidx ${reference_list.genome_fa}"
        }

        passed_bams = bams.filter { _library, bam -> checkFileSize(bam) }
        failed_bams = bams.filter { _library, bam -> !checkFileSize(bam) }
        failed_library_names = failed_bams.map { library, _bam -> library }

        ////////// Filter libraries with insufficient reads //////////
        failed_library_names.collect().subscribe { names ->
            if (names.size() > 0) {
                def joined_names = names.join('<br>')
                sendMail(
                    to: params.email,
                    subject: 'File Read Check',
                    body: """
                    <html>
                      <body>
                        <p>The following libraries:<br> <strong>${joined_names}</strong> do not have enough reads. <br> Continuing with remaining libraries. </p>
                      </body>
                    </html>
                    """
                )
            }
        }

        ///////// Trim, align and mark duplicates //////////
        // Split each uBAM into byte-range chunks so trimming runs in parallel per chunk.
        chunk_size = params.bamslice_chunk_size as long
        bam_chunks = passed_bams.flatMap { library, bam ->
            def file_size = bam.size()
            def num_chunks = ((file_size + chunk_size - 1) / chunk_size) as long
            (0..<num_chunks).collect { i ->
                def start = i * chunk_size
                def end = Math.min(start + chunk_size, file_size)
                tuple(library, start, end)
            }
        }

        fastp( passed_bams.combine(bam_chunks, by:0) )
        mergeFastpJson( fastp.out.fastp_json.groupTuple() )

        fastq_chunks = fastp.out.trimmed_fastq.map { library, chunk_name, fq_files ->
            def fq_list = fq_files instanceof List ? fq_files : [fq_files]
            if (params.single_end) {
                tuple(library, chunk_name, fq_list[0])
            }
            else {
                def r1 = fq_list.find { fq -> fq.name.contains('.1.trimmed.fastq') }
                def r2 = fq_list.find { fq -> fq.name.contains('.2.trimmed.fastq') }
                tuple(library, chunk_name, [r1, r2])
            }
        }
        alignReads( passed_bams.combine(fastq_chunks, by:0), reference_list.bwa_index )
        mergeAndMarkDuplicates( alignReads.out.bam_files.groupTuple() )
        md_bams = mergeAndMarkDuplicates.out.md_bams

        ///////// Methylation Calling //////////
        methylDackel_extract( md_bams, genome_fa, genome_fai )
        methylDackel_mbias( md_bams, genome_fa, genome_fai )

        //////// Intersect methylKit files with target BED file ////////
        if (reference_list.target_bed && !params.skip_target_bed) {
            extract_cytosine_report( md_bams, genome_fa, genome_fai )
            convert_methylkit_to_bed( extract_cytosine_report.out.report, genome_fa, genome_fai )
            prepare_target_bed( reference_list.target_bed, params.target_bed_slop, genome_fa, genome_fai )
            intersect_bed_with_methylkit(
                convert_methylkit_to_bed.out.methylkit_bed,
                prepare_target_bed.out.bed,
                genome_fa,
                genome_fai
            )

            group_bed_intersections( intersect_bed_with_methylkit.out.tsv )

            concat_intersections(
                group_bed_intersections.out.intersections.collect(),
                "intersections"
            )
            concat_positional_summaries(
                group_bed_intersections.out.positional_summary.collect(),
                "positional_summaries"
            )
            concat_region_summaries(
                intersect_bed_with_methylkit.out.region_summary.collect(),
                "region_summaries"
            )
            
        }
        
        ///////// Collect statistics ///////
        gc_bias(  md_bams, genome_fa, genome_fai )
        idx_stats(  md_bams )
        flagstats( md_bams )
        fastqc( md_bams )
        picard_metrics( md_bams, genome_fa, genome_fai )
        tasmanian( md_bams, genome_fa, genome_fai )
        combine_nonconverted_counts( alignReads.out.nonconverted_counts.groupTuple() )
        find_switchback_reads( md_bams, genome_fa )

        //////// Collect files for internal summaries //////////
        agg_opts = [
        ['--workflow', "${params.workflow}"],
        ['--workflow_version', 'unspecified'],
        ['--bam', mergeAndMarkDuplicates.out.md_bams.map{ row -> tuple(row[0], row[1]) }],
        ['--bai', mergeAndMarkDuplicates.out.md_bams.map{ row -> tuple(row[0], row[2]) }],
        ['--metadata_bam_file', bams],
        ['--fastp', mergeFastpJson.out.merged_json],
        ['--aln', picard_metrics.out.for_agg ],
        ['--gc', gc_bias.out.for_agg ],
        ['--dup', mergeAndMarkDuplicates.out.log],
        ['--idx_stats', idx_stats.out.for_agg],
        ['--flagstat', flagstats.out.for_agg],
        ['--fastqc', fastqc.out.for_agg],
        ['--nonconverted_read_counts', combine_nonconverted_counts.out.for_agg ],
        ['--tasmanian', tasmanian.out.for_agg ],
        ['--combined_mbias_records', methylDackel_mbias.out.for_agg ],
        ['--fgbio_switchbacks', find_switchback_reads.out.for_agg ]
        ]
        
        multiqc_opts = agg_opts.clone()
        if (!params.single_end) {
            insert_size_metrics( md_bams )
            multiqc_opts << ['--insert', insert_size_metrics.out.high_mapq]
            agg_opts << ['--insert', insert_size_metrics.out.for_agg]
        }

        if (params.workflow_name_modifier) {
            agg_opts << ['--workflow_name_modifier', params.workflow_name_modifier]
        }

        if (params.enable_neb_agg) {
            agg_tuple = format_ngs_agg_opts(agg_opts)
            aggregate_results( agg_tuple )
        }

        ////////// MultiQC analysis ///////////
        multiqc_tuple = format_ngs_agg_opts(multiqc_opts)
        versions_file = createVersionsFile(channel.topic('versions'))
        multiqc_files = multiqc_tuple.map{ row -> row[2] }.flatten().concat(versions_file).collect()

        multiqc( multiqc_files )
}
