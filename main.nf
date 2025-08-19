nextflow.enable.dsl=2

include { alignReads; mergeAndMarkDuplicates }                                               from './modules/alignment'
include { methylDackel_mbias; methylDackel_extract; convert_methylkit_to_bed }                          from './modules/methylation'
include { prepare_target_bed; intersect_bed_with_methylkit;
          group_bed_intersections; concatenate_intersections }                                          from './modules/bed_processing'
include { gc_bias; idx_stats; flag_stats; fastqc; insert_size_metrics; picard_metrics; tasmanian }      from './modules/compute_statistics'
include { aggregate_emseq; multiqc }                                                                    from './modules/aggregation'

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}
genome = params.genome
params.reference_list = params.genomes[genome]
bams = Channel.fromPath(params.ubam_dir + '/*.bam', checkIfExists: true)
       .map{it -> tuple(it.baseName, it)}
genome_ch = Channel.fromPath(params.reference_list.genome_fa, checkIfExists: true).map{ fa -> tuple(fa, fa+'.fai') }
target_bed_ch = params.reference_list.target_bed ? Channel.fromPath(params.reference_list.target_bed, checkIfExists: true) : Channel.empty()

def checkFileSize (path) {
    return path.toFile().length() >= 200   // Minimum size in bytes for a read file to be considered valid
}


workflow {
    main:

        passed_bams = bams.filter { library, bam -> checkFileSize(bam) }
        failed_bams = bams.filter { library, bam -> !checkFileSize(bam) }
        failed_library_names = failed_bams.map { library, bam -> library }

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

        // align and mark duplicates
        alignedReads = alignReads( passed_bams, params.reference_list.bwa_index )
        markDup      = mergeAndMarkDuplicates( alignedReads.bam_files )
        extract      = methylDackel_extract( markDup.md_bams, genome_ch )
        mbias        = methylDackel_mbias( markDup.md_bams, genome_ch )
        
        // intersect methylKit files with target BED file if provided //
        if (target_bed_ch) {

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
        tasmanian    = tasmanian( markDup.md_bams, genome_ch )


        // channel for internal summaries
        grouped_library_results = bams
	        .join( alignedReads.bam_files )
            .join( alignedReads.nonconverted_counts )
            .join( alignedReads.fastp_reports )
            .join( markDup.log )
            .join( gcbias.for_agg )
            .join( idxstats.for_agg )
            .join( flagstats.for_agg )
            .join( fastqc.for_agg )
            .join( tasmanian.for_agg )
            .join( mbias.for_agg )
            .join( metrics.for_agg )

        if (params.enable_neb_agg.toString().toUpperCase() == "TRUE") {
            aggregate_emseq( grouped_library_results
                                .join( insertsize.for_agg )
            )
        }
       
        // channel for multiqc analysis
        all_results = grouped_library_results
         .join(insertsize.high_mapq_insert_size_metrics)
         .map{[it[2..-1]]}
         .flatten()
         .collect()

        multiqc( all_results )
}
