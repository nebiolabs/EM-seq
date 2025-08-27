nextflow.enable.dsl=2

include { fastp }                          from './modules/fastp'
include { alignReads }                     from './modules/align_reads'
include { mergeAndMarkDuplicates }         from './modules/merge_and_mark_duplicates'
include { methylDackel_mbias }             from './modules/methyldackel_mbias'
include { methylDackel_extract }           from './modules/methyldackel_extract'
include { convert_methylkit_to_bed }       from './modules/convert_methylkit_to_bed'
include { prepare_target_bed }             from './modules/prepare_target_bed'
include { intersect_bed_with_methylkit }   from './modules/intersect_bed_with_methylkit'
include { group_bed_intersections }        from './modules/group_bed_intersections'
include { concatenate_intersections }      from './modules/concatenate_intersections'
include { gc_bias }                        from './modules/gc_bias'
include { idx_stats }                      from './modules/idx_stats'
include { flag_stats }                     from './modules/flag_stats'
include { fastqc }                         from './modules/fastqc'
include { insert_size_metrics }            from './modules/insert_size_metrics'
include { picard_metrics }                 from './modules/picard_metrics'
include { tasmanian }                      from './modules/tasmanian'
include { aggregate_emseq }                from './modules/aggregate_emseq.nf'
include { multiqc }                        from './modules/multiqc.nf'

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

bams = Channel.fromPath(params.ubam_dir + '/*.bam', checkIfExists: true)
       .map{it -> tuple(it.baseName, it)}

genome = params.genome
params.reference_list = params.genomes[genome]
genome_fa = Channel.value(params.reference_list.genome_fa)
genome_fai = Channel.value("${params.reference_list.genome_fa}.fai")

def checkFileSize (path) {
    return path.toFile().length() >= 200   // Minimum size in bytes for a read file to be considered valid
}


workflow {
    main:

        passed_bams = bams.filter { library, bam -> checkFileSize(bam) }
        failed_bams = bams.filter { library, bam -> !checkFileSize(bam) }
        failed_library_names = failed_bams.map { library, bam -> library }

        ////////// Filter libraries with insufficient reads //////////
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

        ///////// Trim, align and mark duplicates //////////
        fastp( passed_bams )
        fastq_chunks = fastp.out.trimmed_fastq
            .map { library, fq1_files, fq2_files -> 
                [fq1_files, fq2_files].transpose().collect { fq1, fq2 -> 
                    [library, fq1.baseName.split(".1.trimmed")[0], fq1, fq2] 
                }
            }
            .flatten()
            .collate(4)
        
        alignReads( fastq_chunks, params.reference_list.bwa_index )
        mergeAndMarkDuplicates( alignReads.out.bam_files.groupTuple() )
        md_bams = mergeAndMarkDuplicates.out.md_bams

        ///////// Methylation Calling //////////
        methylDackel_extract( md_bams, genome_fa, genome_fai )
        methylDackel_mbias( md_bams, genome_fa, genome_fai )

        //////// Intersect methylKit files with target BED file ////////
        if (params.reference_list.target_bed && !params.skip_target_bed) {

            convert_methylkit_to_bed( methylDackel_extract.out, genome_fa, genome_fai )
            prepare_target_bed( params.reference_list.target_bed, params.target_bed_slop, genome_fa, genome_fai )
            intersect_bed_with_methylkit(
                convert_methylkit_to_bed.out.methylkit_bed,
                prepare_target_bed.out,
                genome_fa,
                genome_fai
            )

            group_bed_intersections( intersect_bed_with_methylkit.out )

            concatenate_intersections(
                group_bed_intersections.out.results.collect(),
                group_bed_intersections.out.summary.collect()
            )
        }
        
        ///////// Collect statistics ///////
        gc_bias(  md_bams, genome_fa, genome_fai )
        idx_stats(  md_bams )
        flag_stats( md_bams )
        fastqc( md_bams )
        insert_size_metrics( md_bams ) 
        picard_metrics( md_bams, genome_fa, genome_fai )
        tasmanian( md_bams, genome_fa, genome_fai )

        //////// Collect files for internal summaries //////////
        grouped_library_results = bams
	        .join( alignReads.out.bam_files )
            .join( fastp.out.fastp_json )
            .join( mergeAndMarkDuplicates.out.log )
            .join( picard_metrics.out.for_agg )
            .join( gc_bias.out.for_agg )
            .join( idx_stats.out.for_agg )
            .join( flag_stats.out.for_agg )
            .join( fastqc.out.for_agg )
            .join( alignReads.out.nonconverted_counts )
            .join( tasmanian.out.for_agg )
            .join( methylDackel_mbias.out.for_agg )
            
        if (params.enable_neb_agg) {
            aggregate_emseq( grouped_library_results
                                .join( insert_size_metrics.out.for_agg )
            )
        }
       
        ////////// MultiQC analysis ///////////
        all_results = grouped_library_results
         .join(insert_size_metrics.out.high_mapq_insert_size_metrics)
         .map{[it[2..-1]]}
         .flatten()
         .collect()

        multiqc( all_results )
}
