#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Require Nextflow version 24.04 or later for topic channels support
// Topic channels are used for version tracking in subprocesses
if( !nextflow.version.matches('>=24.04') ) {
    error "This workflow requires Nextflow version 24.04 or later. Please upgrade your Nextflow installation:\n" +
          "  Current version: ${nextflow.version}\n" +
          "  Required version: 24.04+\n" +
          "  Upgrade with: nextflow self-update\n" +
          "  Or install latest: curl -s https://get.nextflow.io | bash"
}

// Enable topic preview flag for versions < 25.04 (graduated from preview in 25.04)
if( nextflow.version.matches('>=24.04') && nextflow.version.matches('<25.04') ) {
    nextflow.preview.topic = true
}


include { download_epd_promoters;
          download_liftover_chain;
          crossmap_epd_promoters;
          download_cpg_islands;
          download_refseq_gtf;
          download_assembly_report;
          normalize_cpg_islands;
          normalize_refseq_features;
          download_dfam_annotations }                     from './modules/prepare_features'
include { calculate_feature_methylation;
          combine_methylation }                           from './modules/feature_methylation'
include { calculate_feature_counts;
          combine_counts }                                from './modules/feature_counts'

params.genome = null
params.methylation_bed = null
params.bam_files_glob = '*.bam'  // Indexes auto-detected: .bai, .csi, .crai
params.tmp_dir = '/tmp'
params.outputDir = 'output'

params.ucsc_cpg_islands_url = null
params.refseq_gtf_url = null
params.ncbi_assembly_report_url = null
params.dfam_url = null
params.epd_promoters_url = null

if (!params.genome) {
    exit 1, "Error: --genome parameter is required"
}
if (!params.bam_files_glob) {
    exit 1, "Error: --bam_files_glob parameter is required"
}

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

genome = params.genome
params.reference_list = params.genomes ? params.genomes[genome] : [:]
genome_fa = params.reference_list.genome_fa ?: params.genome_fa
genome_fai = params.reference_list.genome_fai ?: params.genome_fai

ucsc_cpg_islands_url = params.ucsc_cpg_islands_url ?: params.reference_list.ucsc_cpg_islands_url
refseq_gtf_url = params.refseq_gtf_url ?: params.reference_list.refseq_gtf_url
ncbi_assembly_report_url = params.ncbi_assembly_report_url ?: params.reference_list.ncbi_assembly_report_url
dfam_url = params.dfam_url ?: params.reference_list.dfam_url
epd_promoters_url = params.epd_promoters_url ?: params.reference_list.epd_promoters_url
epd_liftover_url = params.reference_list.epd_liftover_url ?: null

refseq_feature_types = Channel.of(
    'promoter',
    'transcriptional_cis_regulatory_region',
    'enhancer',
    'silencer',
    'mobile_genetic_element',
    'lnc_RNA',
    'snRNA',
    'snoRNA',
    'primary_transcript',
    'tRNA',
    'exon',
    'mRNAexon1',
    'mRNAexon'
)

workflow {
    main:
        bam_with_index = Channel.fromPath(params.bam_files_glob, checkIfExists: true)
            .map { alignment ->
                def base = alignment.toString().replaceAll(/\.(bam|cram)$/, '')
                def idx = ["${alignment}.bai", "${base}.bai",
                           "${alignment}.csi", "${base}.csi",
                           "${alignment}.crai", "${base}.crai"]
                    .collect { file(it) }
                    .find { it.exists() }
                if (!idx) error "No index found for: ${alignment}"
                tuple(alignment, idx)
            }
            .collect()

        if (epd_promoters_url) {
            download_epd_promoters(epd_promoters_url)

            if (epd_liftover_url) {
                download_liftover_chain(epd_liftover_url)
                crossmap_epd_promoters(
                    download_epd_promoters.out.gtf,
                    download_liftover_chain.out.chain_file
                )

                epd_gtf = crossmap_epd_promoters.out.gtf
                    .map { gtf -> tuple('epd_promoter', gtf) }
            } else {
                epd_gtf = download_epd_promoters.out.gtf
                    .map { gtf -> tuple('epd_promoter', gtf) }
            }
        } else {
            epd_gtf = Channel.empty()
        }

        if (ncbi_assembly_report_url) {
            download_assembly_report(ncbi_assembly_report_url)
            assembly_report = download_assembly_report.out.report
        } else {
            assembly_report = Channel.empty()
        }

        if (ucsc_cpg_islands_url && ncbi_assembly_report_url) {
            download_cpg_islands(ucsc_cpg_islands_url)
            normalize_cpg_islands(download_cpg_islands.out.ucsc_file, assembly_report)
            cpg_gtf = normalize_cpg_islands.out.gtf
                .map { gtf -> tuple('cpg_island', gtf) }
        } else {
            cpg_gtf = Channel.empty()
        }

        if (refseq_gtf_url && ncbi_assembly_report_url) {
            download_refseq_gtf(refseq_gtf_url)
            refseq_combined = download_refseq_gtf.out.gtf
                .combine(assembly_report)
                .combine(Channel.value(file(genome_fai, checkIfExists: true)))
                .combine(refseq_feature_types)

            normalize_refseq_features(refseq_combined)
            refseq_gtfs = normalize_refseq_features.out.gtf
        } else {
            refseq_gtfs = Channel.empty()
        }

        if (dfam_url) {
            download_dfam_annotations(dfam_url)
            dfam_gtf = download_dfam_annotations.out.gtf
                .map { gtf -> tuple('dfam_repeat', gtf) }
        } else {
            dfam_gtf = Channel.empty()
        }

        if (params.methylation_bed) {
            methylation_bed = Channel.fromPath(params.methylation_bed, checkIfExists: true)

            all_gtfs = epd_gtf
                .mix(cpg_gtf)
                .mix(refseq_gtfs)
                .mix(dfam_gtf)

            calculate_feature_methylation(all_gtfs, methylation_bed)

            combine_methylation(
                calculate_feature_methylation.out.methylation
                    .map { feature_type, tsv -> tsv }
                    .collect()
            )
        }

        all_features_for_counts = epd_gtf
            .map { type, file -> tuple(type, file, 'gene_id') }
            .mix(cpg_gtf.map { type, file -> tuple(type, file, 'transcript_id') })
            .mix(refseq_gtfs.map { type, file -> tuple(type, file, 'gene_id') })
            .mix(dfam_gtf.map { type, file -> tuple(type, file, 'gene_id') })

        calculate_feature_counts(
            all_features_for_counts.map { type, file, attr -> tuple(type, file) },
            bam_with_index,
            all_features_for_counts.map { type, file, attr -> attr }
        )

        combine_counts(
            calculate_feature_counts.out.counts
                .map { feature_type, tsv -> tsv }
                .collect()
        )
}
