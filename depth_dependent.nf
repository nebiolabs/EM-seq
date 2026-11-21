// Run the depth-dependent analyses on libraries that are already aligned.
//
// main.nf runs the same subworkflows after methylation calling. Use this entry point on
// aligned, duplicate-marked BAMs from earlier runs or other pipelines:
//
//   nextflow run depth_dependent.nf --genome T2T_chm13v2.0+meth_controls --bam_files 'bams/*.bam'
//
// The BAMs go through the depth check first. Once enough libraries pass, MethylDackel extracts CpG
// methylKit files from them and the analyses run. Give --mk_files as well to reuse methylKit
// files from an earlier run instead of extracting again. --mk_files alone runs the analyses
// without the depth check, since the check needs the BAMs.

if (!nextflow.version.matches('>=25.04')) {
    error "This pipeline requires Nextflow >=25.04. Found: ${nextflow.version}."
}

include { check_library_depth; depth_dependent_analysis } from './subworkflows/depth_dependent_analysis.nf'
include { methylDackel_extract }               from './modules/methyldackel_extract'
include { hasContigGroups }                    from './lib/contig_groups.nf'

// Glob of aligned, duplicate-marked, indexed BAMs (<name>.bam.bai or <name>.bai). Library
// names come from the BAM name with trailing lowercase dot-tags removed
// (1ng_x_1.ds.md.bam -> 1ng_x_1).
params.bam_files = null
// Glob of <library>_CpG.methylKit.gz files. Library names must match the BAM names.
params.mk_files = null
// Optional output directory of feature_cov_meth.nf for the same libraries. Its feature
// coverage and depth plots are added to the report.
params.feature_cov_dir = null

workflow {
    if (!params.genome || !params.genomes.containsKey(params.genome)) {
        error "Set --genome to one of: ${params.genomes.keySet().join(', ')}"
    }
    if (!params.bam_files && !params.mk_files) {
        error "Set --bam_files (aligned BAMs), --mk_files (CpG methylKit files), or both."
    }
    def reference_list = params.genomes[params.genome]
    def contig_groups = hasContigGroups(reference_list) ? file("${reference_list.gc_groups_dir}/contig_groups.tsv") : []
    def genome_fa = channel.value(reference_list.genome_fa)

    bams = params.bam_files
        ? channel.fromPath(params.bam_files, checkIfExists: true).map { bam ->
              def candidates = [file("${bam}.bai"), file(bam.toString().replaceFirst(/\.bam$/, '.bai'))]
              def bai = candidates.find { f -> f.exists() }
              if (!bai) {
                  error "No index for ${bam}: expected ${candidates.join(' or ')}"
              }
              tuple(bam.name.replaceFirst(/(\.[a-z]+)*\.bam$/, ''), bam, bai)
          }
        : channel.empty()

    if (params.bam_files) {
        check_library_depth( bams, contig_groups )
        ready = check_library_depth.out.ready
        depth_check = check_library_depth.out.table
    } else {
        ready = channel.value(true)
        depth_check = []
    }

    if (params.mk_files) {
        methylkits = channel.fromPath(params.mk_files, checkIfExists: true)
            .map { f -> tuple(f.name.replaceFirst(/_CpG\.methylKit(\.gz)?$/, ''), f) }
    } else {
        // Extract only once enough libraries pass the depth check, so shallow sets cost one mosdepth pass each.
        ready_bams = bams.combine(ready.filter { is_ready -> is_ready })
            .map { library, bam, bai, _ready -> tuple(library, bam, bai) }
        methylDackel_extract( ready_bams, genome_fa, channel.value(reference_list.genome_fai), 'cpg' )
        methylkits = methylDackel_extract.out.cpg
    }

    depth_dependent_analysis(
        methylkits,
        genome_fa,
        contig_groups,
        params.feature_cov_dir
            ? ['combined_feature_counts.tsv', 'feature_cov_meth_summary.tsv'].collect { f ->
                  file("${params.feature_cov_dir}/${f}", checkIfExists: true) }
            : [],
        bams,
        ready,
        depth_check
    )
}
