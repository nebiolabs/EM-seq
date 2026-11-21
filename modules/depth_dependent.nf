// Depth-dependent analyses of CpG-context Cs, one per strand. See
// subworkflows/depth_dependent_analysis.nf for how these processes connect and how the depth
// check decides whether they run.

// Mean depth per contig, for the depth check. The filters match what MethylDackel counts:
// MAPQ >= --methyl_quality_threshold, and no duplicates, secondary, supplementary or QC-failed
// reads. mosdepth counts the overlapping part of a read pair once, as MethylDackel does.
process depth_dependent_contig_depth {
    label 'medium_cpu'
    tag "${library}"
    publishDir "${params.outputDir}/depth_dependent/contig_depth", mode: 'copy', pattern: '*.mosdepth.summary.txt'
    conda "bioconda::mosdepth=0.3.14"

    input:
        tuple val(library), path(bam), path(bai)

    output:
        path("${library}.mosdepth.summary.txt"), emit: summary
        tuple val("${task.process}"), val('mosdepth'), eval('mosdepth --version | sed "s/^mosdepth //"'), topic: versions

    script:
    """
    mosdepth --threads ${task.cpus} --no-per-base --mapq ${params.methyl_quality_threshold} \\
        --flag 3844 ${library} ${bam}
    """
}

// Per-library pass or fail against --depth_dependent_min_contig_depth. The workflow reads
// depth_dependent_check.tsv to decide whether to run the depth-dependent analyses.
process depth_dependent_check {
    label 'process_single'
    publishDir "${params.outputDir}/depth_dependent", mode: 'copy'
    conda "conda-forge::python=3.12"

    input:
        path(summaries)
        path(contig_groups)   // [] when the genome has no contig groups

    output:
        path('depth_dependent_check.tsv'), emit: table
        path('depth_dependent_contig_depth.tsv'), emit: contig_depth

    script:
    def groups_arg = contig_groups ? "--contig-groups ${contig_groups} --contig-types ${params.depth_dependent_check_contig_types}" : ''
    """
    : library_order.py  # imported module, named so edits to it invalidate the cache
    depth_dependent_check.py ${summaries} --min-contig-depth ${params.depth_dependent_min_contig_depth} ${groups_arg}
    """
}

process count_reference_cpgs {
    label 'process_single'
    publishDir "${params.outputDir}/depth_dependent", mode: 'copy'
    conda "conda-forge::python=3.12 conda-forge::polars=1.31.0 conda-forge::bokeh=3.7.3 conda-forge::numpy=2.5.4 conda-forge::scipy=1.18.1 conda-forge::pyarrow=26.0.1"

    input:
        val(genome_fa)
        path(contig_groups)   // [] when the genome has no contig groups

    output:
        path('reference_cpgs.tsv'), emit: sites

    script:
    def groups_arg = contig_groups ? "--contig-groups ${contig_groups}" : ''
    """
    count_reference_cpgs.py ${genome_fa} ${groups_arg} -o reference_cpgs.tsv
    """
}

process depth_dependent_prep {
    label 'depth_dependent'
    tag "${library}"
    publishDir "${params.outputDir}/depth_dependent/depth", mode: 'copy', pattern: '*.depth_hist.tsv'
    conda "conda-forge::python=3.12 conda-forge::polars=1.31.0 conda-forge::bokeh=3.7.3 conda-forge::numpy=2.5.4 conda-forge::scipy=1.18.1 conda-forge::pyarrow=26.0.1"

    input:
        tuple val(library), path(cpg_methylkit)
        path(reference_cpgs)

    output:
        tuple val(library), path("${library}.cs.parquet"), path("${library}.depth_hist.tsv"), emit: cpgs
        tuple val("${task.process}"), val('polars'), eval('python -c "import polars; print(polars.__version__)"'), topic: versions

    script:
    """
    depth_dependent_prep.py ${cpg_methylkit} ${reference_cpgs} --library ${library}
    """
}

// Read counts for the report's overview. Usable reads are the ones MethylDackel uses with
// this pipeline's settings: mapped primary reads that are not duplicates, QC fails or
// supplementary, with MAPQ >= --methyl_quality_threshold, and properly paired when paired
// (MethylDackel ignores discordant pairs and singletons by default).
process depth_dependent_read_stats {
    tag "${library}"
    cpus 4
    memory '4 GB'
    conda "bioconda::samtools=1.21"

    input:
        tuple val(library), path(bam), path(bai)

    output:
        path("${library}.read_stats.tsv"), emit: stats

    script:
    """
    samtools flagstat -@ ${task.cpus} -O tsv ${bam} > flagstat.tsv
    usable=\$(samtools view -c -@ ${task.cpus} -q ${params.methyl_quality_threshold} -F 0xF04 \\
        -e '!flag.paired || flag.proper_pair' ${bam})
    # The flowcell is the most common one among the first reads. It is matched by pattern, not
    # by field position. Illumina and Element names end in flowcell:lane:tile:x:y, and MGI
    # names start with the flowcell followed by L<lane>C<col>R<row>. Other names give none.
    # head closes the pipe early, so the pipeline's status is ignored.
    flowcell=\$( (samtools view ${bam} | head -n 10000 | cut -f1 |
        sed -nE -e 's/^.*:([^:]+)(:[0-9]+){4}\$/\\1/p;t' -e 's/^([A-Z]+[0-9]+)L[0-9]+C[0-9]+R[0-9]+.*\$/\\1/p' |
        sort | uniq -c | sort -rn | head -n 1 | awk '{print \$2}') || true )
    # flagstat columns: QC-passed, QC-failed, category
    awk -F'\\t' -v lib=${library} -v usable="\$usable" -v flowcell="\$flowcell" '
        \$3 == "primary"            { primary = \$1 + \$2 }
        \$3 == "primary mapped"     { mapped = \$1 + \$2 }
        \$3 == "primary duplicates" { dups = \$1 + \$2 }
        END {
            print "library\\tprimary_reads\\tprimary_mapped\\tprimary_duplicates\\tusable_reads\\tmin_mapq\\tflowcell"
            print lib "\\t" primary "\\t" mapped "\\t" dups "\\t" usable "\\t${params.methyl_quality_threshold}\\t" flowcell
        }' flagstat.tsv > ${library}.read_stats.tsv
    """
}

process depth_dependent_merge {
    label 'depth_dependent'
    publishDir "${params.outputDir}/depth_dependent", mode: 'copy'
    conda "conda-forge::python=3.12 conda-forge::polars=1.31.0 conda-forge::bokeh=3.7.3 conda-forge::numpy=2.5.4 conda-forge::scipy=1.18.1 conda-forge::pyarrow=26.0.1"

    input:
        path(parquets)
        path(reference_cpgs)

    output:
        path('depth_dependent_*.tsv'), emit: stats
        path('depth_dependent_matrix.parquet'), optional: true, emit: matrix

    script:
    def matrix_arg = params.depth_dependent_write_matrix ? '--matrix' : ''
    """
    : library_order.py feature_depth_bokeh.py  # imported modules, named so edits to them invalidate the cache
    depth_dependent_merge.py ${reference_cpgs} ${parquets} \\
        --min-depth ${params.depth_dependent_min_depth} ${matrix_arg} \\
        --report-depths ${params.depth_dependent_report_depths} \\
        --core-jaccard ${params.depth_dependent_core_jaccard} --core-depth ${params.depth_dependent_core_depth}
    """
}

process depth_dependent_report {
    // Embedding feature_cov_meth.nf's per-feature depth plot peaked at 6.4 GB on 8 human libraries.
    label 'process_high_memory'
    // The report sits at the top of the output next to the MultiQC report. Its tables go
    // with the other depth-dependent tables.
    publishDir "${params.outputDir}", mode: 'copy', pattern: '*.html'
    publishDir "${params.outputDir}/depth_dependent", mode: 'copy', pattern: '*.tsv'
    conda "conda-forge::python=3.12 conda-forge::polars=1.31.0 conda-forge::bokeh=3.7.3 conda-forge::numpy=2.5.4 conda-forge::scipy=1.18.1 conda-forge::pyarrow=26.0.1"

    input:
        path(depth_hists)
        path(merge_stats)   // [] when there is only one library
        path(feature_files) // [] or feature_cov_meth.nf outputs: combined_feature_counts.tsv, feature_cov_meth_summary.tsv
        path(read_stats)    // [] or <library>.read_stats.tsv files for the overview tab
        path(depth_check) // [] or depth_dependent_check.tsv

    output:
        path('*depth_dependent_report.html'), emit: html
        path('depth_dependent_depth_summary.tsv'), emit: depth_summary
        path('depth_dependent_clusters.tsv'), optional: true, emit: clusters
        tuple val("${task.process}"), val('bokeh'), eval('python -c "import bokeh; print(bokeh.__version__)"'), topic: versions

    script:
    def merge_arg = merge_stats ? '--merge-dir .' : ''
    def names = (feature_files instanceof List ? feature_files : [feature_files]).findAll { f -> f }*.name
    def feature_args = (names.contains('combined_feature_counts.tsv') ? '--feature-counts combined_feature_counts.tsv ' : '') +
                       (names.contains('feature_cov_meth_summary.tsv') ? '--feature-cov-summary feature_cov_meth_summary.tsv' : '')
    // Nextflow hashes the bin scripts a task's script text names. depth_dependent_report.py
    // imports these modules, so they are named here to make edits to them invalidate the cache.
    def stats_arg = read_stats ? "--read-stats ${read_stats}" : ''
    def check_arg = depth_check ? "--depth-check ${depth_check}" : ''
    """
    : feature_depth_bokeh.py feature_cov_meth_bokeh.py library_order.py report_overlay.py report_style.py depth_dependent_merge.py  # imported modules, named so edits to them invalidate the cache
    depth_dependent_report.py ${depth_hists} ${merge_arg} ${feature_args} ${stats_arg} ${check_arg} \\
        --min-depth ${params.depth_dependent_min_depth}
    """
}
