// Depth-dependent analyses of CpG-context Cs, run when enough libraries pass a depth check.
// Each C is one strand of a CpG, with the reads on that strand as its depth.
//
// check_library_depth measures the mean depth of every contig in each library with mosdepth. A
// library passes when every checked contig (autosomes by default) has a mean depth above
// params.depth_dependent_min_contig_depth. The analyses run when params.depth_dependent_min_libraries
// libraries pass. Otherwise nothing in depth_dependent_analysis runs.
//
// Once they run, every library is analysed, shallow ones included. Each gets a depth
// histogram against the reference C count, and two or more are merged into one per-C matrix.
// The merge produces pairwise correlation, upset patterns, library-specific Cs and clustering
// at each depth in params.depth_dependent_report_depths.

include { depth_dependent_contig_depth; depth_dependent_check; count_reference_cpgs; depth_dependent_prep;
          depth_dependent_merge; depth_dependent_read_stats; depth_dependent_report } from '../modules/depth_dependent'

workflow check_library_depth {
    take:
        bams            // tuple(library, bam, bai)
        contig_groups   // value: contig_groups.tsv, or [] to choose checked contigs by name

    main:
        depth_dependent_contig_depth( bams )
        depth_dependent_check( depth_dependent_contig_depth.out.summary.collect(), contig_groups )
        ready = depth_dependent_check.out.table.map { tsv -> depthDependentReady(tsv, params.depth_dependent_min_libraries as int) }

    emit:
        ready = ready                                   // value: true when the analyses should run
        table = depth_dependent_check.out.table         // depth_dependent_check.tsv
}

workflow depth_dependent_analysis {
    take:
        cpg_methylkits   // tuple(library, <library>_CpG.methylKit.gz)
        genome_fa        // value: reference FASTA path
        contig_groups    // value: contig_groups.tsv, or [] to fall back to name-based exclusion
        feature_files    // value: feature_cov_meth.nf outputs to add to the report, or []
        bams             // tuple(library, bam, bai) for the overview's read counts; may be empty
        ready            // value: check_library_depth.out.ready, or channel.value(true) to skip the check
        depth_check      // value: check_library_depth.out.table for the report, or []

    main:
        // Every input waits on ready, so no task starts here when too few libraries pass.
        is_ready = ready.filter { r -> r }
        ready_fa = is_ready.combine(genome_fa).map { _o, fa -> fa }.first()
        ready_methylkits = cpg_methylkits.combine(is_ready).map { library, mk, _o -> tuple(library, mk) }
        ready_bams = bams.combine(is_ready).map { library, bam, bai, _o -> tuple(library, bam, bai) }

        count_reference_cpgs( ready_fa, contig_groups )
        reference_cpgs = count_reference_cpgs.out.sites

        depth_dependent_prep( ready_methylkits, reference_cpgs )

        // The cross-library analyses need two libraries.
        parquets = depth_dependent_prep.out.cpgs
            .map { library, parquet, _hist -> tuple(library, parquet) }
            .toSortedList { a, b -> naturalCompare(a[0], b[0]) }
            .filter { rows -> rows.size() >= 2 }
            .map { rows -> rows.collect { row -> row[1] } }

        depth_dependent_merge( parquets, reference_cpgs )

        depth_hists = depth_dependent_prep.out.cpgs.map { _library, _parquet, hist -> hist }.collect()
        // ifEmpty turns "merge never ran" into an empty list, so the report still runs
        // with only the depth tab.
        merge_stats = depth_dependent_merge.out.stats.ifEmpty([])
        depth_dependent_read_stats( ready_bams )
        read_stats = depth_dependent_read_stats.out.stats.collect().ifEmpty([])
        depth_dependent_report( depth_hists, merge_stats, feature_files, read_stats, depth_check )

    emit:
        html = depth_dependent_report.out.html
        depth_summary = depth_dependent_report.out.depth_summary
}

// Natural order for library names, matching bin/library_order.py: numbers compare by value
// and a mass unit converts to ng, so 0.1ng < 1ng < 10ng < 100ng and 100pg < 1ng. Only a
// number with a unit can have a decimal point, so S1.10 sorts after S1.2.
def naturalTokens(name) {
    return (name =~ /(?i)(\d+(?:\.\d+)?)\s*(pg|ng|ug|µg|μg)(?![a-z])|(\d+)|\D+/).collect { g ->
        if (g[1] != null) {
            def scale = [pg: 0.001, ng: 1, ug: 1000, 'µg': 1000, 'μg': 1000][g[2].toLowerCase()]
            return [1, '', (g[1] as BigDecimal) * scale]
        }
        if (g[3] != null) {
            return [1, '', g[3] as BigDecimal]
        }
        return [0, g[0].toLowerCase(), 0]
    }
}

def naturalCompare(a, b) {
    def ta = naturalTokens(a as String)
    def tb = naturalTokens(b as String)
    def first = (0..<Math.min(ta.size(), tb.size())).findResult { i ->
        def c = (ta[i][0] <=> tb[i][0]) ?: (ta[i][1] <=> tb[i][1]) ?: (ta[i][2] <=> tb[i][2])
        return c ? c : null
    }
    return first ?: (ta.size() <=> tb.size())
}

// Log each library's shallowest checked contig and return true when at least minLibraries pass.
def depthDependentReady(checkTsv, minLibraries) {
    def lines = checkTsv.readLines()
    def header = lines[0].split('\t') as List
    def rows = lines.drop(1).collect { line -> [header, line.split('\t') as List].transpose().collectEntries() }
    def passing = rows.count { row -> row.passes == 'yes' }
    def ready = passing >= minLibraries
    def threshold = rows ? rows[0].min_contig_depth : '?'
    def table = rows.collect { row ->
        String.format('  %-60s shallowest %-8s %7.1fx  %s', row.library, row.shallowest_contig,
            row.shallowest_mean_depth as double, row.passes == 'yes' ? 'passes' : 'too shallow')
    }
    log.info "Depth check: a library passes when all ${rows ? rows[0].n_checked_contigs : 0} " +
        "checked contigs have mean depth > ${threshold}x; the analyses need ${minLibraries} passing:\n" +
        table.join('\n')
    log.info ready
        ? "Depth-dependent: ${passing} of ${rows.size()} libraries pass; running the depth-dependent analyses on all ${rows.size()}."
        : "Depth-dependent: ${passing} of ${rows.size()} libraries pass; skipping the depth-dependent analyses."
    return ready
}
