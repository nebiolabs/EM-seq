// Picard normalizes coverage against the window-GC distribution of the reference it is handed
//
// Assignment is by exact contig name, read from the built contig_groups.tsv

def hasContigGroups(reference_list) {
    return reference_list?.get('gc_groups_dir') as boolean
}

// One map per group, contigs in the TSV's order, which is composite .fai order.
def readContigGroups(genome, reference_list) {
    def dir = reference_list.gc_groups_dir
    def tsv = file("${dir}/contig_groups.tsv")
    if (!tsv.exists()) {
        error "Genome '${genome}' sets gc_groups_dir='${dir}' but ${tsv} is missing."
    }

    // [:] is a LinkedHashMap, so group order follows first appearance in the TSV.
    def grouped = [:]
    tsv.readLines().drop(1).findAll { it?.trim() }.each { line ->
        // -1 keeps trailing empty fields: an ungrouped contig's row ends in an empty
        // contig_group, and the default split would drop it and make the row look truncated.
        def fields = line.split('\t', -1)
        def group = fields[3]?.trim()
        if (!group) {
            return  // ungrouped contigs get no curve of their own
        }
        grouped.computeIfAbsent(group) { [contigs: [], bp: 0L] }
        grouped[group].contigs << fields[0]
        grouped[group].bp += fields[1] as long
    }

    def groups = grouped.collect { name, data ->
        [
            name       : name,
            contigs    : data.contigs,
            bp         : data.bp,
            fasta      : file("${dir}/${name}.fa"),
            fai        : file("${dir}/${name}.fa.fai"),
            dict       : file("${dir}/${name}.dict"),
            bed        : file("${dir}/${name}.bed"),
            contig_list: file("${dir}/${name}.contigs.txt"),
        ]
    }

    validateContigGroupArtifacts(genome, dir, groups)
    return groups
}

// Checked up front because a missing subset reference otherwise surfaces as an opaque picard
// failure mid-run.
def validateContigGroupArtifacts(genome, dir, groups) {
    def missing = groups.collectMany { group ->
        ['fasta', 'fai', 'dict', 'bed', 'contig_list']
            .findAll { !group[it].exists() }
            .collect { group[it].toString() }
    }
    if (missing) {
        error "Genome '${genome}' is missing pre-built per-organism GC bias references in " +
              "${dir}:\n  ${missing.join('\n  ')}"
    }

    def group_summary = groups.collect { "${it.name} (${it.contigs.size()} contigs, ${it.bp} bp)" }
    log.info "Per-organism GC bias enabled for '${genome}': ${group_summary.join(', ')}"
}
