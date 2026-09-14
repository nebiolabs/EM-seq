// Picard normalizes coverage against the window-GC distribution of the reference it is handed, so
// each group gets its own. Contigs are assigned by exact name from the built contig_groups.tsv.

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

    // Columns are read positionally, so a file whose columns have been reordered must fail loudly
    // here rather than silently grouping by the wrong field.
    def lines = tsv.readLines()
    def required = ['contig_name', 'sequence_length', 'frac_gc', 'contig_group']
    def header = lines ? lines[0].split('\t', -1)*.trim() : []
    if (header.size() < required.size() || header[0..<required.size()] != required) {
        error "${tsv}: first four columns must be ${required.join(', ')}; got ${header}"
    }

    // Validated here rather than assumed: a custom genome's TSV is hand-written, and an illegal
    // contig_type is otherwise only caught at aggregation load time, after the run has finished.
    def legal_types = ['autosomal', 'sex_chromosome', 'organelle', 'unplaced', 'control', 'unknown']

    // [:] is a LinkedHashMap, so group order follows first appearance in the TSV.
    def grouped = [:]
    def group_of = [:]   // every contig the TSV names, grouped or not, for the .fai cross-check
    def length_of = [:]
    lines.drop(1).eachWithIndex { line, idx ->
        def lineno = idx + 2
        if (!line?.trim() || line.startsWith('#')) {
            return
        }
        // -1 keeps trailing empty fields: an ungrouped contig's row ends in an empty
        // contig_group, and the default split would drop it and make the row look truncated.
        def fields = line.split('\t', -1)
        if (fields.size() < required.size()) {
            error "${tsv} line ${lineno}: expected at least ${required.size()} tab-separated " +
                  "fields, got ${fields.size()}: ${fields}"
        }
        def contig = fields[0]?.trim()
        if (!contig) {
            error "${tsv} line ${lineno}: contig_name must not be empty"
        }
        if (group_of.containsKey(contig)) {
            error "${tsv} line ${lineno}: contig '${contig}' already appears earlier in the file; " +
                  "a contig belongs to exactly one group"
        }
        if (!(fields[1]?.trim() ==~ /\d+/)) {
            error "${tsv} line ${lineno}: sequence_length '${fields[1]}' is not a whole number"
        }
        def type = fields.size() > 4 && fields[4]?.trim() ? fields[4].trim() : 'unknown'
        if (!(type in legal_types)) {
            error "${tsv} line ${lineno}: contig_type '${type}' is not one of ${legal_types.join(', ')}"
        }
        def group = fields[3]?.trim()
        group_of[contig] = group
        length_of[contig] = fields[1].trim() as long
        if (!group) {
            return  // ungrouped contigs get no curve of their own
        }
        // The group becomes a filename and an unquoted sed replacement in
        // gc_bias_by_contig_group, so a name outside this set could inject sed or shell syntax.
        if (!(group ==~ /[A-Za-z0-9_.-]+/)) {
            error "${tsv} line ${lineno}: contig_group '${group}' must match [A-Za-z0-9_.-]+"
        }
        grouped.computeIfAbsent(group) { [contigs: [], bp: 0L, windows: 0L] }
        grouped[group].contigs << contig
        grouped[group].bp += length_of[contig]
        // Picard slides its 100 bp window per position rather than tiling, so a contig of length L
        // yields L-101 windows and one <= 101 bp yields none.
        grouped[group].windows += Math.max(0L, length_of[contig] - 101L)
    }

    if (!group_of) {
        error "${tsv} has a header but no data rows."
    }

    def groups = grouped.collect { name, data ->
        [
            name       : name,
            contigs    : data.contigs,
            bp         : data.bp,
            windows    : data.windows,
            fasta      : file("${dir}/${name}.fa"),
            fai        : file("${dir}/${name}.fa.fai"),
            dict       : file("${dir}/${name}.dict"),
        ]
    }

    // Every contig ungrouped would leave groupTuple(size: 0) waiting forever, which surfaces as
    // a run that "succeeds" having aggregated nothing.
    if (groups.isEmpty()) {
        error "Genome '${genome}' sets gc_groups_dir='${dir}' but ${tsv} assigns no contig to a group."
    }

    validateAgainstReferenceIndex(genome, tsv, reference_list, group_of, length_of)
    validateContigGroupArtifacts(genome, dir, groups)
    return groups
}

// A TSV carried over from a similar reference otherwise surfaces as a picard dictionary error once
// per library per group, after alignment has already run.
def validateAgainstReferenceIndex(genome, tsv, reference_list, group_of, length_of) {
    // genome_fai is optional in a genome entry -- main.nf only ever checks for '<genome_fa>.fai',
    // which is what samtools faidx produces -- so fall back to that rather than dying on a null.
    def fai_path = reference_list?.get('genome_fai') ?: "${reference_list.genome_fa}.fai"
    def fai = file(fai_path)
    if (!fai.exists()) {
        error "Genome '${genome}' needs a FASTA index to check ${tsv} against the reference, but " +
              "${fai} is missing. Run: samtools faidx ${reference_list.genome_fa}"
    }

    def fai_length = [:]
    fai.readLines().each { line ->
        if (!line?.trim()) {
            return
        }
        def fields = line.split('\t')
        fai_length[fields[0]] = fields[1] as long
    }

    def listed = { items -> items.take(10).join('\n  ') + (items.size() > 10 ? "\n  ..." : '') }

    def missing = group_of.keySet().findAll { !fai_length.containsKey(it) }
    if (missing) {
        error "${tsv} names ${missing.size()} contig(s) that are not in ${fai}:\n  ${listed(missing)}"
    }

    def mismatched = group_of.keySet()
        .findAll { length_of[it] != fai_length[it] }
        .collect { "${it}: TSV says ${length_of[it]}, ${fai.name} says ${fai_length[it]}" }
    if (mismatched) {
        error "${tsv} disagrees with ${fai} on sequence_length, so it was built against a " +
              "different reference:\n  ${listed(mismatched)}"
    }

    def absent = fai_length.keySet().findAll { !group_of.containsKey(it) }
    if (absent) {
        log.warn "${absent.size()} contig(s) in ${fai} are absent from ${tsv} and get no GC curve " +
                 "of their own: ${absent.take(10).join(', ')}${absent.size() > 10 ? ', ...' : ''}"
    }
}

// Checked up front because a missing subset reference otherwise surfaces as an opaque picard
// failure mid-run.
def validateContigGroupArtifacts(genome, dir, groups) {
    def missing = groups.collectMany { group ->
        ['fasta', 'fai', 'dict']
            .findAll { !group[it].exists() }
            .collect { group[it].toString() }
    }
    if (missing) {
        error "Genome '${genome}' is missing pre-built per-organism GC bias references in " +
              "${dir}:\n  ${missing.join('\n  ')}"
    }

    // gc_bias_by_contig_group works off the subset .fai; the TSV is what sizes the group for
    // --multiqc_gc_group. Disagreement means the directory was built against a different TSV.
    def divergent = groups.findAll { group ->
        def fai_contigs = group.fai.readLines().findAll { it?.trim() }.collect { it.split('\t')[0] }
        fai_contigs != group.contigs
    }.collect {
        "${it.fai.name} does not list the ${it.contigs.size()} contig(s) assigned to " +
        "'${it.name}', in that order"
    }
    if (divergent) {
        error "Genome '${genome}': per-organism GC bias references in ${dir} disagree with " +
              "contig_groups.tsv:\n  ${divergent.join('\n  ')}"
    }

    groups.findAll { it.windows == 0L }.each {
        log.warn "Group '${it.name}' yields no picard GC windows (every contig is <= 101 bp); " +
                 "its GC bias curve will be empty."
    }

    def group_summary = groups.collect {
        "${it.name} (${it.contigs.size()} contigs, ${it.bp} bp, ${it.windows} windows)"
    }
    log.info "Per-organism GC bias enabled for '${genome}': ${group_summary.join(', ')}"
}

// MultiQC has to be handed exactly one curve, and the whole-reference one is the wrong choice: the
// spike-in controls (Xp12 is 68% GC) shift picard's GC_DROPOUT ~24% at a 1% spike-in.
//
// Defaults to the group with the most windows -- the host organism in every '+meth_controls'
// composite -- rather than a fixed name, since most of the configured genomes are not human.
def selectMultiqcGroup(groups, requested) {
    if (requested) {
        def match = groups.find { it.name == requested }
        if (!match) {
            error "--multiqc_gc_group '${requested}' is not a contig group of this genome; " +
                  "available groups: ${groups*.name.join(', ')}"
        }
        log.info "MultiQC GC bias curve: group '${match.name}' (${match.windows} picard windows), " +
                 "set by --multiqc_gc_group."
        return match.name
    }

    def widest = groups.max { it.windows }
    log.info "MultiQC GC bias curve: group '${widest.name}' (${widest.windows} picard windows, the " +
             "largest). Override with --multiqc_gc_group."
    return widest.name
}
