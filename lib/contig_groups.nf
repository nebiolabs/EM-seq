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

    // Columns are read positionally, so a file whose columns have been reordered must fail loudly
    // here rather than silently grouping by the wrong field.
    def lines = tsv.readLines()
    def required = ['contig_name', 'sequence_length', 'frac_gc', 'contig_group']
    def header = lines ? lines[0].split('\t', -1)*.trim() : []
    if (header.size() < required.size() || header[0..<required.size()] != required) {
        error "${tsv}: first four columns must be ${required.join(', ')}; got ${header}"
    }

    // A generated gc_groups_dir is well-formed by construction, but a custom genome's TSV is
    // written by hand against a reference this pipeline never built, so the rules below are
    // enforced here rather than assumed.
    // These are the contig types the aggregation database accepts; anything else is rejected only
    // at load time, after the run has already finished.
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
            bed        : file("${dir}/${name}.bed"),
            contig_list: file("${dir}/${name}.contigs.txt"),
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

// A gc_groups_dir built from this same FASTA already agrees with it, but a TSV written by hand, or
// carried over from a similar reference, has no such guarantee: picard compares each subset BAM's
// header against the subset reference's dictionary, so a stale or mistyped contig otherwise
// surfaces as a dictionary error once per library per group, after alignment has already run.
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
        ['fasta', 'fai', 'dict', 'bed', 'contig_list']
            .findAll { !group[it].exists() }
            .collect { group[it].toString() }
    }
    if (missing) {
        error "Genome '${genome}' is missing pre-built per-organism GC bias references in " +
              "${dir}:\n  ${missing.join('\n  ')}"
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

// MultiQC plots one GC bias curve per library, and its picard/gcbias module collapses a multi-block
// file to a single series, so it has to be handed exactly one curve. The whole-reference curve is
// the wrong one to show: it counts the spike-in controls, and pUC19 (51% GC), lambda (50%) and
// especially Xp12 (68%) pull the visible curve away from the autosome-only curve EM-seq reported
// before per-organism curves existed -- a ~24% shift in picard's GC_DROPOUT at a 1% spike-in.
//
// Defaults to the group with the most picard windows rather than a fixed name like
// 'human_autosome', because five of the configured genomes are not human; the largest group is the
// host organism for every '+meth_controls' composite. --multiqc_gc_group names one explicitly.
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
