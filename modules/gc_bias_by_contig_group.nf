// One organism's GC bias curve from a composite alignment, computed without ever writing the
// subset BAM to disk.
//
// CollectGcBiasMetrics is a single-pass reader and needs no index, so samtools can stream the
// group's reads into it on stdin -- the same idiom gc_bias itself used before per-organism curves
// existed. Materializing the subset BAM instead costs roughly one extra copy of the aligned BAM
// per library (the groups partition the reference), plus its index and an indexing pass, for a
// file nothing downstream reads.
//
// picard_reference MUST be this group's subset reference, never the composite: CollectGcBiasMetrics
// takes its window denominator from the whole reference it is handed, so the composite's other
// organisms would drag this curve toward their own GC.
process gc_bias_by_contig_group {
    tag { "${library}:${group}" }
    label 'medium_cpu'
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.22 conda-forge::gawk=5.3.1 conda-forge::sed=4.9"
    // pattern keeps the multiqc/ copy below out of the published directory, where it would
    // collide with the whole-reference ${library}.gc_metrics from gc_bias.
    publishDir "${params.outputDir}/stats/gc_bias", pattern: "*.gc_metrics"

    input:
        tuple val(library), path(bam), path(bai), val(group), path(group_bed), path(group_contigs), val(picard_reference)

    output:
        tuple val(library), val(group), path("${library}.${group}.gc_metrics"), emit: for_agg
        tuple val(library), val(group), path("multiqc/${library}.gc_metrics"), emit: for_multiqc
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions
        tuple val("${task.process}"), val('picard'), eval('picard CollectGcBiasMetrics --version 2>&1 | cut -f 2 -d ":"'), topic: versions

    script:
    def prefix = "${library}.${group}"
    // The picard flags here must stay in step with modules/gc_bias.nf, which runs the same tool for
    // the whole-reference curve.
    """
    # picard compares the BAM header's sequence dictionary against the reference's and demands the
    # same names, lengths AND order. Checking that up front turns a bad contig list into a clear
    # message instead of a stack trace inside CollectGcBiasMetrics; counting @SQ lines alone would
    # accept a list that merely disagrees on order, which picard then rejects.
    samtools view -H ${bam} \\
      | awk -v contigs=${group_contigs} '
          BEGIN { while ((getline line < contigs) > 0) { keep[line] = 1 } }
          /^@SQ/ {
              name = ""
              for (i = 1; i <= NF; i++) { if (\$i ~ /^SN:/) { name = substr(\$i, 4) } }
              if (name in keep) { print name }
          }
        ' FS='\\t' > bam_group_contigs.txt
    cut -f 1 ${picard_reference}.fai > reference_contigs.txt
    if ! cmp -s bam_group_contigs.txt reference_contigs.txt; then
        echo "ERROR: ${library} group '${group}': BAM header contigs do not match ${picard_reference}.fai" >&2
        echo "  '<' = subset reference, '>' = BAM header" >&2
        diff reference_contigs.txt bam_group_contigs.txt >&2 || true
        exit 1
    fi

    # -M -L rather than naming regions as arguments: a draft assembly member can have >7000 contigs
    # and would blow the argv limit. -L keeps the full composite header, so the @SQ filter is
    # required rather than cosmetic -- picard rejects a header listing every composite contig
    # against a subset reference.
    samtools view -h -M -L ${group_bed} ${bam} \\
      | awk -v contigs=${group_contigs} '
          BEGIN { while ((getline line < contigs) > 0) { keep[line] = 1 } }
          /^@SQ/ {
              name = ""
              for (i = 1; i <= NF; i++) { if (\$i ~ /^SN:/) { name = substr(\$i, 4) } }
              if (name in keep) { print }
              next
          }
          { print }
        ' FS='\\t' OFS='\\t' \\
      | picard -Xmx${task.memory.toGiga()}g CollectGcBiasMetrics \\
          --IS_BISULFITE_SEQUENCED true --VALIDATION_STRINGENCY SILENT \\
          -I /dev/stdin -O ${prefix}.gc_metrics -S ${prefix}.gc_summary_metrics \\
          --CHART /dev/null -R ${picard_reference}

    # The group rides in ACCUMULATION_LEVEL because GcBias.parse in ngs-aggregate_results already
    # reads that column and maps picard's own labels to 'all'. Anchoring on 'All Reads' leaves the
    # ACCUMULATION_LEVEL header line alone. The group needs no escaping: readContigGroups in
    # lib/contig_groups.nf constrains contig_group to [A-Za-z0-9_.-]+, so it cannot reach here
    # containing sed metacharacters.
    sed -i 's/^All Reads\\t/${group}\\t/' ${prefix}.gc_metrics

    # MultiQC's picard/gcbias module is configured with use_filename_as_sample_name, so the file
    # name becomes the series label. Copying to the name the whole-reference curve has always used
    # keeps the report's sample labels identical to historic runs -- only the curve's content
    # changes -- so the two are directly comparable.
    mkdir -p multiqc
    cp ${prefix}.gc_metrics multiqc/${library}.gc_metrics
    """
}
