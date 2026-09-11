// One organism's GC bias curve, streamed out of a composite alignment so no subset BAM is written.
//
// picard_reference MUST be this group's subset reference, never the composite: CollectGcBiasMetrics
// takes its window denominator from the reference it is handed.
process gc_bias_by_contig_group {
    tag { "${library}:${group}" }
    label 'medium_cpu'
    conda "bioconda::picard=3.3.0 bioconda::samtools=1.22 conda-forge::gawk=5.3.1 conda-forge::sed=4.9"
    // pattern keeps the multiqc/ copy out of the published dir, where it would collide with the
    // whole-reference ${library}.gc_metrics from gc_bias.
    publishDir "${params.outputDir}/stats/gc_bias", pattern: "*.gc_metrics"

    input:
        tuple val(library), path(bam), path(bai), val(group), val(picard_reference)

    output:
        tuple val(library), val(group), path("${library}.${group}.gc_metrics"), emit: for_agg
        tuple val(library), val(group), path("multiqc/${library}.gc_metrics"), emit: for_multiqc
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions
        tuple val("${task.process}"), val('picard'), eval('picard CollectGcBiasMetrics --version 2>&1 | cut -f 2 -d ":"'), topic: versions

    script:
    def prefix = "${library}.${group}"
    // Keep the picard flags below in step with modules/gc_bias.nf.
    """
    # The contigs to keep and the regions to slice are just the subset reference's .fai, reshaped,
    # so they cannot disagree with the FASTA picard normalizes against.
    cut -f 1 ${picard_reference}.fai > group_contigs.txt
    awk -v OFS='\\t' '{ print \$1, 0, \$2 }' ${picard_reference}.fai > group.bed

    # picard requires the BAM header's contigs to match the reference dictionary's in name and
    # order, so comparing names catches a stale list here with a clear message rather than a picard
    # stack trace. Lengths need no check: the subset FASTA is cut from the composite.
    samtools view -H ${bam} \\
      | awk -v contigs=group_contigs.txt '
          BEGIN { while ((getline line < contigs) > 0) { keep[line] = 1 } }
          /^@SQ/ {
              name = ""
              for (i = 1; i <= NF; i++) { if (\$i ~ /^SN:/) { name = substr(\$i, 4) } }
              if (name in keep) { print name }
          }
        ' FS='\\t' > bam_group_contigs.txt
    if ! cmp -s bam_group_contigs.txt group_contigs.txt; then
        echo "ERROR: ${library} group '${group}': BAM header contigs do not match ${picard_reference}.fai" >&2
        echo "  '<' = subset reference, '>' = BAM header" >&2
        diff group_contigs.txt bam_group_contigs.txt >&2 || true
        exit 1
    fi

    # -M -L rather than region arguments: a draft assembly member can have >7000 contigs and would
    # blow the argv limit. -L keeps the whole composite header, hence the @SQ filter.
    samtools view -h -M -L group.bed ${bam} \\
      | awk -v contigs=group_contigs.txt '
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

    # The group rides in ACCUMULATION_LEVEL, which GcBias.parse in ngs-aggregate_results already
    # reads. Safe unquoted: readContigGroups constrains contig_group to [A-Za-z0-9_.-]+.
    sed -i 's/^All Reads\\t/${group}\\t/' ${prefix}.gc_metrics

    # MultiQC's picard/gcbias module uses the filename as the sample name, so this copy keeps the
    # report's labels identical to historic runs while the curve's content changes.
    mkdir -p multiqc
    cp ${prefix}.gc_metrics multiqc/${library}.gc_metrics
    """
}
