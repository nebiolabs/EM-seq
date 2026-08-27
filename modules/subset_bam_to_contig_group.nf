// The @SQ filter matches exact contig names from the pre-built list, never a prefix (see
// lib/contig_groups.nf).
process subset_bam_to_contig_group {
    tag { "${library}:${group}" }
    label 'process_single'
    conda "bioconda::samtools=1.22"

    input:
        tuple val(library), path(bam), path(bai), val(group), path(group_bed), path(group_contigs)

    output:
        tuple val(library), val(group), path("${library}.${group}.bam"), path("${library}.${group}.bam.bai"), emit: bam
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions

    script:
    """
    # -M -L rather than regions as arguments: a draft assembly member can have >7000 contigs and
    # would blow the argv limit.
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
      | samtools view -b -o ${library}.${group}.bam -
    samtools index ${library}.${group}.bam

    # The subset header must contain exactly this group's contigs, or picard's dictionary
    # comparison against ${group}.fa fails deep inside the next process.
    expected=\$(sort -u ${group_contigs} | wc -l | tr -d ' ')
    actual=\$(samtools view -H ${library}.${group}.bam | grep -c '^@SQ' || true)
    if [ "\$expected" != "\$actual" ]; then
        echo "ERROR: ${library}.${group}.bam header has \$actual @SQ lines, expected \$expected" >&2
        exit 1
    fi
    """
}
