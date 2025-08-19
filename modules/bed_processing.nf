/*
 * BED Processing Module
 *
 * This module standardizes BED files to BED6 format (chr, start, end, name, score, strand)
 * to simplify downstream processing and avoid dynamic column detection issues.
 *
 * Input BED files are validated and standardized as follows:
 * - Minimum 3 columns required (chr, start, end)
 * - Column 4 (name): uses existing value or generates chr:start-end
 * - Column 5 (score): uses existing value or defaults to "0"
 * - Column 6 (strand): uses existing value or defaults to "."
 * - Validates coordinates and strand values
 */

process prepare_target_bed {
    label 'low_cpu'
    tag "${target_bed.baseName}"
    conda "bioconda::bedtools=2.31.1"

    input:
        path(target_bed)
        tuple path(genome_fa), path(genome_fai)

    output:
        path('*_slop_sorted.bed'), emit: prepared_bed

    script:
    """
    # Apply slop to target BED file (50bp on each side)
    slop_len=50
    target_basename=\$(basename "${target_bed}" .bed)

    echo "Standardizing and applying \${slop_len} bp slop to \${target_basename}..."

    # First, standardize to BED6 format (chr, start, end, name, score, strand)
    awk 'BEGIN { OFS="\\t" }
    !/^#/ {
        # Validate minimum required columns (chr, start, end)
        if (NF < 3) {
            print "Error: BED file must have at least 3 columns (chr, start, end)" > "/dev/stderr"
            exit 1
        }

        # Validate coordinates
        if (\$2 !~ /^[0-9]+\$/ || \$3 !~ /^[0-9]+\$/) {
            print "Error: Invalid coordinates in line: " \$0 > "/dev/stderr"
            exit 1
        }

        if (\$2 >= \$3) {
            print "Error: Start coordinate must be less than end coordinate in line: " \$0 > "/dev/stderr"
            exit 1
        }

        # Build BED6 format
        chr = \$1
        start = \$2
        end = \$3
        name = (NF >= 4 && \$4 != "") ? \$4 : chr ":" start "-" end
        score = (NF >= 5 && \$5 != "") ? \$5 : "0"
        strand = (NF >= 6 && \$6 != "") ? \$6 : "."

        # Validate strand
        if (strand != "+" && strand != "-" && strand != ".") {
            print "Error: Invalid strand value: " strand > "/dev/stderr"
            exit 1
        }

        print chr, start, end, name, score, strand
    }' "${target_bed}" | \\
    bedtools sort -g "${genome_fai}" -i /dev/stdin | \\
    bedtools slop -g "${genome_fai}" -b \${slop_len} > \${target_basename}_slop_sorted.bed

    # Verify output is not empty
    if [ ! -s \${target_basename}_slop_sorted.bed ]; then
        echo "Error: No valid regions found after processing ${target_bed}" >&2
        exit 1
    fi

    echo "Successfully processed \$(wc -l < \${target_basename}_slop_sorted.bed) regions"
    """
}

process intersect_bed_with_methylkit {
    label 'low_cpu'
    tag "${library}"
    conda "bioconda::bedtools=2.31.1"

    input:
        tuple val(library), path(methylkit_bed)
        path(target_bed_prepared)
        tuple path(genome_fa), path(genome_fai)

    output:
        tuple val(library), path('*_intersect.tsv'), emit: intersections

    script:
    """
    methylkit_basename=\$(basename "${methylkit_bed}" .bed)
    target_basename=\$(basename "${target_bed_prepared}" _slop_sorted.bed)

    # Uses bedtools intersect with -wa -wb to keep both annotations
    # nonamecheck since control contigs don't conform
    bedtools intersect -nonamecheck \\
        -a "${target_bed_prepared}" \\
        -b "${methylkit_bed}" \\
        -wa -wb -sorted -g "${genome_fai}" > \${methylkit_basename}_vs_\${target_basename}_intersect.tsv
    """
}

process group_bed_intersections {
    label 'low_cpu'
    tag "${library}"
    conda "conda-forge::gawk=5.3.1"
    publishDir "${params.outputDir}/methylKit_intersections", mode: 'symlink'
    input:
        tuple val(library), path(intersect_file)

    output:
        path('*_intersections.tsv'), emit: intersection_results
        path('*_summary.tsv'), emit: intersection_summary
        tuple val(library), path('*_intersections.tsv'), path('*_summary.tsv'), emit: for_agg

    script:
    """
    # Extract base names from the intersect file
    intersect_basename=\$(basename "${intersect_file}" _intersect.tsv)

    output_file="\${intersect_basename}_intersections.tsv"
    summary_file="\${intersect_basename}_summary.tsv"

    # Add headers
    echo -e "methylkit_file\\tchr\\tstart\\tend\\tcontext\\tmethylation\\ttarget_locus\\ttarget_name" > \${output_file}
    echo -e "methylkit_file\\ttarget_length\\tposition\\tcontext\\tmean_methylation\\tn_loci\\tn_measurements" > \${summary_file}

    # Extract methylkit basename from the intersect filename
    methylkit_basename=\$(echo "\${intersect_basename}" | sed 's/_vs_.*//')

    # Process intersection results with awk
    # Target BED is standardized above to 6 columns: chr, start, end, name, score, strand
    # MethylKit BED has 5 columns: chr, start, end, context, methylation_value
    if [ -s "${intersect_file}" ]; then
        awk -v methylkit_basename="\${methylkit_basename}" \\
            -v output_file="\${output_file}" \\
            -v summary_file="\${summary_file}" '
        BEGIN {
            OFS="\\t"
        }
        !/^#/ {
            # Target BED columns (always 6 columns)
            target_chr = \$1
            target_start = \$2
            target_end = \$3
            target_name = \$4
            target_score = \$5
            target_strand = \$6

            # MethylKit BED columns (starting at column 7)
            meth_chr = \$7
            meth_start = \$8
            meth_end = \$9
            meth_context = \$10
            meth_value = \$11

            # Create locus string (convert 0-based BED to 1-based for display)
            locus = target_chr ":" (target_start + 1) "-" target_end

            # Calculate target length and position within target
            target_length = target_end - target_start
            # Position relative to start of target (1-based)
            position_in_target = meth_start - target_start + 1

            # Output detailed results
            print methylkit_basename, meth_chr, meth_start, meth_end, meth_context, meth_value, locus, target_name >> output_file

            # Collect summary data
            key = methylkit_basename "\\t" target_length "\\t" position_in_target "\\t" meth_context
            sum[key] += meth_value
            count[key]++
            loci[key][locus] = 1  # Track unique loci
        }
        END {
            # Write summary data
            for (k in sum) {
                split(k, parts, "\\t")
                methylkit_file = parts[1]
                target_length = parts[2]
                position = parts[3]
                context = parts[4]
                mean_meth = sum[k] / count[k]
                n_measurements = count[k]
                n_loci = length(loci[k])
                print methylkit_file, target_length, position, context, mean_meth, n_loci, n_measurements >> summary_file
            }
        }' "${intersect_file}"
    else
        echo "No intersections found in ${intersect_file}"
        # Create empty files with headers only
        touch \${output_file}
        touch \${summary_file}
    fi
    """
}

process concatenate_intersections {
    label 'low_cpu'
    publishDir "${params.outputDir}/methylKit_intersections", mode: 'copy'

    input:
        path(intersection_files)
        path(summary_files)

    output:
        path('all_intersections_combined.tsv'), emit: combined_intersections
        path('all_summaries_combined.tsv'), emit: combined_summaries

    script:
    """
    # Combine all intersection files
    echo -e "methylkit_file\\tchr\\tstart\\tend\\tcontext\\tmethylation\\ttarget_locus\\ttarget_name" > all_intersections_combined.tsv

    # Add all intersection data (skip headers from individual files)
    for file in ${intersection_files}; do
        if [ -s "\$file" ]; then
            tail -n +2 "\$file" >> all_intersections_combined.tsv
        fi
    done

    # Combine all summary files
    echo -e "methylkit_file\\ttarget_length\\tposition\\tcontext\\tmean_methylation\\tn_loci\\tn_measurements" > all_summaries_combined.tsv

    # Add all summary data (skip headers from individual files)
    for file in ${summary_files}; do
        if [ -s "\$file" ]; then
            tail -n +2 "\$file" >> all_summaries_combined.tsv
        fi
    done

    # Report statistics
    total_intersections=\$(tail -n +2 all_intersections_combined.tsv | wc -l)
    total_summaries=\$(tail -n +2 all_summaries_combined.tsv | wc -l)

    echo "Combined \${total_intersections} intersection records from \$(echo ${intersection_files} | wc -w) files"
    echo "Combined \${total_summaries} summary records from \$(echo ${summary_files} | wc -w) files"
    """
}
