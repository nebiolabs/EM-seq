process group_bed_intersections {
    label 'process_high_memory'
    tag "${library}"
    conda "conda-forge::gawk=5.3.1"
    publishDir "${params.outputDir}/methylKit_intersections", mode: 'symlink'
    input:
        tuple val(library), path(intersect_file)

    output:
        path('*_intersections.tsv'), emit: results
        path('*_summary.tsv'), emit: summary
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
