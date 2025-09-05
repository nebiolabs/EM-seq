process group_bed_intersections {
    label 'process_high_memory'
    tag "${library}"
    conda "conda-forge::gawk=5.3.1"
    publishDir "${params.outputDir}/methylKit_intersections", mode: 'symlink'
    input:
        tuple val(library), path(intersect_file)

    output:
        path('*_intersections.tsv'), emit: intersections
        path('*_positional_summary.tsv'), emit: positional_summary
        tuple val(library), path('*_intersections.tsv'), path('*_summary.tsv'), emit: for_agg

    script:
    """
    intersect_basename=\$(basename "${intersect_file}" _intersect.tsv)

    output_file="\${intersect_basename}_intersections.tsv"
    summary_file="\${intersect_basename}_positional_summary.tsv"

    echo -e "methylkit_file\\tchr\\tstart\\tend\\tcontext\\tmethylation\\ttarget_locus\\ttarget_name" > \${output_file}
    echo -e "methylkit_file\tchr\ttarget_length\tposition\tcontext\tmean_methylation\tn_loci\tn_measurements" > \${summary_file}

    if [ -s "${intersect_file}" ]; then
        awk -v output_file="\${output_file}" \\
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

            locus = target_chr ":" (target_start + 1) "-" target_end
            target_length = target_end - target_start
            position_in_target = meth_start - target_start + 1

            print "${library}", meth_chr, meth_start, meth_end, meth_context, meth_value, locus, target_name >> output_file

            # Collect summary data
            key = "${library}" "\\t" target_chr "\\t" target_length "\\t" position_in_target "\\t" meth_context
            sum[key] += meth_value
            count[key]++
            loci[key][locus] = 1  # Track unique loci
        }
        END {
            # Write summary data
            for (k in sum) {
                split(k, parts, "\\t")
                methylkit_file = parts[1]
                chr = parts[2]
                target_length = parts[3]
                position = parts[4]
                context = parts[5]
                mean_meth = sum[k] / count[k]
                n_measurements = count[k]
                n_loci = length(loci[k])
                print methylkit_file, chr, target_length, position, context, mean_meth, n_loci, n_measurements >> summary_file
            }
        }' "${intersect_file}"
    else
        echo "No intersections found in ${intersect_file}"
        touch \${output_file}
        touch \${summary_file}
    fi
    """
}
