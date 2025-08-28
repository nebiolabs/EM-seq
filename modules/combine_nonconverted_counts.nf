process combine_nonconverted_counts {
    label 'process_single'
    publishDir "${params.outputDir}/stats/nonconverted_counts"

    input:
        tuple val(library), path(nonconverted_counts)

    output:
        tuple val(library), path("${library}.nonconverted_counts.for_agg.tsv"), emit: for_agg

    script:
    """
    # Combine nonconverted counts for each chunk by library

    awk -v lib="${library}" '
        {
            chr=\$1; cnt=\$2;
            sum[chr] += cnt
        }
        END {
            for (c in sum) {
                print lib "\\t" c "\\t" sum[c]
            }
        }
    ' ${nonconverted_counts} > ${library}.nonconverted_counts.for_agg.tsv
    
    """
}
