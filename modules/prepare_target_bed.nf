/*
 * Prepare Target Bed
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
    label 'process_single'
    tag "${target_bed.baseName}"
    conda "bioconda::bedtools=2.31.1"

    input:
        path(target_bed)
        val(slop_len)
        val(genome_fa)
        val(genome_fai)

    output:
        path('*_slop_sorted.bed')

    script:
    """
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
    bedtools slop -g "${genome_fai}" -b ${slop_len} > ${target_bed.baseName}_slop_sorted.bed

    if [ ! -s ${target_bed.baseName}_slop_sorted.bed ]; then
        echo "Error: No valid regions found after processing ${target_bed}" >&2
        exit 1
    fi

    """
}
