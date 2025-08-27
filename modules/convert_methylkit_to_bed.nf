process convert_methylkit_to_bed {
    label 'process_single'
    tag "${library}"
    conda "conda-forge::pigz=2.8 conda-forge::gawk=5.3.1 conda-forge::sed=4.9"

    input:
        tuple val(library), path(methylkit_CHH_gz), path(methylkit_CHG_gz), path(methylkit_CpG_gz)
        val(genome_fa)
        val(genome_fai)

    output:
        tuple val(library), path('*.bed'), emit: methylkit_bed

    script:
    """
    methylkit_basename=\$(basename "${methylkit_CpG_gz}" _CpG.methylKit.gz)

    # create sed scripts to interconvert chr names from genome index with line numbers for efficient sorting
    awk '{print "s/\\t"\$1"\\t/\\t"NR-1"\\t/"}' "${genome_fai}" > convert_chr_to_num.sed
    awk '{print "s/^"NR-1"\\t/"\$1"\\t/"}' "${genome_fai}" > revert_num_to_chr.sed
    # merge initial methylkit files into a single BED format
    # Convert methylKit to BED format (0-based coordinates)
    # methylKit format: chr.base\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT
    # Output format: chr\tstart\tend\tname\tmethylation\tstrand
    sort -m -k2,2n -k3,3n \
            <(pigz -d -c "${methylkit_CHH_gz}" | sed 's/\$/\\tCHH/' | sed -f convert_chr_to_num.sed) \
            <(pigz -d -c "${methylkit_CHG_gz}" | sed 's/\$/\\tCHG/' | sed -f convert_chr_to_num.sed) \
            <(pigz -d -c "${methylkit_CpG_gz}" | sed 's/\$/\\tCpG/' | sed -f convert_chr_to_num.sed) \
    | awk -v OFS="\\t" 'NR>3{print \$2,\$3-1,\$3,\$8,\$6,(\$4 == "F") ? "+" : "-"}' \
    | sed -f revert_num_to_chr.sed \
    > "\${methylkit_basename}.bed"
    """
}
