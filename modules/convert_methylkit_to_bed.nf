process convert_methylkit_to_bed {
    label 'process_single'
    tag "${library}"
    conda "conda-forge::pigz=2.8 conda-forge::gawk=5.3.1 conda-forge::sed=4.9"

    input:
        tuple val(library), path(cytosine_report)
        val(genome_fa)
        val(genome_fai)

    output:
        tuple val(library), path('*.bed'), emit: methylkit_bed

    script:
    // Using cytosine report as a substitute for merging 3 methylkit files. 
    // Reports C's with no coverage, so we will filter those out
    """
    # Input: chromosome, position (1-based), strand, meth_count, unmeth_count, context, trinucleotide
    # Output: chr, start (0-based), end (1-based), context, methylation_freq, strand
    awk '
        BEGIN { OFS = "\t" }
        !( \$4==0 && \$5==0 ){
            chromosome = \$1
            position = \$2
            strand = \$3
            meth_count = \$4
            unmeth_count = \$5
            context = \$6
            
            # Calculate methylation frequency
            total_count = meth_count + unmeth_count
            meth_freq = (total_count > 0) ? (meth_count / total_count * 100) : 0
            
            # Convert to BED format (0-based start, 1-based end)
            start = position - 1
            end = position
            
            print chromosome, start, end, context, meth_freq, strand
        }
    ' ${cytosine_report} > "${library}.bed"
    """
}
