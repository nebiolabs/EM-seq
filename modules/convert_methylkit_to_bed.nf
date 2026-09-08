process convert_methylkit_to_bed {
    label 'cpus_8'
    tag "${library}"
    conda "conda-forge::parallel=20250622 bioconda::mawk=1.3.4"

    input:
        tuple val(library), path(cytosine_report)
        val(genome_fa)
        val(genome_fai)

    output:
        tuple val(library), path("${library}.bed"), emit: methylkit_bed

    script:
    // Using cytosine report as a substitute for merging 3 methylkit files.
    // Reports C's with no coverage, so we will filter those out.
    //
    // Rows are independent, so `parallel --pipepart` seeks into byte ranges of the
    // file with no copy and runs mawk over each concurrently; `--keep-order` keeps
    // the output in the sorted order `bedtools intersect -sorted` needs downstream.
    """
    # A script file, not an inline awk program: `parallel` re-quotes the command it
    # hands each worker, which mangles braces and \$N as its own syntax.
    cat <<'AWK' > convert_methylkit.awk
    BEGIN { OFS = "\t" }
    # Input: chromosome, position (1-based), strand, meth_count, unmeth_count, context, trinucleotide
    !( \$4==0 && \$5==0 ){
        total_count = \$4 + \$5
        meth_freq = (total_count > 0) ? (\$4 / total_count * 100) : 0
        # Output: chr, start (0-based), end (1-based), context, methylation_freq, strand
        print \$1, \$2 - 1, \$2, \$6, meth_freq, \$3
    }
    AWK

    parallel --pipepart --keep-order -j ${task.cpus} -a ${cytosine_report} \\
        mawk -f convert_methylkit.awk > "${library}.bed"
    """
}
