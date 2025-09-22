Channel.fromPath(params.inputdir + "*.methylKit.gz").into{methylkit_files; for_bins}
ref = params.reference_fasta
annotation = params.annotation
ref_len = params.reference_lengths

process make_bed {
    cpus 5
    conda "bioconda::deeptools"
    errorStrategy 'finish'

    input:
        file methylkit from methylkit_files

    output:
        file('*.bed') into for_bigwig

    shell:
    '''
    lib=$(basename !{methylkit} .methylKit.gz)
    zcat !{methylkit} | tail -n +2 | awk '{ if ($5 >= 8) { print } }' \
    | awk '{ printf("%s\\t%s\\t%s\\t%s\\n", $2, $3 - 1, $3, $6); }' \
    | LC_COLLATE=C sort -k1,1 -k2,2n > ${lib}.bed
    '''

}

process make_bigwig {
    cpus 1
    conda "bioconda::ucsc-bedgraphtobigwig"
    publishDir "/mnt/home/mcampbell/20200317_new_emseq_figure", mode: "copy"
    errorStrategy 'finish'

    input:
        file bed from for_bigwig

    output:
        file ('*.bw') into _bigwig

    shell:
    '''
    lib=$(basename !{bed} .bed)
    bedGraphToBigWig !{bed} !{ref_len} ${lib}.bw
    '''

}

process binned_figure {
    cpus 1
    publishDir "/mnt/home/mcampbell/20200317_new_emseq_figure", mode: "copy"
    errorStrategy 'finish'

    input:
        file methylkit from for_bins
    
    output:
        file ('*.tab') into _bins

    shell:
    '''
    lib=$(basename !{methylkit} .methylKit.gz)
    TSS_cpg_bins.py --methylkit !{methylkit} --fasta !{ref} --prefix ${lib} --annotation !{annotation}
    '''
}