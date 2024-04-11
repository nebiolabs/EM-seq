path_to_ngs_agg = "/mnt/bioinfo/prg/ngs-aggregate_results/current/"
//"/Users/aerijman/Documents/new_ngs/newer/ngs-aggregate_results"

process aggregate_emseq {
    cpus 1
    conda "samtools=1.9"
    publishDir "${params.flowcell}/${library}/ngs-agg"

    input:
         tuple val(email), val(library), val(barcodes), path(metadata_fastq), path(nonconverted_counts_tsv), path(fastp), 
               path(bam), path(bai), 
               path(gc_metrics),
               path(idxstat),
               path(flagstat),
               path(fastqc_zip),
               path(insertsize_metrics),
               path(tasmanian),
               path(mbias),
               path(alignment_summary_metrics_txt)            

    output:
        path('ngs_agg.*')

    // barcode should be split by "-" as bc1-bc2 
    shell:
    '''

    # Workaround that has to be deleted later: make workflow different for each library.

    genome_name=$(echo !{params.genome} | awk -F"/" '{print $NF}' | sed 's/.fa|.fasta//')

    # bc = barcode1 + barcode2 if exists.
    bc=$(echo !{barcodes} | tr -d "][" | awk -F"-" '{bc2=""; if (length($2)==length($1)) {bc2="--barcode2 "$2}; print $1" "bc2;}')

    unzip *fastqc.zip

    cat !{nonconverted_counts_tsv} | awk -v l=!{library} '{print l"\t"$0}' > !{library}.nonconverted_counts.for_agg.tsv
    
    export RBENV_VERSION=$(cat !{path_to_ngs_agg}/.ruby-version)
    DATABASE_ENV=production !{path_to_ngs_agg}/bin/bundle exec !{path_to_ngs_agg}/aggregate_results.rb \
    --bam !{bam} \
    --bai !{bai} \
    --name !{library} \
    --barcode1 ${bc} \
    --lane !{params.lane} \
    --contact_email !{params.email} \
    --project !{params.project} \
    --sample !{params.sample} \
    --genome !{params.genome} \
    --gc !{gc_metrics} \
    --idx_stats !{idxstat} \
    --flagstat !{flagstat} \
    --nonconverted_read_counts !{library}.nonconverted_counts.for_agg.tsv \
    --combined_mbias_records !{mbias} \
    --fastqc *_fastqc/fastqc_data.txt \
    --insert !{insertsize_metrics} \
    --tasmanian !{tasmanian} \
    --aln !{alignment_summary_metrics_txt} \
    --metadatafq_file !{metadata_fastq} \
    --fastp !{fastp} \
    --workflow !{params.workflow}_!{library} 2> ngs_agg.err 1> ngs_agg.out
    '''
    // add number of chimeras!
}
