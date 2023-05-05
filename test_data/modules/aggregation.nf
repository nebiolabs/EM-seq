path_to_ngs_agg = "/mnt/bioinfo/prg/ngs-aggregate_results_dev/current/"
//"/Users/aerijman/Documents/new_ngs/newer/ngs-aggregate_results"
library = params.library

process aggregate_emseq {
    cpus 1
    conda "samtools=1.9"
    publishDir "${library}/ngs-agg"

    input:
         tuple val(library), path(mbias), path(bam), path(bai), val(barcodes)

    output:
        path('ngs_agg.*')

    // barcode should be split by "-" as bc1-bc2 
    shell:
    '''
    genome_name=$(echo !{params.genome} | awk -F"/" '{print $NF}' | sed 's/.fa|.fasta//')
    # ADD THIS bc=!{barcodes} AND THEN bc1= split $bc on "_" and take $1

    bc1=$(echo !{barcodes}- | cut -f 1 -d "-")
    bc2=$(echo !{barcodes}- | cut -f 2 -d "-")
    if [[ "${bc2}" == "" ]]; then
        bc2_arg=""
    else
        bc2_arg="--barcode2 ${bc2}"
    fi
    # unzip *fastqc.zip
    export RBENV_VERSION=$(cat !{path_to_ngs_agg}/.ruby-version)
    DATABASE_ENV=production !{path_to_ngs_agg}/bin/bundle exec !{path_to_ngs_agg}/aggregate_results.rb \
    --bam !{bam} \
    --bai !{bai} \
    --name "!{params.library}" \
    --barcode1 ${bc1} ${bc2_arg} \
    --lane !{params.lane} \
    --contact_email !{params.email} \
    --project "!{params.project}" \
    --sample "!{params.sample}" \
    --genome "${genome_name}" \
    --combined_mbias_records !{mbias} 2> ngs_agg.err 1> ngs_agg.out 
    '''
    //--workflow "Automated EM-seq!{dest_modifier}"
}
    



//tuple val(email), 
//val(library), 
//val(project), 
//val(sample), 
//val(barcode), 
//val(lane), 
//val(genome_name), 
//val(genome_fasta), 
//path(bam), 
//path(bai), 
///* path(fastqc), 
//path(flagstats), 
//path(idxstats), 
//path(noncontrol_gc), */
//path(nonconverted), 
////path(noncontrol_stats), 
//path(mbias), 
///* path(tasmanian), 
//path(picard_metrics), 
//path(metadata_fastq) // from aggregates_for_aggregation
//*/


//bc1=$(echo !{barcode}- | cut -f 1 -d "-")
//bc2=$(echo !{barcode}- | cut -f 2 -d "-")
//if [[ "${bc2}" == "" ]]; then
//bc2_arg=""
//else
//bc2_arg="--barcode2 ${bc2}"
//fi
//unzip *fastqc.zip
//export RBENV_VERSION=$(cat !{path_to_ngs_agg}/.ruby-version)
//DATABASE_ENV=production !{path_to_ngs_agg}/bin/bundle exec !{path_to_ngs_agg}/aggregate_results.rb \
//--bam !{bam} \
//--bai *.md.bai \
//--name "!{library}" \
//--barcode1 ${bc1} ${bc2_arg} \
//--lane !{lane} \
//--contact_email !{email} \
//--project "!{project}" \
//--sample "!{sample}" \
//--genome "!{genome_name}" \
////            --gc *gc_metrics \
////            --idx_stats *idxstat \
////            --flagstat *flagstat \
//--nonconverted_read_counts *nonconverted-counts.tsv \
//--combined_mbias_records !{mbias}/*mbias.tsv \
////            --fastqc *_fastqc/fastqc_data.txt \
////            --insert *insertsize_metrics \
////            --tasmanian *.csv \
//--aln *alignment_summary_metrics.txt \
////            --metadatafq_file !{metadata_fastq} \
//--workflow "Automated EM-seq!{dest_modifier}"
