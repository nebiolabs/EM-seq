#!/usr/bin/env nextflow

// methylkit file in percents format:
//chrBase chr     base    strand  coverage        freqC   freqT	
//CM000994.3.3050095      CM000994.3      3050095 F       8         0.00  100.00
params.mk_files = '*.methylKit.gz'
params.bam_files_glob = '*.md.{bam,bam.bai}'
params.tmp_dir = '/tmp/'
params.output_dir = 'cov_vs_meth.output'
local_ref_files_path = '/mnt/home/langhorst/nebnext_projects/em-seq/em-seq_ref_files'
params.count_dup_reads = false
if (params.count_dup_reads) {
    feature_count_dup_option = ''
} else {
    feature_count_dup_option = '--ignoreDup'
}
params.mouse = false
if (params.mouse) {
    params.genome = local_ref_files_path + '/grcm39+meth_controls.fa'
    //CPG Islands from  UCSC table browser 
    //  	Database: mm39    Primary Table: cpgIslandExt    Row Count: 15,9665   Data downloaded: 2023-01-08
    params.ucsc_cpg_islands_gtf = local_ref_files_path + '/grcm39_cpg_islands.gtf.gz'
    cpg_chr_lookup = '$10,$5' //switch from chr to genbank naming
    params.refseq_gff_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/annotation_releases/109/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz'
    params.ncbi_assembly_report_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/annotation_releases/109/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_assembly_report.txt'
    refseq_chr_lookup = '$7,$5' //switch from NC -> genbank naming
    params.epd_promoter_bed_url = 'https://epd.expasy.org/ftp/epdnew/M_musculus/003/Mm_EPDnew_003_mm10.bed'
    params.old_new_chain_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz'
}

params.human_t2t2 = false    
if (params.human_t2t2) {
    params.genome = local_ref_files_path + '/T2T_chm13v2.0+bs_controls.fa'
    //CPG Islands from  UCSC table browser 
    //Database: hub_3671779_hs1    Primary Table: hub_3671779_cpgIslandExt Data last updated: 2022-02-06
    //Item Count: 30,616
    //could not find a stable link to this annotation...  get a gtf formatted file this file interactively from 
    // https://genome.ucsc.edu/cgi-bin/hgTables?https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1849378212_nz80NzvGAMarzhkhsaEcT896jhoM&clade=mammal&org=Human&db=hs1&hgta_group=regulation&hgta_track=hub_3671779_clinVar20220313&hgta_table=0&hgta_regionType=genome

    params.ucsc_cpg_islands_gtf = local_ref_files_path + '/t2t2_ucsc_cpg_islands.gtf.gz'
    cpg_chr_lookup = '$10,$10' //just use column 10 for T2T
    params.refseq_gff_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/110/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz'
    params.ncbi_assembly_report_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/110/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt'
    refseq_chr_lookup = '$7,$10' //switch from NC -> chr1 naming
    params.epd_promoter_bed_url = 'https://epd.expasy.org/ftp/epdnew/human/006/Hs_EPDnew_006_hg38.bed'
    params.old_new_chain_url = 'https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hg38-chm13v2.over.chain.gz'
    //params.dfam_out_file = 'grch38_dfam405_repeat_mask.fa.out' //TODO: add dfam support
}

Channel.fromPath(params.mk_files).map{ it -> tuple(it.baseName, it.baseName.split(".methylKit")[0].split("_")[-1], it)}.filter{it[1] == "CpG"}.set {methylkits}
Channel.fromFilePairs(params.bam_files_glob, checkIfExists: true).into{ bams_for_epd; bams_for_cpgs; bams_for_refseq; bams_for_dfam }
Channel.fromPath(params.ucsc_cpg_islands_gtf, checkIfExists: true).first().set { ucsc_cpg_islands_gtf }
Channel.value(params.ncbi_assembly_report_url).set { ncbi_assembly_report_url }
//Channel.value(file(params.dfam_out_file)).set { dfam_out }
Channel.value(params.refseq_gff_url).set { refseq_gff_url }
Channel.value(params.epd_promoter_bed_url).set { epd_promoter_bed_url }
Channel.value(params.old_new_chain_url).set { old_new_chain_url }
  
// 'mobile_genetic_element' not present in mouse col 3
Channel.from(['transcriptional_cis_regulatory_region', 'enhancer','promoter',
              'primary_transcript', 'snoRNA','snRNA', 'tRNA',
              'lnc_RNA', 'exon', 'mRNAexon1', 'mRNAexon' ]).into { refseq_feature_types; refseq_feature_types_for_gff }

process fetch_chain_file {
  conda "curl"

  input: 
    val url from old_new_chain_url
  output: 
    file('*.chain.gz') into chain_file

  shell:
  '''
    curl -O "!{url}"
  '''
}

process fetch_refseq_assembly_report {
   conda "curl"

   input:
     val url from ncbi_assembly_report_url
   output:
     file('*.txt') into ncbi_assembly_report
   shell:
   '''
     curl -fsSl !{url} > assembly_report.txt
   ''' 
}

// process methylkit_to_bed {
//   conda "gawk gzip" 
  
//   input: 
//     val mk_file from hq_methylkit
//   output:
//     file('*.bed.gz') into hq_methylkit_bed
  
//   //chrBase chr     base    strand  coverage        freqC   freqT	
//   //CM000994.3.3050095      CM000994.3      3050095 F       8         0.00  100.00
//   //format for next step: CM000994.3      3050094    3050095 8   0.00    100.00
//   shell:
//   '''
//     zcat -f !{mk_file} | awk -v FS='\\t' -v OFS='\\t' '{ print $2,$3-1,$3,$5, $6, $7 }'| gzip > hq_methylkit.bed.gz
//   '''

// }

process clean_epd_gtf {
   conda "curl ucsc-bedtogenepred ucsc-genepredtogtf crossmap"

   input: 
       val(url) from epd_promoter_bed_url
       file(assembly_report) from ncbi_assembly_report
       file(chain_file) from chain_file

   output:
       file('epd_promoters.gtf') into epd_promoters_gtf

   shell:
   '''
     
     curl -fsSL "!{url}" \
       | tr ' ' '\t' \
       | bedToGenePred /dev/stdin /dev/stdout \
       | genePredToGtf file /dev/stdin /dev/stdout \
       > epd_promoters_oldref.gtf && \
       CrossMap gff !{chain_file} epd_promoters_oldref.gtf epd_promoters_newref.gtf
       
     awk -v OFS='\\t' -v FS='\\t' 'NR==FNR {dict[$1]=$2; next} {$1=dict[$1]; print}' \
       <(grep -v '^#' !{assembly_report} | awk -v OFS='\\t' -v FS='\\t' '{print $10,$5}' | tr -d '\\r')  \
       <(zcat -f epd_promoters_newref.gtf | grep -v '^#') > epd_promoters.gtf
   '''	
}

process epd_promoter_counts{
    conda "subread=2.0.0"
    cpus 16

    input:
        file(gtf) from epd_promoters_gtf
        path('*') from bams_for_epd.map{ [it[1][0],it[1][1]] }.flatten().toList()

    output:
        file 'epd_promoter_counts.tsv' into epd_promoter_counts

    shell:
    '''
    featureCounts --primary !{feature_count_dup_option} -Q 10 -M -f -o -O --fraction -p -P -B -C \
        -a !{gtf} \
        --tmpDir !{params.tmp_dir} \
        -T !{task.cpus} \
        -o epd_promoter_counts.tsv *.bam
    '''
}

process clean_cpg_islands_gtf {
    conda "gawk gzip"
    input:
        file ucsc_cpg_gtf from ucsc_cpg_islands_gtf
        file(assembly_report) from ncbi_assembly_report
    output:
        file('cpg_islands.uniqname.gtf') into cpg_islands_gtf 

    shell:
    '''
    awk -v OFS='\\t' -v FS='\\t' 'NR==FNR {dict[$1]=$2; next} {$1=dict[$1]; print}' \
      <(grep -v '^#' !{assembly_report} | awk -v OFS='\\t' -v FS='\\t' '{print !{cpg_chr_lookup}}' | tr -d '\\r')  \
      <(zcat -f !{ucsc_cpg_gtf} | grep -v '^#') \
    |  awk -v FS='\t' -v OFS='\t' '{print $1,$1":"$4"-"$5,$3,$4,$5,$6,$7,$8,$9}' \
    > cpg_islands.uniqname.gtf
    '''
}

process cpg_island_counts{
    conda "subread=2.0.0"
    publishDir "$params.output_dir", mode: 'copy'
    cpus 16

    input:
        file gtf from cpg_islands_gtf
        path('*') from bams_for_cpgs.map{ [it[1][0],it[1][1]] }.flatten().toList()

    output:
        file 'cpg_island_counts.tsv' into cpg_island_counts

    shell:
    '''
    featureCounts --primary !{feature_count_dup_option} -Q 10 -M -f -o -O --fraction -p -P -B -C \
        -a !{gtf} \
        --tmpDir !{params.tmp_dir} \
        -T !{task.cpus} \
        -o cpg_island_counts.tsv *.bam
    '''
}

process refseq_feature_download {
   conda "curl"

   input:
     val url from refseq_gff_url
   output:
     file('*.gff') into refseq_gff
   shell:
   '''
     curl -fsSl !{url} | zcat -f > refseq.gff
   ''' 
}

process refseq_feature_gffs {
    tag {feature}
    conda "subread=2.0.0 bedtools=2.29.2"

    input:
        file(gff) from refseq_gff
        file(assembly_report) from ncbi_assembly_report
        val feature from refseq_feature_types_for_gff

    output:
        tuple feature, file("${feature}.gff") into feature_gff_for_meth
        tuple feature, file('*_flat.saf') into feature_saf_for_counts

    shell:
    '''
    awk -v OFS='\\t' -v FS='\\t' 'NR==FNR {dict[$1]=$2; next} {$1=dict[$1]; print}' \
     <(grep -v '^#' !{assembly_report} | awk -v OFS='\\t' -v FS='\\t' '{print !{refseq_chr_lookup}}' | tr -d '\\r')  \
     <(zcat -f !{gff} | grep -v '^#') \
    | grep "GeneID:" \
    | grep -P -v "_alt\\t" \
    | grep -P -v "^na\\t" \
    | sed -r 's/;Dbxref(=[^;]*)GeneID:([^,;]+)([;,])/;gene_id=\\2;Dbxref\\1GeneID:\\2\\3/' \
    | awk  -v OFS='\\t' -v FS='\\t' \
       '($3=="exon") && (index($9,"gbkey=mRNA") > 0) && (index($9,"-1;Parent") > 0) \
          { print($1,$2,"mRNAexon1",$4,$5,$6,$7,$8,$9); next }
        ($3=="exon") && (index($9,"gbkey=mRNA") > 0) \
          { print($1,$2,"mRNAexon",$4,$5,$6,$7,$8,$9); next }
        { print }
       ' > name_converted.gff

    # exons overlap, we want only the longest to avoid 0 cov exons from featureCounts
    flattenGTF -a name_converted.gff -o flat_name_converted.saf -t !{feature}
 
    # need to switch to bed for intersection later
    tail -n +2 flat_name_converted.saf \
     | awk -v OFS='\\t' -v FS='\\t' '{print $2,$3-1,$4,$1,"-",$5}' \
     | bedtools sort -faidx !{params.genome}.fai -i /dev/stdin > !{feature}_flat.bed

    # filters by feature type
    awk -v type=!{feature} -v OFS='\\t' -v FS='\\t' '($3==type) { print}' name_converted.gff > !{feature}.gff

    # only include those entries that intersect with the desired feature type, back to SAF format
    echo "GeneID\tChr\tStart\tEnd\tStrand" > !{feature}_flat.saf
    bedtools intersect -a !{feature}_flat.bed -b !{feature}.gff  -u \
    | awk -v OFS='\\t' -v FS='\\t' '{print $4,$1,$2+1,$3,$6}' >> !{feature}_flat.saf
    '''
}

// We now look at methylation for all methylkit files not just the high quality one
// If you want to use just the "high quality" one for coverage vs methylation analysis you can
feature_gff_for_meth
    .concat(cpg_islands_gtf.map{ it -> tuple( "cpg_islands", it)})
    .concat(epd_promoters_gtf.map { it -> tuple( "epd_promoters", it)})
    .combine(methylkits)
    .set{feature_methylation_ch}

process feature_methylation {
    tag {feature}
    conda "bedtools=2.29.2 htslib=1.9"
    publishDir "$params.output_dir/features/$sample_id/$context", mode: 'copy'

    input:
        tuple val(feature), val(feature_gff), val(sample_id), val(context), path(methylkit) from feature_methylation_ch

    output:
        tuple val(sample_id), val(context), path('*_methylation.tsv') into feature_methylation_out
        
    shell:
    // gff 9 columns:  CM000994.3      cmsearch        exon    3172239 3172348 .       +       .       ID=exon-X...
    // input:
    // chrBase chr     base    strand  coverage        freqC   freqT
    //CM000994.3.3050095      CM000994.3      3050095 F       4         0.00  100.00
    // "bed" (with additional columns): chr, base0, base1, methylation fraction, coverage
    // coverage >= 5
    // chr, start, end, chr:base0-base1, methylation proportion -> groupby column 4, mean of column 5
    '''
    bedtools intersect -nonamecheck \
    -wa -wb -loj \
    -a !{feature_gff} -b <(zcat !{methylkit} | awk -v FS='\\t' -v OFS='\\t' 'NR>1 {print $2, $3-1, $3, $6/100, $5}') \
    | awk -v FS='\\t' -v OFS='\\t' '$14>=5 {print $10,$11,$12,$1":"$4-1"-"$5, $13}' \
    | sort -k4,4 | bedtools groupby -g 4 -o mean -c 5 \
    > !{feature}_methylation.tsv
    '''
}

process combine_feature_methylation {
    publishDir "$params.output_dir/features/$sample_id/$context", mode: 'copy'

    input:
        tuple val(sample_id), val(context), path(methylation_files) from feature_methylation_out.groupTuple(by: [0,1])

    output:
        path('*combined_methylation.tsv') into combined_methylation_out

    shell:
    // makes one combined methylation file for all feature types per library/context
    '''
    echo 'Feature	Locus	Meth' > !{sample_id}_!{context}_combined_methylation.tsv
    #adds a column (tab separated) containing the name of the file being processed (repeated on each line)
    for f in !{methylation_files} ; do
        filebase=$(basename "${f}" _methylation.tsv)
        lines=$(wc -l <(grep -ve '^\\s*$' -e '^#' "$f") | cut -f 1 -d ' ')
        paste <( yes ${filebase} | head -n $lines ) <(grep -ve '^\\s*$' -e '^#' "$f") >> !{sample_id}_combined_methylation.tsv
    done    
    '''  
}

feature_saf_for_counts
    .combine(bams_for_refseq.map{ [it[1][0],it[1][1]] } ) //combination of every saf with every bam/bai pair)
    .groupTuple(by: [0,1]).set{feature_bams_for_refseq}

process refseq_feature_counts {

    conda "subread=2.0.0"
    publishDir "$params.output_dir", mode: 'copy'
    cpus 16

    input:
        tuple (feature, path(feature_saf), path('*'), path('*') ) from feature_bams_for_refseq

    output:
        file '*_counts.tsv' into feature_counts

    shell:
    '''
    featureCounts --primary !{feature_count_dup_option} -Q 10 -M -f -O --fraction -p -P -B -C \
    -a !{feature_saf} -F SAF\
    -t !{feature} \
    -g 'ID' \
    --tmpDir !{params.tmp_dir} \
    -T !{task.cpus} \
    -o !{feature}_counts.tsv *.bam 
    '''
}

process feature_violin {
    cpus 1
    memory '32 GB'
    conda "python=3.11 scipy pandas matplotlib seaborn"
    publishDir "$params.output_dir/features/violin", mode: 'copy'

    input:
    path(methylation) from combined_methylation_out.collect()

    output:
    tuple path("*.png"), path("*.svg"), path("*.tsv")
    
    script:
    """
    #!/usr/bin/env python

    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import glob
    import numpy as np

    data_frames = []
    for file_path in glob.glob("*combined_methylation.tsv"):
        df = pd.read_csv(
        file_path,
        delimiter = "\t",
        header = 0,
        names = ['Feature','Locus','Meth']
        )
        df['Name'] = os.path.basename(file_path).split(".methylKit")[0].replace('_CpG','')
        data_frames.append(df)

    big_df = pd.concat(data_frames).sort_values(by=['Feature','Name'])

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16,6))

    plt.subplot(1, 1, 1)
    v = sns.violinplot(
        y=big_df["Meth"],
        x=big_df["Feature"],
        hue=big_df['Name'],
        palette = dict(zip(big_df['Name'].unique(), sns.color_palette("OrRd", n_colors=len(big_df['Name'].unique())))),
        orient="v",
        bw='scott',
        cut=0,
        legend=False,
        linewidth=0.25,
        scale='width'
    )
    
    ax.set_ylabel("Mean methylation level", fontsize=14)
    ax.tick_params(axis='y', which='major', labelsize=12)
    ax.tick_params(axis='x',labelrotation=45, labelsize=10)
    ax.set_xlabel('')
    plt.tight_layout()
    plt.savefig(
    "feature_methylation_violinplot.png",
    dpi=300
    )
    plt.savefig(
    "feature_methylation_violinplot.svg",
    format="svg"
    )

    means = big_df.groupby(['Feature','Name'],sort=False)['Meth'].mean()
    counts = big_df.groupby(['Feature','Name'],sort=False)['Meth'].count()
    means.to_csv('combined_methylation_feature_means.tsv', sep = '\t')
    counts.to_csv('combined_methylation_feature_counts.tsv', sep = '\t')
    """
}

/*
process dfam_out_to_gtf {
    conda "ucsc-bedtogenepred ucsc-genepredtogtf"
    input: 
       file rm_out from dfam_out
    output:
       file '*.gtf' into (dfam_gtf_for_meth, dfam_gtf_for_counts)

    shell:
    '''
       awk 'OFS="\t" {print($5,$6-1,$7,$11,$1,".")}' !{rm_out} \
       | tail -n +4 \
       | bedToGenePred /dev/stdin /dev/stdout \
       | genePredToGtf file /dev/stdin /dev/stdout \
       | awk -v FS='\t' -v OFS='\t' '{print $1,$1":"$4"-"$5,$3,$4,$5,$6,$7,$8,$9}' \
       > dfam_repeats.gtf
    '''
}

process dfam_feature_methylation {
    conda "bedtools=2.29.2 htslib=1.9"
    publishDir "$params.output_dir", mode: 'copy'

    input:
        file bed from hq_methylkit_bed
        file(gtf) from dfam_gtf_for_meth

    output:
        file '*_methylation.tsv' into dfam_methylation
        
    shell:
    '''
    bedtools intersect -nonamecheck \
    -wa -wb -loj \
    -a !{gtf} -b <(bgzip -d < !{bed} ) \
    | awk -v FS='\\t' -v OFS='\\t' '$14>0 {print $10,$11,$12,$1":"$4-1"-"$5,($15*1.0)/$14 }' \
    | bedtools groupby -g 4 -o mean -c 5 \
    > dfam_methylation.tsv 
    '''
}

process dfam_feature_counts {

    conda "subread=2.0.0"
    publishDir "$params.output_dir", mode: 'copy'
    cpus 16

    input:
        file(gtf) from dfam_gtf_for_counts
        path('*') from bams_for_dfam.map{ [it[1][0],it[1][1]] }.flatten().toList()

    output:
        file '*_counts.tsv' into dfam_feature_counts

    shell:
    '''
        featureCounts --primary !{feature_count_dup_option} -Q 10 -M -f -o -O --fraction -p -P -B -C \
        -a !{gtf} \
        -t transcript \
        -g 'transcript_id' \
        --tmpDir !{params.tmp_dir} \
        -T !{task.cpus} \
        -o dfam_counts.tsv *.bam 
    '''
}
*/

process combine_counts {
    publishDir "$params.output_dir", mode: 'copy'
    //sigpipe errors are expeceted with tail and head commands
    errorStrategy { task.exitStatus=141 ? 'ignore' : 'terminate' }

    input:
        //file(dfam_counts) from dfam_feature_counts
        file(feature_counts) from feature_counts.collect()
        file(cpg_counts) from cpg_island_counts
        file(epd_counts) from epd_promoter_counts
    output:
        file('combined_feature_counts.tsv') into combined_counts

    shell:
    '''
    set +o pipefail
    
    #constructs the header
    echo -n 'File\t' > combined_feature_counts.tsv
    grep -hve '^\\s*$' -e '^#' !{cpg_counts} | head -n 1 >> combined_feature_counts.tsv

    #adds a column (tab separated) containing the name of the file being processed (repeated on each line)
    for f in !{feature_counts} !{cpg_counts} !{epd_counts}; do
        filebase=$(basename "${f}" _counts.tsv)
        lines=$(wc -l <(grep -ve '^\\s*$' -e '^#' "$f") | cut -f 1 -d ' ')
        paste <( yes ${filebase} | head -n $lines ) <(grep -ve '^\\s*$' -e '^#' "$f") | tail -n +2 >> combined_feature_counts.tsv
    done    
    '''  
}

