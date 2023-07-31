#!/usr/bin/env nextflow
 
params.genome = '/mnt/galaxy/data/genome/grcm39+meth_controls/bwameth_index/grcm39+meth_controls/grcm39+meth_controls.fa'

//CPG Islands from  UCSC table browser 
//  	Database: mm39    Primary Table: cpgIslandExt    Row Count: 15,9665   Data downloaded: 2023-01-08
params.ucsc_cpg_islands_gtf = '/mnt/home/langhorst/nebnext_projects/em-seq/mouse/grcm39_cpg_islands.gtf.gz'

params.refseq_gff_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/annotation_releases/109/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz'
// methylkit file in percents format:
    //chrBase chr     base    strand  coverage        freqC   freqT	
    //CM000994.3.3050095      CM000994.3      3050095 F       8         0.00  100.00
params.high_quality_mk_file = ''

params.bam_files_glob = '*.md.{bam,bam.bai}'

params.tmp_dir = '/state/partition1/sge_tmp/'
params.output_dir = 'output'
params.ncbi_assembly_report_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/annotation_releases/109/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_assembly_report.txt'
params.epd_promoter_bed_url = 'https://epd.expasy.org/ftp/epdnew/M_musculus/003/Mm_EPDnew_003_mm10.bed'
params.mm10_mm39_chain_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz'
//params.dfam_out_file = 'grch38_dfam405_repeat_mask.fa.out'

Channel.fromPath(params.high_quality_meth_bed, checkIfExists: true).first().set { hq_methylkit }
Channel.fromFilePairs(params.bam_files_glob, checkIfExists: true).into{ bams_for_epd; bams_for_cpgs; bams_for_refseq; bams_for_dfam }
Channel.fromPath(params.ucsc_cpg_islands_gtf, checkIfExists: true).first().set { ucsc_cpg_islands_gtf }
Channel.value(params.ncbi_assembly_report_url).set { ncbi_assembly_report_url }
//Channel.value(file(params.dfam_out_file)).set { dfam_out }
Channel.value(params.refseq_gff_url).set { refseq_gff_url }
Channel.value(params.epd_promoter_bed_url).set { epd_promoter_bed_url }
Channel.value(params.mm10_mm39_chain_url).set { mm10_mm39_chain_url }
  
// 'mobile_genetic_element' not present in mouse col 3
Channel.from(['transcriptional_cis_regulatory_region', 'enhancer','promoter',
              'primary_transcript', 'snoRNA','snRNA', 'tRNA',
              'lnc_RNA', 'exon', 'mRNAexon1', 'mRNAexon' ]).into { refseq_feature_types; refseq_feature_types_for_gff }

process fetch_chain_file {
  conda "curl"

  input: 
    val url from mm10_mm39_chain_url
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

process methylkit_to_bed {
  conda "gawk gzip" 
  
  input: 
    val mk_file from hq_methylkit
  output:
    file('*.bed.gz') into hq_methylkit_bed
  
  //chrBase chr     base    strand  coverage        freqC   freqT	
  //CM000994.3.3050095      CM000994.3      3050095 F       8         0.00  100.00
  //format for next step: CM000994.3      3050094    3050095 8   0.00    100.00
  shell:
  '''
    zcat -f !{mk_file} | awk -v FS='\\t' -v OFS='\\t' '{ print $2,$3-1,$3,$5, $6, $7 }'| gzip > hq_methylkit.bed.gz
  '''

}

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
       > epd_promoters_mm10.gtf && \
       CrossMap.py gff !{chain_file} epd_promoters_mm10.gtf epd_promoters_mm39.gtf
       
     awk -v OFS='\\t' -v FS='\\t' 'NR==FNR {dict[$1]=$2; next} {$1=dict[$1]; print}' \
       <(grep -v '^#' !{assembly_report} | awk -v OFS='\\t' -v FS='\\t' '{print $10,$5}' | tr -d '\\r')  \
       <(zcat -f epd_promoters_mm39.gtf | grep -v '^#') > epd_promoters.gtf
   '''	
}

process epd_methylation {
    conda "bedtools=2.29.2 htslib=1.9"
    publishDir "$params.output_dir", mode: 'copy'

    input:
        file(gtf) from epd_promoters_gtf
        file(bed) from hq_methylkit_bed

    output:
        file 'epd_promoter_methylation.tsv' into epd_promoter_meth

    shell:
    // coverage > 0
    // %C / 100 (methylation proportion), average over sites within the feature
    '''
    bedtools intersect -nonamecheck \
      -wa -wb -loj \
      -a !{gtf} -b <(bgzip -d < !{bed} ) \
      | awk -v FS='\\t' -v OFS='\\t' '$13>0 {print $10,$11,$12,$1":"$4-1"-"$5,($14/100)}' \
      | sort -k4,4 | bedtools groupby -g 4 -o mean -c 5 \
      > epd_promoter_methylation.tsv 
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
    featureCounts --primary --ignoreDup -Q 10 -M -f -o -O --fraction -p -P -B -C \
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
      <(grep -v '^#' !{assembly_report} | awk -v OFS='\\t' -v FS='\\t' '{print $10,$5}' | tr -d '\\r')  \
      <(zcat -f !{ucsc_cpg_gtf} | grep -v '^#') \
    |  awk -v FS='\t' -v OFS='\t' '{print $1,$1":"$4"-"$5,$3,$4,$5,$6,$7,$8,$9}' \
    > cpg_islands.uniqname.gtf
    '''
}

process cpg_island_methylation {
    conda "bedtools=2.29.2 htslib=1.9"
    publishDir "$params.output_dir", mode: 'copy'

    input:
        file gtf from cpg_islands_gtf
        file bed from hq_methylkit_bed

    output:
        file 'cpg_island_methylation.tsv' into cpg_island_meth

    shell:
    '''
    bedtools intersect -nonamecheck \
    -wa -wb -loj \
    -a !{gtf} -b <(bgzip -d < !{bed} ) \
    | awk -v FS='\\t' -v OFS='\\t' '$13>0 {print $10,$11,$12,$1":"$4-1"-"$5,($14/100) }' \
    | sort -k4,4 | bedtools groupby -g 4 -o mean -c 5 \
    > cpg_island_methylation.tsv 
    '''
}

process cpg_island_counts{
    conda "subread=2.0.0"
    cpus 16

    input:
        file gtf from cpg_islands_gtf
        path('*') from bams_for_cpgs.map{ [it[1][0],it[1][1]] }.flatten().toList()

    output:
        file 'cpg_island_counts.tsv' into cpg_island_counts

    shell:
    '''
    featureCounts --primary --ignoreDup -Q 10 -M -f -o -O --fraction -p -P -B -C \
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
     <(grep -v '^#' !{assembly_report} | awk -v OFS='\\t' -v FS='\\t' '{print $7,$5}' | tr -d '\\r')  \
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


process refseq_feature_methylation {
    tag {feature}
    conda "bedtools=2.29.2 htslib=1.9"
    publishDir "$params.output_dir", mode: 'copy'

    input:
        file bed from hq_methylkit_bed
        tuple feature, file(feature_gff) from feature_gff_for_meth 

    output:
        file '*_methylation.tsv' into feature_methylation
        
    shell:
    // gff 9 columns
    // e.g. CM000994.3      cmsearch        exon    3172239 3172348 .       +       .       ID=exon-XR_004936710.1-1;Parent=rna-XR_004936710.1;gene_id=115487594;Dbxref=GeneID:115487594,RFAM:RF00026,Genbank:XR_004936710.1,MGI:MGI:5455983;gbkey=ncRNA;gene=Gm26206;inference=COORDINATES: profile:INFERNAL:1.1.1;product=U6 spliceosomal RNA;transcript_id=XR_004936710.1
    // hq_methylkit: chr, base-1, base, cov, %C, %T
    // chr, start, end, chr:base-base, methylation proportion -> groupby column 4
    '''
    bedtools intersect -nonamecheck \
    -wa -wb -loj \
    -a !{feature_gff} -b <(bgzip -d < !{bed} ) \
    | awk -v FS='\\t' -v OFS='\\t' '$13>0 {print $10,$11,$12,$1":"$4-1"-"$5,($14/100) }' \
    | sort -k4,4 | bedtools groupby -g 4 -o mean -c 5 \
    > !{feature}_methylation.tsv 
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
    featureCounts --primary --ignoreDup -Q 10 -M -f -O --fraction -p -P -B -C \
    -a !{feature_saf} -F SAF\
    -t !{feature} \
    -g 'ID' \
    --tmpDir !{params.tmp_dir} \
    -T !{task.cpus} \
    -o !{feature}_counts.tsv *.bam 
    '''
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
        featureCounts --primary --ignoreDup -Q 10 -M -f -o -O --fraction -p -P -B -C \
        -a !{gtf} \
        -t transcript \
        -g 'transcript_id' \
        --tmpDir !{params.tmp_dir} \
        -T !{task.cpus} \
        -o dfam_counts.tsv *.bam 
    '''
}
*/

process combine_methylation {
    publishDir "$params.output_dir", mode: 'copy'

    input:
        //file(dfam_meth) from dfam_methylation
        file(feature_meth) from feature_methylation.collect()
        file(cpg_meth) from cpg_island_meth
        file(epd_meth) from epd_promoter_meth
    output:
        file('combined_methylation.tsv') into combined_meth

    shell:
    '''
    echo 'File	Locus	Frac Methylated' > combined_methylation.tsv
    #adds a column (tab separated) containing the name of the file being processed (repeated on each line)
    for f in !{feature_meth} !{cpg_meth} !{epd_meth} ; do
        filebase=$(basename "${f}" _methylation.tsv)
        lines=$(wc -l <(grep -ve '^\\s*$' -e '^#' "$f") | cut -f 1 -d ' ')
        paste <( yes ${filebase} | head -n $lines ) <(grep -ve '^\\s*$' -e '^#' "$f") >> combined_methylation.tsv
    done    
    '''  
}

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

