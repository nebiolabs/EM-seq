#!/usr/bin/env nextflow
 
params.genome = '/mnt/galaxy/data/genome/grch38_core+bs_controls/sam_indexes/grch38_core+bs_controls/grch38_core+bs_controls.fa'

//CPG Islands from  UCSC table browser 
//  	Database: hg38    Primary Table: cpgIslandExt    Row Count: 31,144   Data last updated: 2018-08-10
params.ucsc_cpg_islands_gtf = 'grch38_cpgIsland_ext.gtf.gz'
params.refseq_gtf = 'GRCh38_latest_genomic.gff.gz'

params.high_quality_meth_bed = '../200ng_4cycles_LB_1.downsampled.md.bam_CpG.txt.bgz'

params.bam_files_glob = '../*.downsampled.md.{bam,bam.bai}'

params.tmp_dir = '/state/partition1/sge_tmp/'
params.output_dir = 'output'
params.ncbi_assembly_report = 'GCF_000001405.39_GRCh38.p13_assembly_report.txt'
params.dfam_out_file = 'grch38_dfam405_repeat_mask.fa.out'

Channel.value(file(params.high_quality_meth_bed)).set{ hq_meth_bed }
Channel.fromFilePairs(params.bam_files_glob, checkIfExists: true).into{ bams_for_epd; bams_for_cpgs; bams_for_refseq; bams_for_dfam }
Channel.value(file(params.ucsc_cpg_islands_gtf)).set { ucsc_cpg_islands_gtf }
Channel.value(file(params.refseq_gtf)).set { refseq_gtf }
Channel.value(file(params.ncbi_assembly_report)).set { ncbi_assembly_report }
Channel.value(file(params.dfam_out_file)).set { dfam_out }
Channel.from(['promoter', 'transcriptional_cis_regulatory_region',
              'enhancer', 'mobile_genetic_element', 'primary_transcript', 
              'lnc_RNA', 'exon', 'mRNAexon1', 'mRNAexon' ]).into { refseq_feature_types; refseq_feature_types_for_gtf }


process clean_epd_gtf {
   conda "curl ucsc-bedtogenepred ucsc-genepredtogtf"

   output:
       file('grch38_promoters.gtf') into epd_promoters_gtf

   shell:
   '''
       curl -fsSL "ftp://ccg.epfl.ch/epdnew/H_sapiens/006/Hs_EPDnew_006_hg38.bed" \
       | tr ' ' '\t' \
       | bedToGenePred /dev/stdin /dev/stdout \
       | genePredToGtf file /dev/stdin /dev/stdout > grch38_promoters.gtf
   '''	
}

process epd_methylation {
    conda "bedtools=2.29.2 htslib=1.9"
    publishDir "$params.output_dir", mode: 'copy'

    input:
        file gtf from epd_promoters_gtf
        file bed from hq_meth_bed

    output:
        file 'epd_promoter_methylation.tsv' into epd_promoter_meth

    shell:
    '''
    bedtools intersect -nonamecheck \
    -wa -wb -loj \
    -a !{gtf} -b <(bgzip -d < !{bed} ) \
    | awk -v FS='\\t' -v OFS='\\t' '$14>0 {print $10,$11,$12,$1":"$4-1"-"$5,($15*1.0)/$14 }' \
    | bedtools groupby -g 4 -o mean -c 5 \
    > epd_promoter_methylation.tsv 
    '''
}

process epd_promoter_counts{
    conda "subread=2.0.0"
    cpus 16

    input:
        file gtf from epd_promoters_gtf
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

    input:
        file ucsc_cpg_gtf from ucsc_cpg_islands_gtf
    output:
        file('grch38_cpg_islands.uniqname.gtf') into cpg_islands_gtf 

    shell:
    '''
    zcat !{ucsc_cpg_gtf} \
    |  awk -v FS='\t' -v OFS='\t' '{print $1,$1":"$4"-"$5,$3,$4,$5,$6,$7,$8,$9}' \
    > grch38_cpg_islands.uniqname.gtf
    '''
}

process cpg_island_methylation {
    conda "bedtools=2.29.2 htslib=1.9"
    publishDir "$params.output_dir", mode: 'copy'

    input:
        file gtf from cpg_islands_gtf
        file bed from hq_meth_bed

    output:
        file 'cpg_island_methylation.tsv' into cpg_island_meth

    shell:
    '''
    bedtools intersect -nonamecheck \
    -wa -wb -loj \
    -a !{gtf} -b <(bgzip -d < !{bed} ) \
    | awk -v FS='\\t' -v OFS='\\t' '$14>0 {print $10,$11,$12,$1":"$4-1"-"$5,($15*1.0)/$14 }' \
    | bedtools groupby -g 4 -o mean -c 5 \
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


process refseq_feature_gtfs {
    tag {feature}
    conda "subread=2.0.0 bedtools=2.29.2"

    input:
        file(gtf) from refseq_gtf
        file(assembly_report) from ncbi_assembly_report
        val feature from refseq_feature_types_for_gtf

    output:
        tuple feature, file('*.gtf') into feature_gtf_for_meth
        tuple feature, file('*_flat.saf') into feature_saf_for_counts

    shell:
    '''
    # uses awk to create a hash lookup from the first file (NCBI assembly report) 
    # translating chr name in the second file
    awk -v OFS='\\t' -v FS='\\t' 'NR==FNR {dict[$1]=$2; next} {$1=dict[$1]; print}' \
    <(grep -v '^#' !{assembly_report} | cut -f 7,10 | tr -d '\\r')  \
    <(zcat !{gtf} | grep -v '^#') \
    | grep "GeneID:" \
    | grep -P -v "_alt\\t" \
    | grep -P -v "^na\\t" \
    | sed -r 's/;Dbxref(=[^;]*)GeneID:([^,;]+)([;,])/;gene_id=\\2;Dbxref\\1GeneID:\\2\\3/' \
    | awk  -v OFS='\\t' -v FS='\\t' \
        '($3=="exon") && (index($9,"gbkey=mRNA") > 0) && (index($9,"-1;Parent") > 0) \
           { print($1,$2,"mRNAexon1",$4,$5,$6,$7,$8,$9); next }
         ($3=="exon") && (index($9,"gbkey=mRNA") > 0)  \
           { print($1,$2,"mRNAexon",$4,$5,$6,$7,$8,$9); next }
         { print }  
        ' \
    > name_converted.gff

    #exons overlap, we want only the longest to avoid 0 cov exons from featureCounts
    flattenGTF -a name_converted.gff -o flat_name_converted.saf -t !{feature}

    #need to switch to bed for intersection later
    tail -n +2 flat_name_converted.saf \
      | awk -v OFS='\\t' -v FS='\\t' '{print $2,$3-1,$4,$1,"-",$5}' \
      | bedtools sort -faidx !{params.genome}.fai -i /dev/stdin > !{feature}_flat.bed
    
    #filters by feature type
    awk -v type=!{feature} -v OFS='\\t' -v FS='\\t' '($3==type) { print}' name_converted.gff \
    > !{feature}.gtf

    #only include those entries that intersect with the desired feature type, back to SAF format
    echo "GeneID\tChr\tStart\tEnd\tStrand" > !{feature}_flat.saf
    bedtools intersect -a !{feature}_flat.bed -b !{feature}.gtf  -u \
      | awk -v OFS='\\t' -v FS='\\t' '{print $4,$1,$2+1,$3,$6}' >> !{feature}_flat.saf
    '''
}

process refseq_feature_methylation {
    tag {feature}
    conda "bedtools=2.29.2 htslib=1.9"
    publishDir "$params.output_dir", mode: 'copy'

    input:
        file bed from hq_meth_bed
        tuple feature, file(feature_gtf) from feature_gtf_for_meth 

    output:
        file '*_methylation.tsv' into feature_methylation
        
    shell:
    '''
    bedtools intersect -nonamecheck \
    -wa -wb -loj \
    -a !{feature_gtf} -b <(bgzip -d < !{bed} ) \
    | awk -v FS='\\t' -v OFS='\\t' '$14>0 {print $10,$11,$12,$1":"$4-1"-"$5,($15*1.0)/$14 }' \
    | bedtools groupby -g 4 -o mean -c 5 \
    > !{feature}_methylation.tsv 
    '''
}


feature_saf_for_counts
    .combine(bams_for_refseq.map{ [it[1][0],it[1][1]] } ) //combination of every saf with every bam/bai pair)
    .groupTuple(by: [0,1])
    .set{feature_bams_for_refseq}

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
       > grch38_dfam405_repeat_mask.gtf
    '''
}

process dfam_feature_methylation {
    conda "bedtools=2.29.2 htslib=1.9"
    publishDir "$params.output_dir", mode: 'copy'

    input:
        file bed from hq_meth_bed
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

process combine_methylation {
    publishDir "$params.output_dir", mode: 'copy'

    input:
        file(dfam_meth) from dfam_methylation
        file(feature_meth) from feature_methylation.collect()
        file(cpg_meth) from cpg_island_meth
        file(epd_meth) from epd_promoter_meth
    output:
        file('combined_methylation.tsv') into combined_meth

    shell:
    '''
    echo 'File	Locus	Frac Methylated' > combined_methylation.tsv
    #adds a column (tab separated) containing the name of the file being processed (repeated on each line)
    for f in !{dfam_meth} !{feature_meth} !{cpg_meth} !{epd_meth} ; do
        filebase=$(basename "${f}" _methylation.tsv)
        lines=$(wc -l <(grep -ve '^\\s*$' -e '^#' "$f") | cut -f 1 -d ' ')
        paste <( yes ${filebase} | head -n $lines ) <(grep -ve '^\\s*$' -e '^#' "$f") >> combined_methylation.tsv
    done    
    '''  
}
process combine_counts {
    publishDir "$params.output_dir", mode: 'copy'

    input:
        file(dfam_counts) from dfam_feature_counts
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
    for f in !{dfam_counts} !{feature_counts} !{cpg_counts} !{epd_counts}; do
        filebase=$(basename "${f}" _counts.tsv)
        lines=$(wc -l <(grep -ve '^\\s*$' -e '^#' "$f") | cut -f 1 -d ' ')
        paste <( yes ${filebase} | head -n $lines ) <(grep -ve '^\\s*$' -e '^#' "$f") | tail -n +2 >> combined_feature_counts.tsv
    done    
    '''  
}
