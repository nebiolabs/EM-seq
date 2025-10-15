process download_epd_promoters {
    label 'low_cpu'
    conda "conda-forge::curl conda-forge::gawk bioconda::ucsc-bedtogenepred bioconda::ucsc-genepredtogtf "
    publishDir "${params.outputDir}/features", mode: 'copy'

    input:
        val(epd_url)

    output:
        path('epd_promoter.gtf'), emit: gtf
        tuple val("${task.process}"), val('bedToGenePred'), eval('bedToGenePred 2>&1 | grep "bedToGenePred" | sed \'s/bedToGenePred - //\''), topic: versions

    script:
    """
    curl -fsSL "${epd_url}" \
    | tr ' ' '\t' \
    | bedToGenePred /dev/stdin /dev/stdout \
    | genePredToGtf file /dev/stdin /dev/stdout \
    | gawk -v OFS='\\t' -v FS='\\t' '\$3=="transcript" {\$3="epd_promoter"} {print}' \
    > epd_promoter.gtf
    """
}

process download_liftover_chain {
    label 'low_cpu'
    conda "conda-forge::curl bioconda::htslib=1.22.1"
    publishDir "${params.outputDir}/features", mode: 'copy'

    input:
        val(chain_file_url)

    output:
        path('hg38ToHs1.over.chain.gz'), emit: chain_file
        tuple val("${task.process}"), val('curl'), eval('curl --version | head -n 1 | cut -d " " -f 2'), topic: versions

    script:
    """
    curl -fsSL "${chain_file_url}" | bgzip > hg38ToHs1.over.chain.gz
    """
}

process crossmap_epd_promoters {
    label 'low_cpu'
    conda "bioconda::crossmap=0.7.0"
    publishDir "${params.outputDir}/features", mode: 'copy'

    input:
        path(epd_orig_gtf)
        path(chain_file)

    output:
        path('epd_promoter_lifted.gtf'), emit: gtf
        tuple val("${task.process}"), val('liftOver'), eval('liftOver 2>&1 | grep "liftOver" | head -n 1 || echo "v1"'), topic: versions

    script:
    """
    CrossMap gff ${chain_file} ${epd_orig_gtf} epd_promoter_lifted.gtf
    """
}

process download_cpg_islands {
    label 'low_cpu'
    conda "conda-forge::curl"
    publishDir "${params.outputDir}/features", mode: 'copy'

    input:
        val(cpg_file_url)

    output:
        path('*cpgIslandExt*'), emit: ucsc_file
        tuple val("${task.process}"), val('curl'), eval('curl --version | head -n 1 | cut -d " " -f 2'), topic: versions

    script:
    """
    curl -fsSL -O "${cpg_file_url}"
    """
}

process download_refseq_gtf {
    label 'low_cpu'
    conda "conda-forge::curl"

    input:
        val(refseq_url)

    output:
        path('refseq.gff.gz'), emit: gtf
        tuple val("${task.process}"), val('curl'), eval('curl --version | head -n 1 | cut -d " " -f 2'), topic: versions

    script:
    """
    curl -fsSL "${refseq_url}" -o refseq.gff.gz
    """
}

process download_assembly_report {
    label 'low_cpu'
    conda "conda-forge::curl"

    input:
        val(report_url)

    output:
        path('assembly_report.txt'), emit: report
        tuple val("${task.process}"), val('curl'), eval('curl --version | head -n 1 | cut -d " " -f 2'), topic: versions

    script:
    """
    curl -fsSL "${report_url}" -o assembly_report.txt
    """
}

process normalize_cpg_islands {
    label 'low_cpu'
    conda "bioconda::htslib=1.22.1 bioconda::ucsc-bigbedtobed bioconda::ucsc-bedtogenepred bioconda::ucsc-genepredtogtf"
    publishDir "${params.outputDir}/features", mode: 'copy'

    input:
        path(cpg_file)
        path(assembly_report)

    output:
        path('cpg_island.gtf'), emit: gtf
        tuple val("${task.process}"), val('bigBedToBed'), eval('bigBedToBed 2>&1 | grep "bigBedToBed" | head -n 1 || echo "v4"'), topic: versions

    script:
    def is_bigbed = cpg_file.name.endsWith('.bb')
    def is_txt = cpg_file.name.endsWith('.txt.gz')
    if (is_bigbed)
        """
        # Create chromosome mapping from assembly report (RefSeq accession -> chr name)
        grep -v '^#' ${assembly_report} \
        | awk -F'\t' '\$5 != "na" {print \$5"\\t"\$1}' \
        > chr_mapping.txt

        # Convert bigBed to BED
        bigBedToBed ${cpg_file} cpg_islands_raw.bed

        # Map RefSeq accessions to chr names and convert to GTF
        awk -v FS='\t' -v OFS='\t' 'NR==FNR {map[\$1]=\$2; next}
        {
            chrom = (\$1 in map) ? map[\$1] : \$1
            start = \$2; end = \$3; name = \$4
            id = chrom":"start"-"end
            print chrom, "cpgIslands", "cpg_island", start+1, end, ".", "+", ".", "transcript_id \\""id"\\"; gene_id \\""id"\\";"
        }' chr_mapping.txt cpg_islands_raw.bed > cpg_island.gtf
        """
    else if (is_txt)
        """
        # Handle text format (original logic)
        zcat -f ${cpg_file} \
        | awk -v FS='\t' -v OFS='\t' '{print \$1,"cpgIslands","cpg_island",\$4,\$5,\$6,\$7,\$8,\$9}' \
        > cpg_island.gtf
        """
    else
        error "Unsupported file format for CpG islands: ${cpg_file.name}"
}
process normalize_refseq_features {
    label 'medium_cpu'
    tag { feature }
    conda "bioconda::subread=2.1.1 bioconda::bedtools=2.31.1 conda-forge::grep conda-forge::gawk conda-forge::gzip"
    publishDir "${params.outputDir}/features", mode: 'copy'

    input:
        tuple path(refseq_gtf), path(assembly_report), path(genome_fai), val(feature)

    output:
        tuple val(feature), path("${feature}.gtf"), emit: gtf
        tuple val("${task.process}"), val('bedtools'), eval('bedtools --version | sed \'s/^bedtools v//\''), topic: versions

    script:
    """
    # Translate chromosome names using NCBI assembly report
    gawk -v OFS='\t' -v FS='\t' 'NR==FNR {dict[\$1]=\$2; next} {\$1=dict[\$1]; print}' \
    <(grep -v '^#' ${assembly_report} | cut -f 7,10 | tr -d '\\r')  \
    <(zcat -f < ${refseq_gtf} | grep -v '^#') \
    | grep "GeneID:" \
    | grep -P -v "_alt\\t" \
    | grep -P -v "^na\\t" \
    | sed -r 's/;Dbxref(=[^;]*)GeneID:([^,;]+)([;,])/;gene_id=\\2;Dbxref\\1GeneID:\\2\\3/' \
    | gawk -v OFS='\\t' -v FS='\\t' \
        '(\$3=="exon") && (index(\$9,"gbkey=mRNA") > 0) && (index(\$9,"-1;Parent") > 0) \
           { print(\$1,\$2,"mRNAexon1",\$4,\$5,\$6,\$7,\$8,\$9); next }
         (\$3=="exon") && (index(\$9,"gbkey=mRNA") > 0)  \
           { print(\$1,\$2,"mRNAexon",\$4,\$5,\$6,\$7,\$8,\$9); next }
         { print }
        ' \
    > name_converted.gff

    # Flatten overlapping features to avoid 0 coverage from featureCounts
    flattenGTF -a name_converted.gff -o flat_name_converted.saf -t ${feature}

    # Convert to BED for intersection
    tail -n +2 flat_name_converted.saf \
      | gawk -v OFS='\\t' -v FS='\\t' '{print \$2,\$3-1,\$4,\$1,"-",\$5}' \
      | bedtools sort -faidx ${genome_fai} -i /dev/stdin > ${feature}_flat.bed

    # Filter by feature type
    gawk -v type=${feature} -v OFS='\\t' -v FS='\\t' '(\$3==type) { print }' name_converted.gff \
    > ${feature}.gtf

    # Only include entries that intersect with the desired feature type
    echo "GeneID\tChr\tStart\tEnd\tStrand" > ${feature}_flat.saf
    bedtools intersect -a ${feature}_flat.bed -b ${feature}.gtf -u \
      | gawk -v OFS='\\t' -v FS='\\t' '{print \$4,\$1,\$2+1,\$3,\$6}' >> ${feature}_flat.saf
    """
}

process download_dfam_annotations {
    label 'low_cpu'
    conda "conda-forge::curl conda-forge::grep conda-forge::gawk bioconda::ucsc-bedtogenepred bioconda::ucsc-genepredtogtf"
    publishDir "${params.outputDir}/features", mode: 'copy'

    input:
        val(dfam_url)

    output:
        path('dfam.gtf'), emit: gtf
        tuple val("${task.process}"), val('genePredToGtf'), eval('genePredToGtf 2>&1 | grep "genePredToGtf" | sed \'s/genePredToGtf - //\''), topic: versions

    script:
    """
    curl -fsSL "${dfam_url}" \
    | grep -v '^#' \
    | gawk -v OFS='\\t' -v FS='\\t' '\$17==1 {
        score = int(\$16 * 20)
        if (score < 0) score = 0
        if (score > 1000) score = 1000
        start = \$10 < \$11 ? \$10 - 1 : \$11 - 1
        end = \$10 < \$11 ? \$11 : \$10
        print \$1, start, end, \$3, score, \$9
    }' \
    | bedToGenePred /dev/stdin /dev/stdout \
    | genePredToGtf file /dev/stdin /dev/stdout \
    | gawk -v OFS='\\t' -v FS='\\t' '{print \$1,"dfam_repeats", \$3, \$4, \$5, \$6, \$7, \$8, \$9}' \
    | sed -r 's/\\ttranscript\\t/\\tdfam_repeat\\t/' \
    > dfam.gtf
    """
}
