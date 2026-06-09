nextflow.enable.dsl=2

params.mk_files = '*.methylKit.gz'
params.bam_files_glob = '*.md.{bam,bam.bai}'
params.tmp_dir =  '/tmp/'
params.output_dir = 'feature_cov.output'
params.count_dup_reads = false
params.mouse = false
params.human_t2t2 = false
params.context = 'CpG'  // methylation context(s) to analyse; comma-separated for multiple: 'CpG,CHG,CHH'

// Path to locally cached reference files used by the --mouse and --human_t2t2 shortcuts.
// Override with --local_ref_files_path for a different location.
params.local_ref_files_path = "${HOME}/nebnext_projects/em-seq/em-seq_ref_files"

// Genome-specific params declared once with null defaults.
// Set these on the CLI for a custom assembly; --mouse and --human_t2t2 pre-populate them.
// cpg_chr_lookup / refseq_chr_lookup are awk column specs (src_col,dest_col) translating
// assembly-report accession IDs to the chromosome naming used in BAMs and GTF/SAF files.
params.genome                   = null
params.ucsc_cpg_islands_gtf     = null
params.cpg_chr_lookup           = null  // e.g. '$10,$5' (mouse) or '$10,$10' (T2T)
params.refseq_gff_url           = null
params.ncbi_assembly_report_url = null
params.refseq_chr_lookup        = null  // e.g. '$7,$5' (mouse) or '$7,$10' (T2T)
params.epd_promoter_bed_url     = null
params.old_new_chain_url        = null

feature_count_dup_option = params.count_dup_reads ? '' : '--ignoreDup'

// Preset maps for known assemblies. Using local maps rather than re-assigning params
// avoids Nextflow lint warnings about params being defined multiple times.
def mouse_preset = [
    genome:                   "${params.local_ref_files_path}/grcm39+meth_controls.fa",
    ucsc_cpg_islands_gtf:     "${params.local_ref_files_path}/grcm39_cpg_islands.gtf.gz",
    cpg_chr_lookup:           '$10,$5',  // assembly report col10 (chrN) -> col5 (GenBank accession)
    refseq_gff_url:           'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/annotation_releases/109/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz',
    ncbi_assembly_report_url: 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/annotation_releases/109/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_assembly_report.txt',
    refseq_chr_lookup:        '$7,$5',   // NC_ accession -> GenBank chr name
    epd_promoter_bed_url:     'https://epd.expasy.org/ftp/epdnew/M_musculus/003/Mm_EPDnew_003_mm10.bed',
    old_new_chain_url:        'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz',
]
def t2t_preset = [
    genome:                   "${params.local_ref_files_path}/T2T_chm13v2.0+bs_controls.fa",
    ucsc_cpg_islands_gtf:     "${params.local_ref_files_path}/t2t2_ucsc_cpg_islands.gtf.gz",
    cpg_chr_lookup:           '$10,$10',  // T2T: col10 used for both source and dest
    refseq_gff_url:           'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/110/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz',
    ncbi_assembly_report_url: 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/110/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt',
    refseq_chr_lookup:        '$7,$10',   // NC_ accession -> chr1-style name
    epd_promoter_bed_url:     'https://epd.expasy.org/ftp/epdnew/human/006/Hs_EPDnew_006_hg38.bed',
    old_new_chain_url:        'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz',
]
def preset = params.mouse ? mouse_preset : (params.human_t2t2 ? t2t_preset : [:])

// Resolve effective values: CLI param > preset default > null.
// Collected into a single `ref` map so references throughout the file are
// unambiguous (ref.genome, ref.cpg_chr_lookup, …) and won't shadow
// bash variables inside process script blocks.
def ref = [
    genome:                   params.genome                   ?: preset.genome,
    ucsc_cpg_islands_gtf:     params.ucsc_cpg_islands_gtf     ?: preset.ucsc_cpg_islands_gtf,
    cpg_chr_lookup:           params.cpg_chr_lookup           ?: preset.cpg_chr_lookup,
    refseq_gff_url:           params.refseq_gff_url           ?: preset.refseq_gff_url,
    ncbi_assembly_report_url: params.ncbi_assembly_report_url ?: preset.ncbi_assembly_report_url,
    refseq_chr_lookup:        params.refseq_chr_lookup        ?: preset.refseq_chr_lookup,
    epd_promoter_bed_url:     params.epd_promoter_bed_url     ?: preset.epd_promoter_bed_url,
    old_new_chain_url:        params.old_new_chain_url        ?: preset.old_new_chain_url,
]

// Validate all required genome values are set
def missing = ref.findAll { k, v -> v == null }.keySet().toList()
if (missing) {
    error """\
        Missing required genome parameters: ${missing.join(', ')}
        Use --mouse or --human_t2t2 for preset assemblies, or provide these params directly for a custom assembly.
        """.stripIndent()
}


process fetch_chain_file {
    conda "curl"

    input:
        val url

    output:
        path('*.chain.gz')

    script:
    """
    curl -fsSL -o chain_file.chain.gz "${url}"
    """
}

process fetch_refseq_assembly_report {
    conda "curl"

    input:
        val url

    output:
        path('*.txt')

    script:
    """
    curl -fsSl "${url}" > assembly_report.txt
    """
}

process clean_epd_gtf {
    conda "curl ucsc-bedtogenepred ucsc-genepredtogtf crossmap"

    input:
        val(url)
        path(assembly_report)
        path(chain_file)

    output:
        path('epd_promoters.gtf')

    script:
    """

    curl -fsSL "${url}" \
      | tr ' ' '\t' \
      | bedToGenePred /dev/stdin /dev/stdout \
      | genePredToGtf file /dev/stdin /dev/stdout \
      > epd_promoters_oldref.gtf && \
      CrossMap gff ${chain_file} epd_promoters_oldref.gtf epd_promoters_newref.gtf

    awk -v OFS='\\t' -v FS='\\t' 'NR==FNR {dict[\$1]=\$2; next} {\$1=dict[\$1]; print}' \
      <(grep -v '^#' ${assembly_report} | awk -v OFS='\\t' -v FS='\\t' '{print ${ref.cpg_chr_lookup}}' | tr -d '\\r')  \
      <(zcat -f epd_promoters_newref.gtf | grep -v '^#') > epd_promoters.gtf
    """
}

process epd_promoter_counts{
    conda "subread=2.0.0"
    cpus 16

    input:
        path(gtf)
        path('*')

    output:
        path 'epd_promoter_counts.tsv'

    script:
    """
    featureCounts --primary ${feature_count_dup_option} -Q 10 -M -f -o -O --fraction -p -P -B -C \
        -a ${gtf} \
        --tmpDir ${params.tmp_dir} \
        -T ${task.cpus} \
        -o epd_promoter_counts.tsv *.bam
    """
}

process clean_cpg_islands_gtf {
    conda "gawk gzip"

    input:
        path ucsc_cpg_gtf
        path(assembly_report)

    output:
        path('cpg_islands.uniqname.gtf')

    script:
    """
    awk -v OFS='\\t' -v FS='\\t' 'NR==FNR {dict[\$1]=\$2; next} {\$1=dict[\$1]; print}' \
      <(grep -v '^#' ${assembly_report} | awk -v OFS='\\t' -v FS='\\t' '{print ${ref.cpg_chr_lookup}}' | tr -d '\\r')  \
      <(zcat -f ${ucsc_cpg_gtf} | grep -v '^#') \
      |  awk -v FS='\t' -v OFS='\t' '{print \$1,\$1":"\$4"-"\$5,\$3,\$4,\$5,\$6,\$7,\$8,\$9}' \
      > cpg_islands.uniqname.gtf
    """
}

process cpg_island_counts{
    conda "subread=2.0.0"

    cpus 16

    input:
        path gtf
        path('*')

    output:
        path 'cpg_island_counts.tsv'

    script:
    """
    featureCounts --primary ${feature_count_dup_option} -Q 10 -M -f -o -O --fraction -p -P -B -C \
        -a ${gtf} \
        --tmpDir ${params.tmp_dir} \
        -T ${task.cpus} \
        -o cpg_island_counts.tsv *.bam
    """
}

process refseq_feature_download {
    conda "curl"

    input:
        val url

    output:
        path('*.gff')

    script:
    """
    curl -fsSl ${url} | zcat -f > refseq.gff
    """
}

process refseq_feature_gffs {
    tag {feature}
    conda "subread=2.0.0 bedtools=2.29.2"

    input:
        path(gff)
        path(assembly_report)
        val(feature)

    output:
        tuple val(feature), path("${feature}.gff")
        tuple val(feature), path('*_flat.saf')

    script:
    """
    awk -v OFS='\\t' -v FS='\\t' 'NR==FNR {dict[\$1]=\$2; next} {\$1=dict[\$1]; print}' \
     <(grep -v '^#' ${assembly_report} \
       | awk -v OFS='\\t' -v FS='\\t' '{print ${ref.refseq_chr_lookup}}' \
       | tr -d '\\r')  \
     <(zcat -f ${gff} | grep -v '^#') \
    | grep "GeneID:" \
    | grep -P -v "_alt\\t" \
    | grep -P -v "^na\\t" \
    | sed -r 's/;Dbxref(=[^;]*)GeneID:([^,;]+)([;,])/;gene_id=\\2;Dbxref\\1GeneID:\\2\\3/' \
    | awk  -v OFS='\\t' -v FS='\\t' \
       '(\$3=="exon") && (index(\$9,"gbkey=mRNA") > 0) && (index(\$9,"-1;Parent") > 0) \
          { print(\$1,\$2,"mRNAexon1",\$4,\$5,\$6,\$7,\$8,\$9); next }
        (\$3=="exon") && (index(\$9,"gbkey=mRNA") > 0) \
          { print(\$1,\$2,"mRNAexon",\$4,\$5,\$6,\$7,\$8,\$9); next }
        { print }
       ' > name_converted.gff

    # exons overlap, we want only the longest to avoid 0 cov exons from featureCounts
    flattenGTF -a name_converted.gff -o flat_name_converted.saf -t ${feature}

    # need to switch to bed for intersection later
    tail -n +2 flat_name_converted.saf \
     | awk -v OFS='\\t' -v FS='\\t' '{print \$2,\$3-1,\$4,\$1,"-",\$5}' \
     | bedtools sort -faidx ${ref.genome}.fai -i /dev/stdin > ${feature}_flat.bed

    # filters by feature type
    awk -v type=${feature} -v OFS='\\t' -v FS='\\t' '(\$3==type) { print}' name_converted.gff > ${feature}.gff

    # only include those entries that intersect with the desired feature type, back to SAF format
    echo "GeneID\tChr\tStart\tEnd\tStrand" > ${feature}_flat.saf
    bedtools intersect -a ${feature}_flat.bed -b ${feature}.gff  -u \
    | awk -v OFS='\\t' -v FS='\\t' '{print \$4,\$1,\$2+1,\$3,\$6}' >> ${feature}_flat.saf
    """
}

process feature_methylation {
    tag {feature}
    conda "bedtools=2.29.2 htslib=1.9"
    memory '32 GB'

    input:
        tuple val(feature), val(feature_gff), val(sample_id), val(context), path(methylkit)

    output:
        tuple val(sample_id), val(context), path('*_methylation.tsv')

    script:
    // gff 9 columns:  CM000994.3      cmsearch        exon    3172239 3172348 .       +       .       ID=exon-X...
    // input:
    // chrBase chr     base    strand  coverage        freqC   freqT
    //CM000994.3.3050095      CM000994.3      3050095 F       4         0.00  100.00
    // "bed" (with additional columns): chr, base0, base1, methylation fraction, coverage
    // coverage >= 5
    // chr, start, end, chr:base0-base1, methylation proportion -> groupby column 4, mean of column 5
    """
    bedtools intersect -nonamecheck \
    -wa -wb -loj \
    -a ${feature_gff} -b <(zcat ${methylkit} | awk -v FS='\\t' -v OFS='\\t' 'NR>1 {print \$2, \$3-1, \$3, \$6/100, \$5}') \
    | awk -v FS='\\t' -v OFS='\\t' '\$14>=5 {print \$10,\$11,\$12,\$1":"\$4-1"-"\$5, \$13}' \
    | sort -k4,4 | bedtools groupby -g 4 -o mean -c 5 \
    > ${feature}_methylation.tsv
    """
}

process combine_feature_methylation {
    publishDir "${params.output_dir}/features/${sample_id}/${context}", mode: 'copy'

    input:
        tuple val(sample_id), val(context), path(methylation_files)

    output:
        path('*combined_methylation.tsv')

    script:
    """
    echo 'Feature\tLocus\tMeth' > ${sample_id}_${context}_combined_methylation.tsv
    for f in ${methylation_files} ; do
        filebase=\$(basename "\${f}" _methylation.tsv)
        # Use awk to prepend the feature name — avoids SIGPIPE from 'yes | head'
        awk -v name="\${filebase}" -v OFS='\\t' 'NF && !/^[[:space:]]*\$/ && !/^#/ {print name, \$0}' "\${f}" \
            >> ${sample_id}_${context}_combined_methylation.tsv
    done
    """
}

process refseq_feature_counts {
    conda "subread=2.0.0"
    publishDir "${params.output_dir}", mode: 'copy'
    cpus 16

    input:
        tuple val(feature), path(feature_saf)
        path bam_files

    output:
        path '*_counts.tsv'

    script:
    """
    featureCounts --primary ${feature_count_dup_option} -Q 10 -M -f -O --fraction -p -P -B -C \
    -a ${feature_saf} -F SAF\
    -t ${feature} \
    -g 'ID' \
    --tmpDir ${params.tmp_dir} \
    -T ${task.cpus} \
    -o ${feature}_counts.tsv *.bam
    """
}

process feature_depth_bokeh {
    cpus 1
    memory '32 GB'
    conda "bokeh polars numpy"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        path(combined_counts)

    output:
        path('feature_depth_bokeh.html')

    script:
    """
    feature_depth_bokeh.py ${combined_counts} -o feature_depth_bokeh.html
    """
}

process combine_counts {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        path(feature_counts)
        path(cpg_counts)
        path(epd_counts)

   output:
       path('combined_feature_counts.tsv')

   script:
   """
   echo -n 'File\t' > combined_feature_counts.tsv
   # Read the header directly with awk to avoid SIGPIPE from grep|head
   awk '!/^[[:space:]]*\$/ && !/^#/ {print; exit}' ${cpg_counts} >> combined_feature_counts.tsv

   for f in ${feature_counts} ${cpg_counts} ${epd_counts}; do
       filebase=\$(basename "\${f}" _counts.tsv)
       # Count data lines without piping through head/wc on process substitutions
       lines=\$(grep -ve '^\\s*\$' -e '^#' "\${f}" | awk 'END{print NR}')
       # Use awk to repeat filebase n times instead of 'yes | head' (avoids SIGPIPE)
       paste <(awk -v n="\${lines}" -v name="\${filebase}" 'BEGIN{for(i=0;i<n;i++) print name}') \\
             <(grep -ve '^\\s*\$' -e '^#' "\${f}") | tail -n +2 >> combined_feature_counts.tsv
   done
   """
}


workflow {

    // Define channels
    methylkits = Channel.fromPath(params.mk_files)
        .map{ it -> tuple(it.baseName.split(".methylKit")[0], it.baseName.split(".methylKit")[0].split("_")[-1], it)}
        .filter{ params.context.tokenize(',').contains(it[1]) }

    bams = Channel.fromFilePairs(params.bam_files_glob, checkIfExists: true)

    // Run processes
    chain_file        = fetch_chain_file(Channel.value(ref.old_new_chain_url))
    ncbi_assembly_report = fetch_refseq_assembly_report(Channel.value(ref.ncbi_assembly_report_url))
    epd_promoters_gtf = clean_epd_gtf(Channel.value(ref.epd_promoter_bed_url), ncbi_assembly_report, chain_file)
    epd_counts        = epd_promoter_counts(epd_promoters_gtf, bams.map{ it -> it[1] }.collect())
    cpg_islands_gtf   = clean_cpg_islands_gtf(Channel.fromPath(ref.ucsc_cpg_islands_gtf, checkIfExists: true).first(), ncbi_assembly_report)
    cpg_counts        = cpg_island_counts(cpg_islands_gtf, bams.map{ it -> it[1] }.collect())
    refseq_gff        = refseq_feature_download(Channel.value(ref.refseq_gff_url))

    // Define feature types
    refseq_feature_types = Channel.of('exon', 'CDS', 'gene', 'mRNA', 'mRNAexon', 'mRNAexon1')

    // Process RefSeq features
    refseq_outputs = refseq_feature_gffs(refseq_gff, ncbi_assembly_report, refseq_feature_types)
    feature_gff_for_meth = refseq_outputs[0]
    feature_saf_for_counts = refseq_outputs[1]

    // Feature methylation
    feature_methylation_ch = feature_gff_for_meth
        .concat(cpg_islands_gtf.map{ it -> tuple("cpg_islands", it) })
        .concat(epd_promoters_gtf.map{ it -> tuple("epd_promoters", it) })
        .combine(methylkits)

    feature_methylation_out = feature_methylation(feature_methylation_ch)
    combined_methylation_out = combine_feature_methylation(feature_methylation_out.groupTuple(by: [0,1]))

    // Feature counts
    all_bam_files = bams.map{ it -> it[1] }.flatten().collect()

    feature_counts = refseq_feature_counts(feature_saf_for_counts, all_bam_files)

    // Combine all counts
    combined_counts = combine_counts(feature_counts.collect(), cpg_counts, epd_counts)

    // Bokeh dot plot of depth per library, per feature category
    feature_depth_bokeh(combined_counts)
}
