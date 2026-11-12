nextflow.enable.dsl = 2

params.mk_files = '*.methylKit.gz'
params.bam_files_glob = '*.md.{bam,bam.bai}'
params.tmp_dir = '/tmp/'
params.output_dir = 'feature_cov.output'
params.count_dup_reads = false
params.mouse = false
params.human_t2t2 = false
params.context = 'CpG' // methylation context(s) to analyse; comma-separated for multiple: 'CpG,CHG,CHH'

// Path to locally cached reference files used by the --mouse and --human_t2t2 shortcuts.
// Override with --local_ref_files_path for a different location.
params.local_ref_files_path = "${System.getenv('HOME')}/nebnext_projects/em-seq/em-seq_ref_files"

// Genome-specific params declared once with null defaults.
// Set these on the CLI for a custom assembly; --mouse and --human_t2t2 pre-populate them.
// cpg_chr_lookup / refseq_chr_lookup are awk column specs (src_col,dest_col) translating
// assembly-report accession IDs to the chromosome naming used in BAMs and GTF/SAF files.
params.genome = null
params.ucsc_cpg_islands_gtf = null
params.cpg_chr_lookup = null // e.g. '$10,$5' (mouse) or '$10,$10' (T2T)
params.refseq_gff_url = null
params.ncbi_assembly_report_url = null
params.refseq_chr_lookup = null // e.g. '$7,$5' (mouse) or '$7,$10' (T2T)
params.epd_promoter_bed_url = null
params.old_new_chain_url = null

params.feature_count_dup_option = params.count_dup_reads ? '' : '--ignoreDup'

params.invert_pattern = null // sample-id regex for chemistries reading modified C as T (TAPS, 5-base)
params.meth_min_coverage = 5 // per-site coverage
params.meth_min_cytosines = 5 // covered sites required per feature
params.repeatmasker_bed = null // UCSC rmsk .bb or table (not RepeatMasker .out), path or URL; repeats skipped if unset
params.repeat_classes = 'SINE/Alu:270,LINE/L1:5400' // class/family:min ALIGNED bp (~90% of consensus, keeps mostly full length subset)

// Resolve effective values: CLI param > preset default > null.
// Collected into a single `params.ref` map so references throughout the file are
// unambiguous (params.ref.genome, params.ref.cpg_chr_lookup, …) and won't shadow
// bash variables inside process script blocks.
def resolveGenomeRefs() {
    // Preset maps for known assemblies. Using local maps rather than re-assigning params
    // avoids Nextflow lint warnings about params being defined multiple times.
    def mouse_preset = [
        genome: "${params.local_ref_files_path}/grcm39+meth_controls.fa",
        ucsc_cpg_islands_gtf: "${params.local_ref_files_path}/grcm39_cpg_islands.gtf",
        cpg_chr_lookup: '$10,$5', // assembly report col10 (chrN) -> col5 (GenBank accession)
        refseq_gff_url: 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/annotation_releases/109/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz',
        ncbi_assembly_report_url: 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/annotation_releases/109/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_assembly_report.txt',
        refseq_chr_lookup: '$7,$5', // NC_ accession -> GenBank chr name
        epd_promoter_bed_url: 'https://epd.expasy.org/ftp/epdnew/M_musculus/003/Mm_EPDnew_003_mm10.bed',
        old_new_chain_url: 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz',
    ]
    def t2t_preset = [
        genome: "${params.local_ref_files_path}/T2T_chm13v2.0+bs_controls.fa",
        ucsc_cpg_islands_gtf: "${params.local_ref_files_path}/t2t2_ucsc_cpg_islands.gtf.gz",
        cpg_chr_lookup: '$10,$10', // T2T: col10 used for both source and dest
        refseq_gff_url: 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/110/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz',
        ncbi_assembly_report_url: 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/110/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt',
        refseq_chr_lookup: '$7,$10', // NC_ accession -> chr1-style name
        epd_promoter_bed_url: 'https://epd.expasy.org/ftp/epdnew/human/006/Hs_EPDnew_006_hg38.bed',
        old_new_chain_url: 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz',
        repeatmasker_bed: 'https://hgdownload.soe.ucsc.edu/gbdb/hs1/t2tRepeatMasker/chm13v2.0_rmsk.bb',
    ]
    def preset = params.mouse ? mouse_preset : (params.human_t2t2 ? t2t_preset : [:])

    def ref = [
        genome: params.genome ?: preset.genome,
        ucsc_cpg_islands_gtf: params.ucsc_cpg_islands_gtf ?: preset.ucsc_cpg_islands_gtf,
        cpg_chr_lookup: params.cpg_chr_lookup ?: preset.cpg_chr_lookup,
        refseq_gff_url: params.refseq_gff_url ?: preset.refseq_gff_url,
        ncbi_assembly_report_url: params.ncbi_assembly_report_url ?: preset.ncbi_assembly_report_url,
        refseq_chr_lookup: params.refseq_chr_lookup ?: preset.refseq_chr_lookup,
        epd_promoter_bed_url: params.epd_promoter_bed_url ?: preset.epd_promoter_bed_url,
        old_new_chain_url: params.old_new_chain_url ?: preset.old_new_chain_url,
        repeatmasker_bed: params.repeatmasker_bed ?: preset.repeatmasker_bed,
    ]

    // Validate all required genome values are set (repeatmasker_bed is optional: repeats are skipped without it)
    def missing = ref.findAll { k, v -> v == null && k != 'repeatmasker_bed' }.keySet().toList()
    if (missing) {
        error(
            """\
            Missing required genome parameters: ${missing.join(', ')}
            Use --mouse or --human_t2t2 for preset assemblies, or provide these params directly for a custom assembly.
            """.stripIndent()
        )
    }
    return ref
}

params.ref = resolveGenomeRefs()


process fetch_chain_file {
    conda "curl"

    input:
    val url

    output:
    path '*.chain.gz'

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
    path '*.txt'

    script:
    """
    curl -fsSl "${url}" > assembly_report.txt
    """
}

process clean_epd_gtf {
    conda "curl ucsc-bedtogenepred ucsc-genepredtogtf crossmap"

    input:
    val url
    path assembly_report
    path chain_file

    output:
    path 'epd_promoters.gtf'

    script:
    """

    curl -fsSL "${url}" \
      | tr ' ' '\t' \
      | bedToGenePred /dev/stdin /dev/stdout \
      | genePredToGtf file /dev/stdin /dev/stdout \
      > epd_promoters_oldref.gtf && \
      CrossMap gff ${chain_file} epd_promoters_oldref.gtf epd_promoters_newref.gtf

    awk -v OFS='\\t' -v FS='\\t' 'NR==FNR {dict[\$1]=\$2; next} {\$1=dict[\$1]; print}' \
      <(grep -v '^#' ${assembly_report} | awk -v OFS='\\t' -v FS='\\t' '{print ${params.ref.cpg_chr_lookup}}' | tr -d '\\r')  \
      <(zcat -f epd_promoters_newref.gtf | grep -v '^#') > epd_promoters.gtf
    """
}

process epd_promoter_counts {
    conda "subread=2.0.0"
    cpus 16

    input:
    path gtf
    path '*'

    output:
    path 'epd_promoter_counts.tsv'

    script:
    """
    featureCounts --primary ${params.feature_count_dup_option} -Q 10 -M -f -o -O --fraction -p -P -B -C \
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
    path assembly_report

    output:
    path 'cpg_islands.uniqname.gtf'

    script:
    """
    awk -v OFS='\\t' -v FS='\\t' 'NR==FNR {dict[\$1]=\$2; next} {\$1=dict[\$1]; print}' \
      <(grep -v '^#' ${assembly_report} | awk -v OFS='\\t' -v FS='\\t' '{print ${params.ref.cpg_chr_lookup}}' | tr -d '\\r')  \
      <(zcat -f ${ucsc_cpg_gtf} | grep -v '^#') \
      |  awk -v FS='\t' -v OFS='\t' '{print \$1,\$1":"\$4"-"\$5,\$3,\$4,\$5,\$6,\$7,\$8,\$9}' \
      > cpg_islands.uniqname.gtf
    """
}

process cpg_island_counts {
    conda "subread=2.0.0"

    cpus 16

    input:
    path gtf
    path '*'

    output:
    path 'cpg_island_counts.tsv'

    script:
    """
    featureCounts --primary ${params.feature_count_dup_option} -Q 10 -M -f -o -O --fraction -p -P -B -C \
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
    path '*.gff'

    script:
    """
    curl -fsSl ${url} | zcat -f > refseq.gff
    """
}

process refseq_feature_gffs {
    tag { feature }
    conda "subread=2.0.0 bedtools=2.29.2"

    input:
    path gff
    path assembly_report
    val feature

    output:
    tuple val(feature), path("${feature}.gff")
    tuple val(feature), path('*_flat.saf')
    tuple val(feature), path('*_flat.bed')

    script:
    """
    awk -v OFS='\\t' -v FS='\\t' 'NR==FNR {dict[\$1]=\$2; next} {\$1=dict[\$1]; print}' \
     <(grep -v '^#' ${assembly_report} \
       | awk -v OFS='\\t' -v FS='\\t' '{print ${params.ref.refseq_chr_lookup}}' \
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
     | bedtools sort -faidx ${params.ref.genome}.fai -i /dev/stdin > ${feature}_flat.bed

    # filters by feature type
    awk -v type=${feature} -v OFS='\\t' -v FS='\\t' '(\$3==type) { print}' name_converted.gff > ${feature}.gff

    # only include those entries that intersect with the desired feature type, back to SAF format
    echo "GeneID\tChr\tStart\tEnd\tStrand" > ${feature}_flat.saf
    bedtools intersect -a ${feature}_flat.bed -b ${feature}.gff  -u \
    | awk -v OFS='\\t' -v FS='\\t' '{print \$4,\$1,\$2+1,\$3,\$6}' >> ${feature}_flat.saf
    """
}

process fetch_repeatmasker {
    conda "curl ucsc-bigbedtobed"

    input:
    val src

    output:
    path 'repeatmasker.tsv.gz'

    script:
    // src: URL or absolute local path; bigBed (.bb) or text table (optionally gzipped)
    """
    if [[ "${src}" == *://* ]]; then
        curl -fsSL -o rmsk.src "${src}"
    else
        ln -s "${src}" rmsk.src
    fi
    if [[ "${src}" == *.bb ]]; then
        bigBedToBed rmsk.src /dev/stdout | gzip -c > repeatmasker.tsv.gz
    else
        zcat -f rmsk.src | gzip -c > repeatmasker.tsv.gz
    fi
    """
}

process repeat_feature_gffs {
    tag { feature }
    conda "bedtools=2.29.2"

    input:
    tuple val(feature), val(rclass), val(minaln), path(rmsk)

    output:
    tuple val(feature), path("${feature}_flat.saf")
    tuple val(feature), path("${feature}_flat.bed")

    script:
    """
    zcat -f ${rmsk} \
    | awk -v OFS='\t' -v cls='${rclass}' -v minaln=${minaln} -v feat='${feature}' \
        'BEGIN{FS="\t"} \$1 !~ /^#/ {
            h = index(\$4, "#"); if (h == 0 || substr(\$4, h+1) != cls) next
            # aligned bases only: blockStart == -1 marks an unaligned block
            ns = split(\$11, bs, ","); split(\$12, bt, ",")
            aln = 0
            for (i = 1; i <= ns; i++)
                if (bs[i] != "" && bt[i] != "" && bt[i]+0 >= 0) aln += bs[i]
            if (aln >= minaln)
                print \$1, \$2, \$3, feat"_"\$1"_"\$2"_"\$3, ".", \$6
         }' \
    | bedtools sort -faidx ${params.ref.genome}.fai -i /dev/stdin \
    > ${feature}_flat.bed

    echo "GeneID\tChr\tStart\tEnd\tStrand" > ${feature}_flat.saf
    awk -v OFS='\t' '{print \$4,\$1,\$2+1,\$3,\$6}' ${feature}_flat.bed >> ${feature}_flat.saf

    echo "${feature}: \$(wc -l < ${feature}_flat.bed) intervals" >&2
    """
}

process feature_methylation {
    tag "${sample_id} ${context}"
    conda "gawk"
    memory '48 GB'

    input:
    tuple val(names), path(gffs), val(fmts), val(sample_id), val(context), path(methylkit)

    output:
    tuple val(sample_id), val(context), path('*_methylation.tsv')

    script:
    // One streaming pass over the methylKit file for all feature classes
    // (replaces per-class bedtools intersect, which could not handle CHG/CHH).
    def invert = params.invert_pattern && (sample_id ==~ params.invert_pattern)
    def manifest = [names, gffs.collect { it.name }, fmts].transpose().collect { n, f, t -> "${n}\t${f}\t${t}" }.join('\n')
    """
    cat > manifest.tsv <<'MANIFEST'
${manifest}
MANIFEST

    # every class must produce a file, even an empty one, or the output glob
    # silently loses a class
    cut -f1 manifest.tsv | while read -r c; do : > "\${c}_methylation.tsv"; done

    zcat ${methylkit} | gawk -v FS='\t' -v OFS='\t' \
        -v MANIFEST=manifest.tsv -v BUCKET=1000 \
        -v MINCOV=${params.meth_min_coverage} -v MINC=${params.meth_min_cytosines} \
        -v INV=${invert ? 1 : 0} '
      BEGIN {
        while ((getline line < MANIFEST) > 0) {
          split(line, a, "\\t"); cls = a[1]; gff = a[2]; fmt = a[3]
          while ((getline g < gff) > 0) {
            if (substr(g, 1, 1) == "#") continue
            if (split(g, f, "\\t") < 5) continue
            ne++
            # BED is 0-based half-open; GFF/GTF is 1-based inclusive
            st = (fmt == "bed") ? f[2] + 0 : f[4] - 1
            en = (fmt == "bed") ? f[3] + 0 : f[5] + 0
            CL[ne] = cls; S[ne] = st; E[ne] = en
            # printed locus is 1-based to match combined_feature_counts.tsv
            ID[ne] = f[1] ":" (st + 1) "-" en
            lo = int(st / BUCKET); hi = int((en - 1) / BUCKET)
            for (k = lo; k <= hi; k++) idx[f[1] ":" k] = idx[f[1] ":" k] ne ","
          }
          close(gff)
          printf("  indexed %-14s cumulative %d intervals\\n", cls, ne) > "/dev/stderr"
        }
        close(MANIFEST)
      }
      FNR > 1 && \$5 + 0 >= MINCOV {
        p = \$3 - 1
        key = \$2 ":" int(p / BUCKET)
        if (!(key in idx)) next
        # +0: MethylDackel pads frequencies with leading spaces
        m =\$6 + 0
        if (INV == 1) m = 100 - m
        n = split(idx[key], c, ",")
        for (i = 1; i <= n; i++) {
          if (c[i] == "") continue
          j = c[i] + 0
          if (p >= S[j] && p < E[j]) { sum[j] += m; cnt[j]++ }
        }
      }
      END {
        for (j = 1; j <= ne; j++)
          if (cnt[j] >= MINC)
            printf("%s\\t%.4f\\t%d\\n", ID[j], sum[j] / cnt[j], cnt[j]) \
              > (CL[j] "_methylation.tsv")
      }'

    for f in *_methylation.tsv; do
        echo "${sample_id} ${context} invert=${invert} \$f: \$(wc -l < "\$f") features" >&2
    done
    """
}

process combine_feature_methylation {
    publishDir { "${params.output_dir}/features/${sample_id}/${context}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(context), path(methylation_files)

    output:
    path '*combined_methylation.tsv'

    script:
    """
    echo 'Feature\tLocus\tMeth\tnCytosines' > ${sample_id}_${context}_combined_methylation.tsv
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
    featureCounts --primary ${params.feature_count_dup_option} -Q 10 -M -f -O --fraction -p -P -B -C \
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
    path combined_counts

    output:
    path 'feature_depth_bokeh.html'

    script:
    """
    feature_depth_bokeh.py ${combined_counts} -o feature_depth_bokeh.html
    """
}

process combine_counts {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path feature_counts
    path cpg_counts
    path epd_counts

    output:
    path 'combined_feature_counts.tsv'

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

    methylkits = channel.fromPath(params.mk_files)
        .map { f ->
            def stem = f.baseName.replaceFirst(/\.methylKit.*$/, '')
            def m = (stem =~ /_(CpG|CHG|CHH)(?:\..*)?$/)
            def ctx = m ? m[0][1] : stem.split('_')[-1]
            def sid = stem.replaceFirst(/_(CpG|CHG|CHH)(?:\..*)?$/, '')
            tuple(sid, ctx, f)
        }
        .filter { _sample_id, context, _methylkit -> params.context.tokenize(',').contains(context) }

    bams = channel.fromFilePairs(params.bam_files_glob, checkIfExists: true)

    // Run processes
    chain_file = fetch_chain_file(channel.value(params.ref.old_new_chain_url))
    ncbi_assembly_report = fetch_refseq_assembly_report(channel.value(params.ref.ncbi_assembly_report_url))
    epd_promoters_gtf = clean_epd_gtf(channel.value(params.ref.epd_promoter_bed_url), ncbi_assembly_report, chain_file)
    epd_counts = epd_promoter_counts(epd_promoters_gtf, bams.map { it -> it[1] }.collect())
    cpg_islands_gtf = clean_cpg_islands_gtf(channel.fromPath(params.ref.ucsc_cpg_islands_gtf, checkIfExists: true).first(), ncbi_assembly_report)
    cpg_counts = cpg_island_counts(cpg_islands_gtf, bams.map { it -> it[1] }.collect())
    refseq_gff = refseq_feature_download(channel.value(params.ref.refseq_gff_url))

    // Define feature types
    refseq_feature_types = channel.of('exon', 'CDS', 'gene', 'mRNA', 'mRNAexon', 'mRNAexon1')

    // Process RefSeq features
    refseq_outputs = refseq_feature_gffs(refseq_gff, ncbi_assembly_report, refseq_feature_types)
    feature_saf_for_counts = refseq_outputs[1]

    // Repeats, fed into the same counts and methylation channels
    if (params.ref.repeatmasker_bed) {
        repeat_specs = channel.of(params.repeat_classes.tokenize(','))
            .flatten()
            .map { spec ->
                def f = spec.trim().tokenize(':')
                tuple(f[0].replaceAll('[^A-Za-z0-9]', '_'), f[0], (f.size() > 1 ? f[1] as Integer : 0)) // name, class/family, min aligned bases
            }
        // local paths are made absolute: the fetch task runs in its own work dir
        def rmsk_src = params.ref.repeatmasker_bed ==~ /^[a-z]+:\/\/.*/
            ? params.ref.repeatmasker_bed
            : file(params.ref.repeatmasker_bed, checkIfExists: true).toString()
        rmsk = fetch_repeatmasker(channel.value(rmsk_src))
        repeat_outputs = repeat_feature_gffs(repeat_specs.combine(rmsk))
        feature_saf_for_counts = feature_saf_for_counts.concat(repeat_outputs[0])
    }

    all_feature_intervals = refseq_outputs[2]
        .concat(params.ref.repeatmasker_bed ? repeat_outputs[1] : channel.empty())
        .map { name, f -> tuple(name, f, 'bed') }
        .concat(cpg_islands_gtf.map { it -> tuple("cpg_island", it, 'gff') })
        .concat(epd_promoters_gtf.map { it -> tuple("epd_promoter", it, 'gff') })

    feature_methylation_ch = all_feature_intervals
        .toList()
        .map { t -> tuple(t.collect { it[0] }, t.collect { it[1] }, t.collect { it[2] }) }
        .combine(methylkits)

    feature_methylation_out = feature_methylation(feature_methylation_ch)
    // no groupTuple: each task already emits every class
    combine_feature_methylation(feature_methylation_out)

    // Feature counts
    all_bam_files = bams.map { it -> it[1] }.flatten().collect()

    feature_counts = refseq_feature_counts(feature_saf_for_counts, all_bam_files)

    // Combine all counts
    combined_counts = combine_counts(feature_counts.collect(), cpg_counts, epd_counts)

    // Bokeh dot plot of depth per library, per feature category
    feature_depth_bokeh(combined_counts)
}
