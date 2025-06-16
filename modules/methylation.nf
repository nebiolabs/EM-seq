process methylDackel_mbias {
    label 'medium_cpu'
    errorStrategy 'retry'
    tag "${library}"
    conda "bioconda::methyldackel=0.6.1 bioconda::samtools=1.21 conda-forge::pigz=2.8 conda-forge::sed=4.9"
    publishDir "${params.outputDir}/methylDackelExtracts/mbias"

    input:
        tuple val(library), path(md_bam), path(md_bai)
        tuple path(genome_fa), path(genome_fai)

    output:
        path('*.svg'), emit: mbias_output_svg
        path('*.tsv'), emit: mbias_output_tsv
        tuple val(library), path('*.tsv'), emit: for_agg

    script:
    """
    echo -e "chr\tcontext\tstrand\tRead\tPosition\tnMethylated\tnUnmethylated\tnMethylated(+dups)\tnUnmethylated(+dups)" > ${library}_${barcodes}_combined_mbias.tsv
    chrs=(`samtools view -H "${md_bam}" | grep @SQ | cut -f 2 | sed 's/SN://'| grep -v _random | grep -v chrUn | sed 's/|/\\|/'`)

    for chr in \${chrs[*]}; do
        for context in CHH CHG CpG; do
            arg=''
            if [ "\$context" = 'CHH' ]; then
            arg='--CHH --noCpG'
            elif [ "\$context" = 'CHG' ]; then
            arg='--CHG --noCpG'
            fi
            # need two calls to add columns containing the counts without filtering duplicate reads (for rrEM-seq where start/end is constrained)
            # not sure why we need both --keepDupes and -F, probably a bug in mbias
            join -t \$'\t' -j1 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -a 1 -e 0 \
            <( \
                MethylDackel mbias --noSVG \$arg -@ ${task.cpus} -r \$chr ${genome_fa} "${md_bam}" | \
                tail -n +2 | awk '{print \$1"-"\$2"-"\$3"\t"\$0}' | sort -k 1b,1
            ) \
            <( \
                MethylDackel mbias --noSVG --keepDupes -F 2816 \$arg -@ ${task.cpus} -r \$chr ${genome_fa} "${md_bam}" | \
                tail -n +2 | awk '{print \$1"-"\$2"-"\$3"\t"\$0}' | sort -k 1b,1
            ) \
            | sed "s/^/\${chr}\t\${context}\t/" \
            >> ${library}_combined_mbias.tsv
        done
    done
    # makes the svg files for trimming checks
    MethylDackel mbias -@ ${task.cpus} --noCpG --CHH --CHG -r \${chrs[0]} ${genome_fa} "${md_bam}" ${library}_chn
    for f in *chn*.svg; do sed -i "s/Strand<\\/text>/Strand \$f \${chrs[0]} CHN <\\/text>/" \$f; done;

    MethylDackel mbias -@ ${task.cpus} -r \${chrs[0]} ${genome_fa} "${md_bam}" ${library}_cpg
    for f in *cpg*.svg; do sed -i "s/Strand<\\/text>/Strand \$f \${chrs[0]} CpG<\\/text>/" \$f; done;
    """
}


process methylDackel_extract {
    label 'high_cpu'
    tag "${library}"
    publishDir "${params.outputDir}/methylDackelExtracts", mode: 'copy'
    conda "bioconda::methyldackel=0.6.1 bioconda::samtools=1.21 conda-forge::pigz=2.8"

    input:
        tuple val(library), path(md_bam), path(md_bai)
        tuple path(genome_fa), path(genome_fai)

    output:
        tuple val(library), path('*CHG.methylKit.gz'), path('*CHH.methylKit.gz'),path('*CpG.methylKit.gz'), emit: extract_output 

    script:
    """
    MethylDackel extract --methylKit -q 20 --nOT 0,0,0,5 --nOB 0,0,5,0 -@ ${task.cpus} \
        --CHH --CHG -o ${library} ${genome_fa} "${md_bam}" 
    pigz -p ${task.cpus} *.methylKit 
    """
}

process convert_methylkit_to_bed {
    label 'low_cpu'
    tag "${library}"
    conda "conda-forge::pigz=2.8 conda-forge::gawk=5.3.1 conda-forge::sed=4.9"

    input:
        tuple val(library), path(methylkit_CHH_gz),path(methylkit_CHG_gz),path(methylkit_CpG_gz),path(genome_fa),path(genome_fai)

    output:
        tuple val(library), path('*.bed'), emit: methylkit_bed

    script:
    """
    methylkit_basename=\$(basename "${methylkit_CpG_gz}" _CpG.methylKit.gz)
    
    # create sed scripts to interconvert chr names from genome index with line numbers for efficient sorting
    awk '{print "s/\\t"\$1"\\t/\\t"NR-1"\\t/"}' "${genome_fai}" > convert_chr_to_num.sed
    awk '{print "s/^"NR-1"\\t/"\$1"\\t/"}' "${genome_fai}" > revert_num_to_chr.sed 
    # merge initial methylkit files into a single BED format
    # Convert methylKit to BED format (0-based coordinates)
    # methylKit format: chr.base\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT
    # Output format: chr\tstart\tend\tname\tmethylation\tstrand
    sort -m -k2,2n -k3,3n \
            <(pigz -d -c "${methylkit_CHH_gz}" | sed 's/\$/\\tCHH/' | sed -f convert_chr_to_num.sed) \
            <(pigz -d -c "${methylkit_CHG_gz}" | sed 's/\$/\\tCHG/' | sed -f convert_chr_to_num.sed) \
            <(pigz -d -c "${methylkit_CpG_gz}" | sed 's/\$/\\tCpG/' | sed -f convert_chr_to_num.sed) \
    | awk -v OFS="\\t" 'NR>3{print \$2,\$3-1,\$3,\$8,\$6,(\$4 == "F") ? "+" : "-"}' \
    | sed -f revert_num_to_chr.sed \
    > "\${methylkit_basename}.bed"
    """
}
