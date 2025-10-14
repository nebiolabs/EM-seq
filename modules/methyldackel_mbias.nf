process methylDackel_mbias {
    label 'medium_cpu'
    errorStrategy 'retry'
    tag "${library}"
    conda "bioconda::methyldackel=0.6.1 bioconda::samtools=1.21 conda-forge::pigz=2.8 conda-forge::sed=4.9"
    publishDir "${params.outputDir}/methylDackelExtracts/mbias", mode: 'copy'

    input:
        tuple val(library), path(md_bam), path(md_bai)
        val(genome_fa)
        val(genome_fai)

    output:
        path('*.svg'), emit: mbias_output_svg
        path('*.tsv'), emit: mbias_output_tsv
        tuple val(library), path("${library}.combined_mbias.tsv"), emit: for_agg
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions
        tuple val("${task.process}"), val('methyldackel'), eval('MethylDackel --version 2>&1 | cut -f 2 -d ":"'), topic: versions

    script:
    """
    echo -e "chr\tcontext\tstrand\tRead\tPosition\tnMethylated\tnUnmethylated\tnMethylated(+dups)\tnUnmethylated(+dups)" > ${library}.combined_mbias.tsv
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
            >> ${library}.combined_mbias.tsv
        done
    done
    # makes the svg files for trimming checks
    MethylDackel mbias -@ ${task.cpus} --noCpG --CHH --CHG -r \${chrs[0]} ${genome_fa} "${md_bam}" ${library}_chn
    for f in *chn*.svg; do sed -i.bak "s/Strand<\\/text>/Strand \$f \${chrs[0]} CHN <\\/text>/" \$f; done;

    MethylDackel mbias -@ ${task.cpus} -r \${chrs[0]} "${genome_fa}" "${md_bam}" ${library}_cpg
    for f in *cpg*.svg; do sed -i.bak "s/Strand<\\/text>/Strand \$f \${chrs[0]} CpG<\\/text>/" \$f; done;
    """
}
