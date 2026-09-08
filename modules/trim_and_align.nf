process trimAndAlign {
    label 'high_cpu'
    tag { library }
    conda "conda-forge::python=3.10 bioconda::bwameth=0.2.7 bioconda::mark-nonconverted-reads=1.2 bioconda::samtools=1.22 bioconda::bamslice=0.2.1 bioconda::fgumi=0.7.0"
    publishDir "${params.outputDir}/bwameth_align", mode: 'symlink'

    input:
        tuple val(library), path(bam), val(start_offset), val(end_offset)
        path(adapter_fasta)
        val(bwa_index)
        path(genome_fa)
        path(genome_fai)
        path(genome_dict)

    output:
        tuple val(library), path("${chunk}.aln.bam"), emit: bam_files
        tuple val(library), val(chunk), path("${chunk}.fastp.json"), emit: fastp_json
        tuple val(library), path("${chunk}.nonconverted_counts.tsv"), emit: nonconverted_counts
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions
        tuple val("${task.process}"), val('bwameth'), eval("bwameth.py --version | sed 's/^bwa-meth.py //'"), topic: versions
        tuple val("${task.process}"), val('python'), eval('python --version | sed \'s/^Python //\''), topic: versions
        tuple val("${task.process}"), val('bamslice'), eval('bamslice --version | cut -f 2 -d " "'), topic: versions
        tuple val("${task.process}"), val('fgumi'), eval('fgumi --version | cut -f 2 -d " "'), topic: versions
        tuple val("${task.process}"), val('bamadap'), eval("${params.bamadap_bin} --version | cut -f 2 -d ' '"), topic: versions

    script:
    chunk = "${library}_${start_offset}_${end_offset}"
    def trim_cpus  = 4
    def align_cpus = Math.max(1, (task.cpus * 7).intdiv(8))
    def zip_cpus   = 2
    def sort_cpus  = Math.max(1, task.cpus.intdiv(8))
    def bwameth_pairing = params.single_end ? '' : '-p'
    def slice_path = "${params.tmp_dir}/${chunk}.slice.bam"
    """
    set +o pipefail
    inst_name=\$(samtools view ${bam} | head -n 1 | cut -d ":" -f 1)
    rg_line=\$(samtools view -H ${bam} | grep "^@RG" | sed 's/\\t/\\\\t/g' | head -n1)
    set -o pipefail
    trim_polyg=\$(echo "\${inst_name}" | awk '{if (\$1~/^A0|^NB|^NS|^VH|^LH/) {print "--trim-poly-g"} else {print ""}}')

    # we write this to a temporary file on local storage so 
    # that both bamadap and fgumi can access it simply
    bamslice --input ${bam} --start-offset ${start_offset} --end-offset ${end_offset} \\
        --format bam --output ${slice_path}

    ${params.bamadap_bin} --input ${slice_path} --format fastq --no-mate-suffix \\
        -l 2 --disable-quality-filtering \\
        --adapter-fasta ${adapter_fasta} \${trim_polyg} \\
        --threads ${trim_cpus} --json "${chunk}.fastp.json" \\
    | bwameth.py -t ${align_cpus} --read-group "\${rg_line}" --reference ${bwa_index} ${bwameth_pairing} - 2> "${library}.log.bwamem" \\
    | mark-nonconverted-reads.py --reference ${bwa_index} 2> "${chunk}.nonconverted_counts.tsv" \\
    | fgumi zipper --unmapped ${slice_path} --reference ${genome_fa} \\
        --exclude-missing-reads true --threads ${zip_cpus} --compression-level 0 \\
    | samtools sort -T ${params.tmp_dir}/samtools_sort_tmp -@ ${sort_cpus} \\
        -m ${(task.memory.toGiga()*5).intdiv(8)}G -o "${chunk}.aln.bam" /dev/stdin

    rm -f ${slice_path}
    """
}
