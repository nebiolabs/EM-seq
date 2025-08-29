process alignReads {
    label 'high_cpu'
    tag { library }
    conda "conda-forge::python=3.10 bioconda::bwameth=0.2.7 bioconda::mark-nonconverted-reads=1.2 bioconda::samtools=1.22"
    publishDir "${params.outputDir}/bwameth_align", mode: 'symlink'

    input:
        tuple val(library), path(bam), val(chunk_name), path(read1), path(read2)
        val(genome_fa)

    output:
        tuple val(library), path("*.aln.bam"), path("*.aln.bam.bai"), emit: bam_files
        tuple val(library), path("*.nonconverted_counts.tsv"), emit: nonconverted_counts
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions
        tuple val("${task.process}"), val('bwameth'), eval('bwameth.py --version | sed \'s/^bwa-meth.py //\''), topic: versions
        tuple val("${task.process}"), val('python'), eval('python --version | sed \'s/^Python //\''), topic: versions

    script:
    """
    echo "Input file size: ${read1.size() / (1024 * 1024 * 1024)} GB"
    echo "Memory allocated: ${task.memory}" 

    get_rg_line() {
        set +o pipefail
        local file=\$1
        rg_line=\$(samtools view -H \$file | grep "^@RG" | sed 's/\\t/\\\\t/g' | head -n1)
        set -o pipefail
        export rg_line
    }
    get_rg_line ${bam}

    bwameth.py -p -t ${Math.max(1,(task.cpus*7).intdiv(8))} --read-group "\${rg_line}" --reference ${genome_fa} ${read1} ${read2} 2> "${library}.log.bwamem" \
    | mark-nonconverted-reads.py --reference ${genome_fa} 2> "${chunk_name}.nonconverted_counts.tsv" \
    | samtools view -u /dev/stdin \
    | samtools sort -T ${params.tmp_dir}/samtools_sort_tmp -@ ${Math.max(1,task.cpus.intdiv(8))} \
       -m ${(task.memory.toGiga()*5).intdiv(8)}G --write-index \
       -o "${chunk_name}.aln.bam##idx##${chunk_name}.aln.bam.bai" /dev/stdin


    """
}
