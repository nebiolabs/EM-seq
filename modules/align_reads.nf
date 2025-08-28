process alignReads {
    label 'high_cpu'
    tag { library }
    conda "conda-forge::python=3.10 bioconda::bwameth=0.2.7 bioconda::fastp=0.26 bioconda::mark-nonconverted-reads=1.2 bioconda::samtools=1.22"
    publishDir "${params.outputDir}/bwameth_align", mode: 'symlink'
    memory {
        try { 
            def fileSize = read1.size() / (1024 * 1024 * 1024)
            if (fileSize < 1.8) return '64 GB'
            else if (fileSize < 6.5) return '128 GB'
            else return '256 GB'
        }
        catch (Exception _e) {
            return '128 GB'  // Default memory if size cannot be determined
        }
    }

    input:
        tuple val(library), val(chunk_name), path(read1), path(read2)
        val(genome_fa)

    output:
        tuple val(library), path("*.aln.bam"), path("*.aln.bam.bai"), emit: bam_files
        tuple val(library), path("*.nonconverted.tsv"), emit: nonconverted_counts

    script:
    """
    echo "Input file size: ${read1.size() / (1024 * 1024 * 1024)} GB"
    echo "Memory allocated: ${task.memory}" 

    barcodes_from_fastq() {
        set +o pipefail
        zcat -f \$1 \
        | head -n10000 \
        | awk '{
            if (NR%4==1) {
                split(\$0, parts, ":"); 
                arr[ parts[ length(parts) ] ]++
            }} END { for (i in arr) {print arr[i]"\\t"i} }' \
        | sort -k1nr | head -n1 | cut -f2 
        set -o pipefail
    }

    get_rg_line() {
        set +o pipefail
        local file=\$1
        barcodes=(\$(barcodes_from_fastq \$file))
        rg_line="@RG\\tID:\${barcodes}\\tSM:${library}\\tBC:\${barcodes}"
        set -o pipefail
        export rg_line
    }

    get_rg_line ${read1}

    bwameth.py -p -t ${Math.max(1,(task.cpus*7).intdiv(8))} --read-group "\${rg_line}" --reference ${genome_fa} ${read1} ${read2} 2> "${library}.log.bwamem" \
    | mark-nonconverted-reads.py --reference ${genome_fa} 2> "${library}.nonconverted.tsv" \
    | samtools view -u /dev/stdin \
    | samtools sort -T ${params.tmp_dir}/samtools_sort_tmp -@ ${Math.max(1,task.cpus.intdiv(8))} \
       -m ${(task.memory.toGiga()*5).intdiv(8)}G --write-index \
       -o "${chunk_name}.aln.bam##idx##${chunk_name}.aln.bam.bai" /dev/stdin


    """
}
