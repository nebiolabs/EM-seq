process tasmanian {
    tag "${library}"
    label 'process_low'
        
	
    conda 'bioconda::samtools=1.21 bioconda::tasmanian-mismatch=2.0.5'

    input:
        tuple val(library), path(bam), path(bai)
        each path(reference_fasta)
        each path (masked_bedgraph)

    output:
        tuple val(library), path("${library}.tasmanian.tsv"), emit: tasmanian_for_aggregate
        tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed \'s/^samtools //\''), topic: versions
        tuple val("${task.process}"), val('tasmanian-mismatch'), val('*should be* 2.0.5'), topic: versions

    script:
    """
    set +e
    set +o pipefail

    ln -s \$(readlink -f ${reference_fasta}) ref.fa
    samtools faidx ref.fa

    tasmanian-mismatch ${bam} ref.fa \
        --position-mode read \
        --min-base-quality 20 \
        --min-map-quality 20 \
        -b ${masked_bedgraph} \
        --bed-filter-mode mask \
        -o ${library}.tasmanian.read.tsv


    tasmanian-mismatch ${bam} ref.fa \
        --position-mode insert \
        --min-base-quality 20 \
        --min-map-quality 20 \
        -b ${masked_bedgraph} \
        --bed-filter-mode mask \
        -o ${library}.tasmanian.insert.tsv

	head -n1 ${library}.tasmanian.read.tsv | awk '{print $0"\tpair_mode"}' > ${library}.tasmanian.tsv 
	tail -n +2 ${library}.tasmanian.read.tsv | awk '{print $0"\tread"}' >> ${library}.tasmanian.tsv
	tail -n +2 ${library}.tasmanian.insert.tsv | awk '{print $0"\tinsert"}' >> ${library}.tasmanian.tsv
    """
}
