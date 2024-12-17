pwd=$(pwd)
tmp="${pwd}/test_data/tmp"
mkdir ${tmp}

# make 151 nt-long bam and fastq files #
# -------------
ln -s ${pwd}/test_data/emseq-test_R*.fastq ${tmp}/

cat <(echo -e "@HD\tVN:1.6\tSO:unsorted\n@RG\tID:test1\tSM:testsample1") <(samtools sort ${tmp}/emseq-test_R1.fastq | samtools view | awk 'BEGIN{OFS="\t"}{$2=77; print $0"\tBC:Z:CGTCAAGA-GGGTTGTT\tRG:Z:NS500.4"}') <(samtools sort ${tmp}/emseq-test_R1.fastq | samtools view | awk 'BEGIN{OFS="\t"}{$2=141; print $0"\tBC:Z:CGTCAAGA-GGGTTGTT\tRG:Z:NS500.4"}') | samtools view -ho ${tmp}/emseq-test.u.bam

# make genome and index #
# ---------------------

echo ">chr_Human_autosome_chr1" > ${tmp}/reference.fasta
samtools view ${tmp}/emseq-test.u.bam | shuf | awk 'BEGIN{ srand() } {
        if ($10 ~ /GC/) {
                while ($10 ~ /GC/) {
                        if (rand() > 0.5) {
                                gsub(/GC/, "GT")
                        }
                }
        }
        seq = seq""$10
}END{ print seq }' | fold -w60 >> ${tmp}/reference.fasta

######################## bwa 


# make 76 nt-long fastq and bam files #
# -------------
for rn in 1 2; do cat ${tmp}/emseq-test_R${rn}.fastq | awk '{if (NR%4==2 || NR%4==0) {print substr($0,1,76)} else {print $0} }' > ${tmp}/emseq-test_76.R${rn}.fastq; done

cat <(echo -e "@HD\tVN:1.6\tSO:unsorted\n@RG\tID:test1\tSM:testsample1") <(samtools sort ${tmp}/emseq-test_76.R1.fastq | samtools view | awk 'BEGIN{OFS="\t"}{$2=77; print $0"\tBC:Z:CGTCAAGA-GGGTTGTT\tRG:Z:NS500.4"}') <(samtools sort ${tmp}/emseq-test_76.R1.fastq | samtools view | awk 'BEGIN{OFS="\t"}{$2=141; print $0"\tBC:Z:CGTCAAGA-GGGTTGTT\tRG:Z:NS500.4"}') | samtools view -ho ${tmp}/emseq-test.76.u.bam

# make gzip fastq files. #
# -------------
gzip -c ${tmp}/emseq-test_76.R1.fastq > ${tmp}/emseq-test_76.R1.fastq.gz
gzip -c ${tmp}/emseq-test_76.R2.fastq > ${tmp}/emseq-test_76.R2.fastq.gz



genome_path=${tmp}/reference.fasta

#. /mnt/home/aerijman/miniconda3/bin/activate
. ~/miniconda3/bin/activate
conda activate nextflow

pushd ${tmp}
# loop for different type of files, fastq, fastq.gz and bam

nextflow run ${pwd}/main.nf  \
  --input_glob "*1.fastq.gz" \
  --genome "${genome_path}" \
  --email "aerijman@neb.com" \
  --max_input_reads 10000 \
  --flowcell "test_pipeline" \
  -with-report  "emseq_metadata_report.html" \
  -with-timeline "emseq_metadata_timeline.html" \
  -with-dag "emseq_metadata_dag.html" \
  -w "${tmp}/work" \
  --read_length 151 \
  --enable_neb_agg "false" \
  -resume

# /mnt/hpc_scratch/" \
