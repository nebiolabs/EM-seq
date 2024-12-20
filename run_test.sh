pwd=$(pwd)
tmp="${pwd}/test_data/tmp"
mkdir ${tmp}

# make 151 nt-long bam and fastq files WITH simulated conversions #
# --------------------------------------------------------------- #

ln -s ${pwd}/test_data/emseq-test_R*.fastq ${tmp}/

mkfifo tmp_fq
paste -d "\n" <(samtools sort ${tmp}/emseq-test_R1.fastq | samtools view | awk 'BEGIN{OFS="\t"}{$2=77; print $0"\tBC:Z:CGTCAAGA-GGGTTGTT\tRG:Z:NS500.4"}') <(samtools sort ${tmp}/emseq-test_R2.fastq | samtools view | awk 'BEGIN{OFS="\t"}{$2=141; print $0"\tBC:Z:CGTCAAGA-GGGTTGTT\tRG:Z:NS500.4"}') > tmp_fq &
tmp_fq_pid=$!

cat <(echo -e "@HD\tVN:1.6\tSO:unsorted\n@RG\tID:test1\tSM:testsample1") <(cat tmp_fq | \
awk '
    function modify_base(base, next_base) {
        if (base == "C") {
            if (next_base == "G") {
                return (rand() < 0.2) ? "T" : "C";
            } else {
                return "T";
            }
        } return base;} 
        BEGIN { srand() } {
        modified_seq ="" 
        for (i = 1; i <= length($10); i++) {
            base = substr($10, i, 1);
            next_base = (i < length($10)) ? substr($10, i + 1, 1) : "";
            modified_seq = modified_seq modify_base(base, next_base);
        }
        $10 = modified_seq
        print $0; }' | tr " " "\t") | \
     samtools view -h -o ${tmp}/emseq-test.u.bam

wait $tmp_fq_pid
rm tmp_fq


# make genome and index #
# ---------------------

gen_rand_seq() { 
    length=$1
    tr -dc 'ACGT' < /dev/urandom | head -c ${length}
}
export -f gen_rand_seq
revcomp() {
    echo $1 | rev | tr "[ATCGNatcgn]" "[TAGCNtagcn]"
}
export -f revcomp


echo ">chr_Human_autosome_chr1" > ${tmp}/reference.fa

samtools view ${tmp}/emseq-test.u.bam | \
    cut -f10 | paste - - | \
    awk 'BEGIN{srand()}{
        min=20; 
        max=700; 
        random_number = int(min + rand() * (max - min + 1)); 
        r=system("gen_rand_seq " random_number);
        read2=system("revcomp " $2);
        # read2=$2
        print $1""r""read2
    }' | tr -d "\n" | fold -w60 >> ${tmp}/reference.fa
    #| shuf | tr -d "\n" | fold -w60 >> ${tmp}/reference.fa



# make 76 nt-long fastq and bam files INCLUDING simulated conversions #
# ------------------------------------------------------------------- #

cat <(samtools view -H ${tmp}/emseq-test.u.bam) <(samtools view ${tmp}/emseq-test.u.bam | awk 'BEGIN{OFS="\t"}{$10=substr($10,1,76); $11=substr($11,1,76); print $0}') > ${tmp}/emseq-test_76.u.bam

BARCODE=$(samtools view -f 77 ${tmp}/emseq-test.u.bam | grep -o "BC:Z:[A-Z,+,-]*" | sort | uniq -c | head | awk '{gsub("BC:Z","N:0",$2); print $2}' )
samtools view -f 77 ${tmp}/emseq-test_76.u.bam | samtools fastq | sed "s/\/1$/ 1:$BARCODE/"  > ${tmp}/emseq-test_76.R1.fastq
samtools view -f 141 ${tmp}/emseq-test_76.u.bam | samtools fastq | sed "s/\/2$/ 2:$BARCODE/" > ${tmp}/emseq-test_76.R2.fastq


# make gzip fastq files. #
# ---------------------- #

gzip -c ${tmp}/emseq-test_76.R1.fastq > ${tmp}/emseq-test_76.R1.fastq.gz
gzip -c ${tmp}/emseq-test_76.R2.fastq > ${tmp}/emseq-test_76.R2.fastq.gz



# run tests #
# --------- #

genome_path=${tmp}/reference.fa

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
