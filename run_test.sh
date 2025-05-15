#!/usr/bin/env bash
set -euo pipefail
set -x

cp nextflow.config.example nextflow.config


pwd=$(pwd)
tmp="${pwd}/test_data/tmp"
[ -d "${tmp}" ] || mkdir -p "${tmp}"

if [ "${GITHUB_ACTIONS:-}" == "true" ]; then
    #exec $SHELL
    echo "gh actions..."
else
    micromamba create --name nextflow.emseq --yes python=3
fi
    micromamba install --name nextflow.emseq --yes bioconda:nextflow bioconda::samtools
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate nextflow.emseq
#fi

# make 151 nt-long bam and fastq files WITH simulated conversions #
# --------------------------------------------------------------- #
echo "generating reads and genome reference data..."


ln -sf ${pwd}/test_data/emseq-test_R*.fastq ${tmp}/

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
echo ">chr_Human_autosome_chr1" > ${tmp}/reference.fa

gen_rand_seq() {
    length=$1
    LC_CTYPE=C tr -dc 'ACGT' < /dev/urandom | head -c ${length}
} 
revcomp() {
    echo $1 | rev | tr "[ATCGNatcgn]" "[TAGCNtagcn]"
} 
export -f gen_rand_seq revcomp

samtools view ${tmp}/emseq-test.u.bam | \
    cut -f10 | paste - - | \
    awk -v gen_rand_seq="gen_rand_seq" -v revcomp="revcomp" 'BEGIN{srand()}{
        min=20; 
        max=700;
        random_number = int(min + rand() * (max - min + 1)); 
        cmd1 = "gen_rand_seq "random_number 
        cmd1 | getline r;
        close(cmd1)
        cmd2 = "revcomp " $2 
        cmd2 | getline read2;
        close(cmd2)
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

echo "running nextflow pipeline..."

genome_path=${tmp}/reference.fa

pushd ${tmp}
# loop for different type of files, fastq, fastq.gz and bam

realpath "*1.fastq.gz"

nextflow run ${pwd}/main.nf \
  --input_glob "*1.fastq.gz" \
  --path_to_genome_fasta ${genome_path} \
  --email "aerijman@neb.com" \
  --max_input_reads 10000 \
  --flowcell "test_pipeline" \
  -with-report  "emseq_metadata_report.html" \
  -with-timeline "emseq_metadata_timeline.html" \
  -with-dag "emseq_metadata_dag.html" \
  -w "${tmp}/work" \
  --read_length 151 \
  --enable_neb_agg "false" \
  -resume 2>&1 > ${pwd}/test.log.out

# Check if Nextflow run was successful
if [ $? -ne 0 ]; then
    echo "Nextflow pipeline failed" >> ${pwd}/test.log.out 
    exit 1
else
    echo "Nextflow pipeline succeeded" >> ${pwd}/test.log.out
fi


# Check results
echo "checking results..."
cat em-seq_output/stats/flagstats/emseq-test_76.flagstat |\
    grep -q "5000 + 0 properly paired" && echo "flagstats OK" >> ${pwd}/test.log.out || echo "flagstats not OK" >> ${pwd}/test.log.out
tail -n2 em-seq_output/stats/picard_alignment_metrics/emseq-test_76.alignment_summary_metrics.txt |\
    awk 'BEGIN{result="alignment metrics not OK"}{if ($1==76 && $2>4996) {result="alignment metrics OK"}}END{print result}' >> ${pwd}/test.log.out

rm -r ${tmp} ${pwd}/nextflow.config 

cat ${pwd}/test.log.out
echo "FINISHED"

