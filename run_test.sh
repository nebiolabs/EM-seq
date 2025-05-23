#!/usr/bin/env bash
set -euo pipefail


pwd=$(pwd)
tmp="${pwd}/test_data/tmp"
[ -d "${tmp}" ] || mkdir -p "${tmp}"


# check conda is installed OR install it #
# -------------------------------------- #
if command -v conda &> /dev/null; then
    activate_conda=$(type -a conda | grep -Eo '/.*conda' | head -n 1 | sed 's/\/bin\|\/condabin/\/bin/' | sed 's/bin.*/bin\/activate/')
else
    echo "please install conda. Then run this script"
    exit 1
fi
. ${activate_conda} || echo "Perhaps conda is not in: $(type -a conda) ?"


if ! conda env list | grep -q "nextflow.emseq"; then 
    conda update conda -y
    conda update --all -y
    conda config --add channels defaults 
    conda config --add channels conda-forge 
    conda config --add channels bioconda 
    conda install -c conda-forge libgcc-ng=12 -y

    conda create --name nextflow.emseq --yes python=3 && conda install --name nextflow.emseq --yes bioconda:nextflow bioconda::samtools
else
    echo "conda env found."
fi

if [ "${GITHUB_ACTIONS:-}" == "true" ]; then
    exec $SHELL
fi

conda init
conda activate nextflow.emseq

# make 151 nt-long bam and fastq files WITH simulated conversions #
# --------------------------------------------------------------- #

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

genome_path=${tmp}/reference.fa

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
cat em-seq_output/stats/flagstats/emseq-test_76..flagstat |\
    grep -q "5000 + 0 properly paired" && echo "flagstats OK" >> ${pwd}/test.log.out || echo "flagstats not OK" >> ${pwd}/test.log.out
tail -n2 em-seq_output/stats/picard_alignment_metrics/emseq-test_76..alignment_summary_metrics.txt |\
    awk 'BEGIN{result="alignment metrics not OK"}{if ($1==76 && $2>4996) {result="alignment metrics OK"}}END{print result}' >> ${pwd}/test.log.out

rm -r ${tmp}
echo "FINISHED"

