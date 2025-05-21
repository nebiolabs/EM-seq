#!/usr/bin/env bash
set -euo pipefail
set -x

# user migth need custom config file
if [ ! -f nextflow.config ]; then
    echo "Copying example nextflow.config to nextflow.config"
    cp nextflow.config.example nextflow.config
fi

# set up tmp folder and copy test data into it
pwd=$(pwd)
tmp="${pwd}/test_data/tmp"
[ -d "${tmp}" ] || mkdir -p "${tmp}"

ln -sf ${pwd}/test_data/emseq-test*.fastq.gz ${tmp}
ln -sf ${pwd}/test_data/reference.fa ${tmp}


if [ "${GITHUB_ACTIONS:-}" == "true" ]; then
    echo "gh actions..."
else
    micromamba create --name nextflow.emseq --yes python=3
fi
    micromamba install --name nextflow.emseq --yes bioconda:nextflow bioconda::samtools
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate nextflow.emseq

# generate test data from minimal set of reads (fq/fq.gz/bam)
 echo "generating reads"
 
 # fastq (perhaps make 76bp long in future test)
 gunzip -c ${tmp}/emseq-testg_R1.fastq.gz > ${tmp}/emseq-test_R1.fastq
 gunzip -c ${tmp}/emseq-testg_R2.fastq.gz > ${tmp}/emseq-test_R2.fastq

 # bam
 paste -d "\n" <(samtools sort ${tmp}/emseq-test_R1.fastq | samtools view | awk 'BEGIN{OFS="\t"}{$2=77; print $0"\tBC:Z:CGTCAAGA-GGGTTGTT\tRG:Z:NS500.4"}') \
	 		   <(samtools sort ${tmp}/emseq-test_R2.fastq | samtools view | awk 'BEGIN{OFS="\t"}{$2=141; print $0"\tBC:Z:CGTCAAGA-GGGTTGTT\tRG:Z:NS500.4"}') \
| samtools view -u -o ${tmp}/emseq-test.u.bam




# run tests #
# --------- #

echo "running nextflow pipeline..."


genome_path=${tmp}/reference.fa

pushd ${tmp}
rm -f ${pwd}/test.log.out
# loop for different type of files, fastq, fastq.gz and bam

echo "check if files are in here"
ls -ltr 
echo "finished listing files"

function test_pipeline {
    local file=$1

    # Run the Nextflow pipeline with the specified input file
    nextflow run ${pwd}/main.nf \
        --input_glob "${file}" \
        --path_to_genome_fasta ${genome_path} \
        --email "aerijman@neb.com" \
        --max_input_reads 10000 \
        --flowcell "test_pipeline" \
        -with-report  "emseq_metadata_report.html" \
        -with-timeline "emseq_metadata_timeline.html" \
        -with-dag "emseq_metadata_dag.html" \
        -w "${tmp}/work" \
        --read_length 151 \
        --enable_neb_agg "false" 2>&1 >> ${pwd}/test.log.out

        # Check if Nextflow run was successful
        if [ $? -ne 0 ]; then
            echo "Nextflow pipeline failed" >> ${pwd}/test.log.out 
            exit 1
        else
            echo "Nextflow pipeline succeeded" >> ${pwd}/test.log.out
        fi


        # Check results
        echo "checking results..."
        cat em-seq_output/stats/flagstats/emseq-testg.flagstat |\
            grep -q "1972 + 0 properly paired" && echo "flagstats OK" >> ${pwd}/test.log.out || echo "flagstats not OK" >> ${pwd}/test.log.out
        tail -n2 em-seq_output/stats/picard_alignment_metrics/emseq-testg.alignment_summary_metrics.txt |\
            awk 'BEGIN{result="alignment metrics not OK"}{if ($1==150 && $3>2238) {result="alignment metrics OK"}}END{print result}' >> ${pwd}/test.log.out
}

test_pipeline "emseq-test*1.fastq.gz"
test_pipeline "emseq-test*1.fastq"
test_pipeline "emseq-test*bam"


# rm -r ${tmp}

cat ${pwd}/test.log.out
echo "FINISHED"

