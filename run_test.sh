#!/usr/bin/env bash
set -euo pipefail


# user might need custom config file, don't overwrite it
if [ ! -f nextflow.config ]; then
    echo "Copying example nextflow.config to nextflow.config"
    cp nextflow.config.example nextflow.config
fi

# set up tmp folder and copy test data into it
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
tmp="${script_dir}/test_data/tmp"
test_log="${script_dir}/test.log.out"

if [ ! -d "${tmp}" ]; then
    mkdir -p "${tmp}"
fi

for fq in "${script_dir}"/test_data/emseq-test*.fastq.gz; do
    basename=$(basename "$fq")
    ln -sf "$fq" "${tmp}/${basename}"
    ln -sf "$fq" "${tmp}/2${basename}"
done


# Check for micromamba
if ! command -v micromamba >/dev/null 2>&1; then
    echo "Error: micromamba not found in PATH. Please install micromamba or adjust your PATH."
    exit 1
fi
eval "$(micromamba shell hook --shell=bash)"
if [ "${GITHUB_ACTIONS:-}" == "true" ]; then
    echo "gh actions..."
else
    #see if nextflow.emseq env exists, if not create it
    micromamba activate nextflow.emseq || \
        micromamba create --name nextflow.emseq --yes python=3 bioconda:nextflow bioconda::samtools
fi
# if running on osX need to set CONDA_OVERRIDE_OSX=<current osx version>
if [[ "$OSTYPE" == "darwin"* ]]; then
    export CONDA_OVERRIDE_OSX=$(sw_vers -productVersion | cut -d '.' -f 1)
    echo "Setting CONDA_OVERRIDE_OSX to ${CONDA_OVERRIDE_OSX}"
fi

micromamba activate nextflow.emseq

# generate test data from minimal set of reads (fq/fq.gz/bam)
echo "Generating test reads..."

# Extract fastq files
echo "Extracting fastq files..."
if ! gunzip -c "${tmp}/emseq-testg_R1.fastq.gz" > "${tmp}/emseq-test_R1.fastq"; then
    echo "Error: Failed to extract R1 fastq"
    exit 1
fi

if ! gunzip -c "${tmp}/emseq-testg_R2.fastq.gz" > "${tmp}/emseq-test_R2.fastq"; then
    echo "Error: Failed to extract R2 fastq"
    exit 1
fi

# Generate unaligned BAM file
echo "Generating unaligned BAM file..."
if ! samtools import \
    -1 "${tmp}/emseq-test_R1.fastq" \
    -2 "${tmp}/emseq-test_R2.fastq" \
    -o "${tmp}/emseq-test.u.bam" \
    --barcode-tag BC \
    --barcode CGTCAAGA-GGGTTGTT \
    --rg-line $'@RG\tID:NS500.4\tSM:emseq-test\tPL:ILLUMINA'; then
    echo "Error: Failed to generate unaligned BAM file"
    exit 1
fi

# Verify the output BAM file is valid
if ! samtools quickcheck -u "${tmp}/emseq-test.u.bam"; then
    echo "Error: Generated BAM file failed validation"
    exit 1
fi


# run tests #
# --------- #

echo "Running nextflow pipeline..."

genome_path=$(realpath "test_data/reference.fa")
target_bed=$(realpath "test_data/emseq_test_regions.bed")

pushd "${tmp}" || {
    echo "Error: Failed to change to tmp directory"
    exit 1
}

# Initialize test log
rm -f "${test_log}"

function test_pipeline() {
    local file="$1"
	local bed_argument="$2"
    echo -n "Testing with file: ${file} ... " | tee -a "${test_log}"

    # Run the Nextflow pipeline with the specified input file
    if nextflow run "${script_dir}/main.nf" \
        --input_glob "${file}" \
        --path_to_genome_fasta "${genome_path}" \
        --email "me@example.com" \
        --max_input_reads 10000 \
        --flowcell "test" \
        -with-report "emseq_metadata_report.html" \
        -with-timeline "emseq_metadata_timeline.html" \
        -with-dag "emseq_metadata_dag.html" \
        -w "${tmp}/work" \
        --read_length 151 \
        --testing_mode "true" \
		${bed_argument} \
        --enable_neb_agg "false" 2>&1 >> "${test_log}"; then

        echo "Nextflow pipeline succeeded" >> "${test_log}"
    else
        echo "Nextflow pipeline failed" >> "${test_log}"
        return 1
    fi

    # Check results
    echo "Test complete. Checking results..." | tee -a "${test_log}"

    return 0
}

# Run tests with different input formats
echo "Running pipeline tests..."

bed_arg1="--target_bed ${target_bed}"
bed_arg2=""

if ! test_pipeline "*emseq-test*1.fastq.gz" "${bed_arg1}"; then
    echo "❌ Test failed for fastq.gz files"
    exit 1
fi

if ! test_pipeline "emseq-test*1.fastq" "${bed_arg1}"; then
    echo "❌ Test failed for fastq files"
    exit 1
fi

if ! test_pipeline "emseq-test*bam" "${bed_arg1}"; then
    echo "❌ Test failed for bam files"
    exit 1
fi

if ! test_pipeline "emseq-test*bam" "${bed_arg2}"; then
    echo "❌ Test failed for bam files"
    exit 1
fi


echo "echo ✅ All tests passed!"

# Optional cleanup (uncomment if desired)
# echo "Cleaning up temporary files..."
# rm -rf "${tmp}"

# Display test results
echo "Test Results:"
cat "${test_log}"
