#!/usr/bin/env bash -e

# This script takes a combined methylation TSV file as input and calculates the percentage of methylation for CpG and non-CpG sites for each chromosome.
# The input file should have the following columns: chromosome, context (CpG or CH), and methylation information.
# The output is a tab-separated file with the following columns: chromosome, percentage of methylation for CpG sites, percentage of methylation for CpG sites including duplicates, percentage of methylation for non-CpG sites, percentage of methylation for non-CpG sites including duplicates.

# Usage: summarize_combined_mbias.sh <combined_meth_tsv>

# Parameters:
#   - combined_meth_tsv: Path to the combined methylation TSV file.

# Example usage:
#   summarize_combined_mbias.sh input.tsv

combined_meth_tsv="$1"
bases_to_skip=5

if [ ! -f "$combined_meth_tsv" ]; then
    echo "Usage: $0 <combined_meth_tsv>"
    exit 1
fi
max_possible_read_length=1000
min_pos=$bases_to_skip
read_len=$(cut -f 5 "$combined_meth_tsv" | head -n $max_possible_read_length | sort -n | tail -n 1)
max_pos=$((read_len - $bases_to_skip))
script_path=$(dirname "$0")

awk -v min_pos="$min_pos" -v max_pos="$max_pos" \
 -f "${script_path}/per_contig_methylation_from_combined_mbias.awk" <(tail -n +2 "$combined_meth_tsv")