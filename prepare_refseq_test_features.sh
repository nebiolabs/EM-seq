#!/bin/bash
# Extract one example of each required feature type from RefSeq GFF file
# for testing the feature_cov workflow

set -euo pipefail

# Feature types we need
FEATURE_TYPES=(
    "promoter"
    "transcriptional_cis_regulatory_region"
    "enhancer"
    "mobile_genetic_element"
    "primary_transcript"
    "lnc_RNA"
    "exon"
)

# Input/output files
REFSEQ_URL="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/GCF_009914755.1-RS_2023_10/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz"
OUTPUT_FILE="tests/fixtures/features/test_refseq.gff"
TMP_FILE=$(mktemp)

echo "Downloading and filtering RefSeq annotations..."
echo "Looking for feature types: ${FEATURE_TYPES[*]}"

# Download and extract header
curl -fsSL "$REFSEQ_URL" | gunzip | grep '^#' > "$OUTPUT_FILE"

# Extract first occurrence of each feature type
for feature in "${FEATURE_TYPES[@]}"; do
    echo "  Extracting first '$feature'..."
    curl -fsSL "$REFSEQ_URL" \
        | gunzip \
        | grep -v '^#' \
        | awk -F'\t' -v feat="$feature" '$3 == feat' \
        | head -1 \
        >> "$TMP_FILE"
done

# Sort by chromosome and position, append to output
sort -k1,1 -k4,4n "$TMP_FILE" >> "$OUTPUT_FILE"

# Compress
gzip -f "$OUTPUT_FILE"

rm "$TMP_FILE"

echo "Created tests/fixtures/features/test_refseq.gff.gz"
echo "Feature counts:"
gunzip -c tests/fixtures/features/test_refseq.gff.gz | grep -v '^#' | awk -F'\t' '{print $3}' | sort | uniq -c
