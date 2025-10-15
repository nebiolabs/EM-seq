#this window of refseq contains all the features we're interested in
REFSEQ_URL="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/GCF_009914755.1-RS_2023_10/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz"
curl -fsSL -O "$REFSEQ_URL"
cat < GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz | head -n 85000 | tail -n 40000 | bgzip > test_refseq.gff.gz
