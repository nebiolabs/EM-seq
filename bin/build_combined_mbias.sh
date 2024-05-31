#!/usr/bin/env bash -e

# This script generates a combined m-bias file for a given BAM file
# The m-bias file contains information about the methylation status of each read position in the BAM file
# broken down by chromosome and context (CpG, CHG, CHH)
# The script requires the presence of the executables 'samtools' and 'MethylDackel' in the system PATH.

# Usage: build_combined_mbias.sh <bam_file>

# Parameters:
# - bam_file: The input BAM file for which the m-bias file needs to be generated.

# Note: The script assumes that there is a reference genome file (.fa or .fasta) specified in the BAM file headers.

# Example usage:
# $ ./build_combined_mbias.sh input.bam

# Output:
# The script generates the combined m-bias data and prints it to the standard output.

# Dependencies:
# - samtools: Required for extracting information from BAM file headers.
# - MethylDackel: Required for calculating m-bias for each read position.


bam=$1

if [ -f "$bam" ]; then
    echo "Usage: $0 <bam_file>"
    exit 1
fi
SAMTOOLS=$(command -v samtools)
METHYLDACKEL=$(command -v MethylDackel)

if [[ -z $SAMTOOLS || -z $METHYLDACKEL ]]; then
    echo "Error: Required executables (samtools or MethylDackel) not found."
    exit 1
fi
#depends on there being a .fa or .fasta file in the headers of the bam file
ref=(`$SAMTOOLS view -H $bam | sed -n -e 's!^.* \([^ ]*\.fa[^ ]*\).*!\1!p'`)
cpus=2

chrs=(`$SAMTOOLS view -H $bam | grep @SQ | cut -f 2 | sed 's/SN://'| grep -v _random | grep -v chrUn | sed 's/|/\|/'`)

echo -e "chr\tcontext\tstrand\tRead\tPosition\tnMethylated\tnUnmethylated\tnMethylated(+dups)\tnUnmethylated(+dups)"
for chr in ${chrs[*]}; do
    for context in CHH CHG CpG; do
        arg=''
        if [ $context = 'CHH' ]; then
           arg='--CHH --noCpG'
        elif [ $context = 'CHG' ]; then
           arg='--CHG --noCpG'
        fi
        join -t $'\t' -j1 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -a 1 -e 0 \
        <(MethylDackel mbias --noSVG                     $arg -@ $cpus -r $chr $ref $bam | tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1)\
        <(MethylDackel mbias --noSVG --keepDupes -F 2816 $arg -@ $cpus -r $chr $ref $bam | tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1)\
        | sed "s/^/${chr}\t${context}\t/"
    done
done