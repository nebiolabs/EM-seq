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
read_len=$(cut -f 5 "$combined_meth_tsv" | head -n "$max_possible_read_length" | sort -n | tail -n 1)
max_pos=$((read_len - $bases_to_skip))

awk -v min_pos="$min_pos" -v max_pos="$max_pos" '
BEGIN{
    OFS="\t"; total_count=0;
    print "chr", "CG %meth", "CG %meth(+dups)", "CH %meth", "CH %meth(+dups)"
}
{
    total_count+=$6+$7;
    if(chr != $1 && NR != 1){
        if (sumMeth["CH"]/total_count >= 0.00000001 ) { # ultra low frequency methylation sites are not reported
            print chr, sumMeth["CG"]/total["CG"]*100, sumMethDups["CG"]/totalDups["CG"]*100, 
                       sumMeth["CH"]/total["CH"]*100, sumMethDups["CH"]/totalDups["CH"]*100
        }
        sumMeth["CG"]=0; sumMethDups["CG"]=0; total["CG"]=0; totalDups["CG"]=0;
        sumMeth["CH"]=0; sumMethDups["CH"]=0; total["CH"]=0; totalDups["CH"]=0;
    }
    chr=$1; context=($2=="CpG" ? "CG" : "CH");
    if ($5 > min_pos && $5 < max_pos) {
        sumMeth[context]+=$6; sumMethDups[context]+=$8; total[context]+=$6+$7; totalDups[context]+=$8+$9
    }
}
END{
    if (sumMeth["CH"]/total_count >= 0.00000001) { # ultra low frequency methylation sites are not reported
        print chr, sumMeth["CG"]/total["CG"]*100, sumMethDups["CG"]/totalDups["CG"]*100, 
                   sumMeth["CH"]/total["CH"]*100, sumMethDups["CH"]/totalDups["CH"]*100
    }
}
' <(tail -n +2 "$combined_meth_tsv")