function print_contig_meth(chr, sumMeth, sumMethDups, total, totalDups, total_count, min_pos, max_pos) {
    if (sumMeth["CH"]/total_count >= 0.00000001 ) { # ultra low frequency methylation sites are not reported
        print chr, sumMeth["CG"]/total["CG"]*100, sumMethDups["CG"]/totalDups["CG"]*100, 
                    sumMeth["CH"]/total["CH"]*100, sumMethDups["CH"]/totalDups["CH"]*100
    }
}
BEGIN{
    OFS="\t"; total_count=0;
    print "chr", "CG %meth", "CG %meth(+dups)", "CH %meth", "CH %meth(+dups)"
}
{

    total_count+=$6+$7;
    if(chr != $1 && NR != 1){
        print_contig_meth(chr, sumMeth, sumMethDups, total, totalDups, total_count, min_pos, max_pos);
        sumMeth["CG"]=0; sumMethDups["CG"]=0; total["CG"]=0; totalDups["CG"]=0;
        sumMeth["CH"]=0; sumMethDups["CH"]=0; total["CH"]=0; totalDups["CH"]=0;
    }
    chr=$1; context=($2=="CpG" ? "CG" : "CH");
    if ($5 > min_pos && $5 < max_pos) {
        sumMeth[context]+=$6; sumMethDups[context]+=$8; total[context]+=$6+$7; totalDups[context]+=$8+$9
    }
}
END{
    print_contig_meth(chr, sumMeth, sumMethDups, total, totalDups, total_count, min_pos, max_pos);
}