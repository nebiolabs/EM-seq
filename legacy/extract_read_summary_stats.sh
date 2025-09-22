#!/bin/bash

# this program extracts field sfrom  samtools flagstat to a tabular file suitable for tableau reporting

#9791    173     total (QC-passed reads + QC-failed reads)
#9787    151     primary
#4       22      secondary
#0       0       supplementary
#5496    62      duplicates
#5496    62      primary duplicates
#9790    173     mapped
#99.99%  100.00% mapped %
#9786    151     primary mapped
#99.99%  100.00% primary mapped %
#9787    151     paired in sequencing
#4901    79      read1
#4886    72      read2
#9719    0       properly paired
#99.31%  0.00%   properly paired %
#9785    151     with itself and mate mapped
#1       0       singletons
#0.01%   0.00%   singletons %
#50      29      with mate mapped to a different chr
#2       0       with mate mapped to a different chr (mapQ>=5)


echo -e 'file\tprimary_reads\tprimary_duplicates\tprimary_mapped\tusable_reads'
for f in "$@"; do 

  echo -ne "$f\t"
  samtools flagstat -O tsv $f | grep -E -e 'primary$' -e 'primary duplicates$' -e 'primary mapped$' \
    | awk -v FS='\t' -v OFS='\t' '{print $1+$2}' | paste - - - \
    | tr '\n' '\t'
  samtools view -F 0xF00 -q 10 -c $f 
done 

