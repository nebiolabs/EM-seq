genome=$1
flowcell=$2
insert_read1=$3
insert_read2=$4
library=$5
cpus=$6
tile=$7
lane=$8

rg_id="@RG\\tID:${fastq_barcode}\\tSM:${library}" 
bwa_mem_log_filename="${library}_${fastq_barcode}${flowcell}_${lane}_${tile}.log.bwamem"
bam_filename="${library}_${fastq_barcode}_${flowcell}_${lane}_${tile}.aln.bam"

inst_name=$(zcat -f ${insert_read1} | head -n 1 | cut -f 1 -d ':' | sed 's/^@//')
fastq_barcode=$(zcat -f ${insert_read1} | head -n 1 | sed -r 's/.*://')

trim_polyg=$(echo ${inst_name} | grep -q "^A0\|^NB\|^NS\|^VH" && echo "--trim_poly_g" || echo "")
echo ${trim_polyg} | \
awk '{
    if (length($1)>0) { print "2-color instrument: poly-g trim mode on" } 
}'

echo "$rg_id $inst_name $bam_filename $bwa_mem_log_filanme"

#seqtk mergepe <(zcat -f ${insert_read1}) <(zcat -f ${insert_read2}) \
#| fastp --stdin --stdout -l 2 -Q ${trim_polyg} --interleaved_in --overrepresentation_analysis -j ${library}_fastp.json 2> fastp.stderr \
#| bwameth.py -p -t ${task.cpus} --read-group ${rg_id} --reference ${genome} /dev/stdin 2> ${bwa_mem_log_filename} \
#| mark-nonconverted-reads.py 2> "${library}_${fastq_barcode}_${flowcell}_${lane}_${tile}.nonconverted.tsv" \
#| sambamba view -t 2 -S -f bam -o ${bam_filename} /dev/stdin 2> sambamba.stderr;
