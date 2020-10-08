// this consumes fastq files, adapter trims, then aligns reads to the specified reference using bwa-meth
// mode can be set to tile_fastqs or run_fastqs depending on whether the system should map each tile's reads in distinct jobs then combine (tile_fastqs)
// or all reads for a library in a single job (run_fastqs)

flowcell = params.flowcell
genome = params.genome
tmp_dir =  params.tmp_dir ?: ENV['tmp']
outputPath = params.outdir
fastq_mode = params.fastq_mode ?: 'run_fastqs'
println "Processing " + flowcell + "... => " + outputPath


Channel.from(params.fastqs_by_library).set{fq_set_channel}
    
process mapping {
    cpus fastq_mode == 'tile-fastq' ? 4 : 16
    errorStrategy 'retry'
    tag { [flowcell, fq_set.library] }
    conda "bwameth=0.2.2 seqtk=1.3 sambamba=0.7.0 fastp=0.20.0 mark-nonconverted-reads=1.1"

    input:
        val fq_set from fq_set_channel

    output:
        set val(fq_set.library), file("*.aln.bam") into aligned_files
        set val(fq_set.library), file("*.nonconverted.tsv") into nonconverted_counts

    shell:
    '''
    bwakit_path=\$(dirname \$(readlink -f \$(which run-bwamem)))
    #todo: check for presence of executables?
    export PATH="\${bwakit_path}:\${PATH}"

    seqtk mergepe <(zcat -f "!{fq_set.insert_read1}") <(zcat -f "!{fq_set.insert_read2}") \
    | trimadap -5 AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -3 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -3 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -3 ATCTCGTATGCCGTCTTCTGCTTG -3 CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -3 CTGTCTCTTATACACATCTGACGCTGCCGACGA 2> "!{fq_set.library}_!{fq_set.barcode}!{fq_set.flowcell}_!{fq_set.lane}_!{fq_set.tile}.log.trim" \
    | bwameth.py -p -t !{task.cpus} --read-group "@RG\tID:!{fq_set.barcode}\tSM:!{fq_set.library}" --reference !{genome} /dev/stdin 2>  "!{fq_set.library}_!{fq_set.barcode}!{fq_set.flowcell}_!{fq_set.lane}_!{fq_set.tile}.log.bwamem" \
    | mark-nonconverted-reads.py 2> "!{fq_set.library}_!{fq_set.barcode}!{fq_set.flowcell}_!{fq_set.lane}_!{fq_set.tile}.nonconverted.tsv" \
    | sambamba view -t 2 -S -f bam -o "!{fq_set.library}_!{fq_set.barcode}!{fq_set.flowcell}_!{fq_set.lane}_!{fq_set.tile}.aln.bam" /dev/stdin;
    '''

}

process mergeAndMarkDuplicates {
    cpus 8
    errorStrategy 'retry'
    tag { library }
    publishDir "${outputDir}", mode: 'copy', pattern: '*.{md.bam}*'
    conda "samtools=1.9 samblaster=0.1.24 sambamba=0.7.0"

    input:
        set val(library), file(libraryBam) from aligned_files.groupTuple()

    output:
        set val(library), file('*.md.bam'), file('*.md.bam.bai') into md_bams
        file('*.samblaster') into samblaster_logs

    shell:
    '''
    samtools cat  -b <( find . -name '*.aln.bam' ) \
    | samtools view -h /dev/stdin \
    | samblaster 2> !{library}.log.samblaster \
    | sambamba view -t 2 -l 0 -S -f bam /dev/stdin \
    | sambamba sort --tmpdir=!{tmp_dir} -t !{task.cpus} -m 20GB -o !{library}.md.bam /dev/stdin

    '''
}

    // we need to send the same file into multiple channels for downstream 
    // analysis whether they came from aligned bams or fastqs
    md_bams.into {  md_files_for_mbias; md_files_for_extract; md_files_for_fastqc; 
                    md_files_for_samstats; md_files_for_picard; 
                    md_files_for_goleft; md_files_for_picard_gc; md_files_for_samflagstats; 
                    md_files_for_aggregate; md_files_for_human_reads;
                 }
                 

    process methylDackel_mbias {
        cpus 8
        errorStrategy 'retry'
        tag {library}
        conda "methyldackel=0.4.0 samtools=1.9"

        input:
            tuple email, library, file(md_file), file(md_bai) from md_files_for_mbias.groupTuple(by: [0,1])

        output:
            tuple email, file('*.svg') into mbias_output_svg
            tuple email, file('*.tsv') into mbias_output_tsv
            tuple email, library, file('*.tsv') into mbias_for_aggregate

        shell:
        '''
        echo -e "chr\tcontext\tstrand\tRead\tPosition\tnMethylated\tnUnmethylated\tnMethylated(+dups)\tnUnmethylated(+dups)" > !{library}_combined_mbias.tsv
        chrs=(`samtools view -H !{md_file} | grep @SQ | cut -f 2 | sed 's/SN://'| grep -v _random | grep -v chrUn | sed 's/|/\\|/'`)

        for chr in ${chrs[*]}; do
            for context in CHH CHG CpG; do
                arg=''
                if [ $context = 'CHH' ]; then
                arg='--CHH --noCpG'
                elif [ $context = 'CHG' ]; then
                arg='--CHG --noCpG'
                fi
                # need two calls to add columns containing the counts without filtering duplicate reads (for rrEM-seq where start/end is constrained)
                # not sure why we need both --keepDupes and -F, probably a bug in mbias
                join -t $'\t' -j1 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -a 1 -e 0 \
                <( \
                    MethylDackel mbias --noSVG $arg -@ !{task.cpus} -r $chr !{genome} !{md_file} | \
                    tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
                ) \
                <( \
                    MethylDackel mbias --noSVG --keepDupes -F 2816 $arg -@ !{task.cpus} -r $chr !{genome} !{md_file} | \
                    tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
                ) \
                | sed "s/^/${chr}\t${context}\t/" \
                >> !{library}_combined_mbias.tsv
            done
        done
        # makes the svg files for trimming checks
        MethylDackel mbias -@ !{task.cpus} --noCpG --CHH --CHG -r ${chrs[0]} !{genome} !{md_file} !{library}_chn
        for f in *chn*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CHN <\\/text>/" $f; done;

        MethylDackel mbias -@ !{task.cpus} -r ${chrs[0]} !{genome} !{md_file} !{library}_cpg
        for f in *cpg*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CpG<\\/text>/" $f; done;

        '''

    }

    process methylDackel_extract {
        cpus 8
        tag {library}
        publishDir "${outputPath}", mode: 'copy'
        conda "methyldackel=0.4.0 pigz=2.4"

        input:
            tuple email, library, file(md_file), file(md_bai) from md_files_for_extract.groupTuple(by: [0,1])

        output:
            tuple email, library, file('*.methylKit.gz') into extract_output

        shell:
        '''
        MethylDackel extract --methylKit --OT 0,0,0,95 --OB 0,0,5,0 -@ !{task.cpus} --CHH --CHG -o !{library} !{genome} !{md_file}
        pigz -p !{task.cpus} *.methylKit
        '''

    }

    process select_human_reads {
        cpus 8
        tag {library}
        conda "sambamba=0.7.1i bedtools="

        input:
            tuple email, library, file(md_file), file(md_bai) from md_files_for_human_reads.groupTuple(by: [0,1])

        output:
            tuple email, library, file('*.human.bam') into human_bams_gc
            tuple email, library, file('*.human.bam') into human_bams_inserts

        shell:
        '''
        sambamba view -t 8 -l 0 -f bam !{md_file} chr1 chr2 chr3 chr4 chr5 chr6 \
                                                  chr7 chr8 chr9 chr10 chr11 chr12 \
                                                  chr13 chr14 chr15 chr16 chr17 chr18 \
                                                  chr19 chr20 chr21 chr22 chrX chrY \
        > !{md_file}.human.bam
        '''

    }

    process runFastQC {
        cpus 1
        errorStrategy 'retry'
        tag { library }
        conda "fastqc=0.11.8"

        input:
            tuple email, library, file(md_file), file(md_bai) from md_files_for_fastqc.groupTuple(by: [0,1])

        output:
            tuple email, file('*_fastqc.zip') into fastqc_results
            tuple email, library, file('*_fastqc.zip') into fastqc_results_for_aggregate

        shell:
        '''
        fastqc -f bam !{md_file}
        '''

    }
    process sum_nonconverted_reads {	

        input:	
            tuple email, library, file(count_files) from nonconverted_counts.groupTuple(by: [0,1])	

        output:	
            tuple email, file('*-nonconverted-counts.tsv') into cat_nonconversions	
            tuple email, library, file('*-nonconverted-counts.tsv') into nonconverted_counts_for_aggregate

        shell:	
        '''	
        files=(*.tsv)	
        paste *.tsv | awk -v numFiles=${#files[@]} -v OFS='\t' '	
        {	
        row = sep = ""	
        for(i=1; i < NF/numFiles; ++i) { row = row sep $i; sep = OFS }	
        sum = $(NF/numFiles) # last header col. / (1st) data col. to sum	
        for(i=2; i<=numFiles; ++i) sum += $(NF/numFiles * i) # add other cols.	
        printf "%s%s%s\\n", row, OFS, sum	
        }' > tmp-counts.tsv	
        awk '{print "!{library}\t" $0}' tmp-counts.tsv > !{library}-nonconverted-counts.tsv	
        '''    	
    }	

    process combine_nonconversion {	
        publishDir "${outputPath}", mode: 'copy'	

        input:	
            tuple email, file ('*') from cat_nonconversions.groupTuple()

        output:	
            file ("combined-nonconverted.tsv")	

        shell:	
        '''	
        cat *.tsv > combined-nonconverted.tsv	
        '''	

    }

    process samtools_flagstats {
        cpus 2
        errorStrategy 'retry'
        tag { library }
        conda "samtools=1.9"

        input:
            tuple email, library, file(md_file),file(md_bai) from md_files_for_samflagstats.groupTuple(by: [0,1])

        output:
            tuple email, file('*.flagstat') into flagstats
            tuple email, file('*.idxstat') into idxstats
            tuple email, library, file('*.flagstat') into flagstats_for_aggregate
            tuple email, library, file('*.idxstat') into idxstats_for_aggregate
            

        shell:
        '''
        samtools flagstat -@!{task.cpus} !{md_file} > !{md_file}.flagstat
        samtools idxstats !{md_file} > !{md_file}.idxstat
        '''
    }

    process samtools_stats {
        cpus 2
        errorStrategy 'retry'
        tag { library }
        conda "samtools=1.9"

        input:
            tuple email, library, file(md_file),file(md_bai) from md_files_for_samstats.groupTuple(by: [0,1])

        output:

            tuple email, file('*.samstat') into samstats

        shell:
        '''
        samtools stats -@!{task.cpus} !{md_file} > !{md_file}.samstat
        '''
    }

    process picard_gc_bias {
        cpus 1
        errorStrategy 'retry'
        tag { library }
        conda "picard=2.20.7"

        input:
            tuple email, library, file(md_file), file(md_bai) from md_files_for_picard_gc.groupTuple(by: [0,1])

        output:
            tuple email, file('*gc_metrics') into picard_gc_stats

        shell:
        '''
        picard -Xmx4g CollectGcBiasMetrics IS_BISULFITE_SEQUENCED=true VALIDATION_STRINGENCY=LENIENT I=!{md_file} O=!{md_file}.gc_metrics S=!{md_file}.gc_summary_metrics CHART=!{md_file}.gc.pdf R=!{genome}
        '''
    }

    process picard_stats {

        cpus 4
        errorStrategy 'retry'
        tag { library }
        conda "picard=2.20.7"

        input:
            tuple email, library, file(md_file), file(md_bai) from md_files_for_picard.groupTuple(by: [0,1])

        output:
            tuple email, file('*_metrics') into picard_stats

        shell:
        '''
        picard -Xmx16g CollectInsertSizeMetrics VALIDATION_STRINGENCY=LENIENT  I=!{md_file} O=!{md_file}.insertsize_metrics MINIMUM_PCT=0 HISTOGRAM_FILE=/dev/null
        '''
    }

    process human_gc_bias {
        cpus 1
        errorStrategy 'retry'
        tag { library }
        conda "picard=2.20.7"

        input:
            tuple email, library, file(md_file) from human_bams_gc.groupTuple(by: [0,1])

        output:
            tuple email, library, file('*gc_metrics') into human_gc_stats_for_aggregate

        shell:
        '''
        picard -Xmx4g CollectGcBiasMetrics IS_BISULFITE_SEQUENCED=true VALIDATION_STRINGENCY=LENIENT I=!{md_file} O=!{md_file}.gc_metrics S=!{md_file}.gc_summary_metrics CHART=!{md_file}.gc.pdf R=!{genome}
        '''
    }

    process human_insert_size {

        cpus 4
        errorStrategy 'retry'
        tag { library }
        conda "picard=2.20.7"

        input:
            tuple email, library, file(md_file) from human_bams_inserts.groupTuple(by: [0,1])

        output:
            tuple email, library, file('*insertsize_metrics') into human_stats_for_aggregate

        shell:
        '''
        picard -Xmx16g CollectInsertSizeMetrics VALIDATION_STRINGENCY=LENIENT  I=!{md_file} O=!{md_file}.insertsize_metrics MINIMUM_PCT=0.0001 HISTOGRAM_FILE=/dev/null
        '''
    }

    process goleft {
    cpus 1
    conda 'goleft=0.2.0'

    input:
        tuple email, library, file(md_file), file(md_bai) from md_files_for_goleft.groupTuple(by: [0,1])

    output:
        tuple email, file("${library}/*-indexcov.ped") into goleft_ped
        tuple email, file("${library}/*-indexcov.roc") into goleft_roc

    shell:
    '''
        goleft indexcov --directory !{library} *.bam
    '''
    }

    joined_for_multiqc = fastqc_results.groupTuple(by: 0).join(flagstats.groupTuple(by: 0), by: 0)
        .join(idxstats.groupTuple(by: 0), by: 0)
        .join(samstats.groupTuple(by: 0), by: 0)
        .join(picard_stats.groupTuple(by: 0), by: 0)
        .join(picard_gc_stats.groupTuple(by: 0), by: 0)
        .join(goleft_ped.groupTuple(by: 0), by: 0)
        .join(goleft_roc.groupTuple(by: 0), by: 0)
        .join(samblaster_logs.groupTuple(by: 0))
        .join(fastp_log_files.groupTuple(by: 0))

    process multiqc {
        cpus 1
        publishDir "${outputPath}", mode: 'copy'
        conda "multiqc=1.7"

        input:
           tuple email, file(fastqc), file(flagstats), file(idxstats), file(samstats), file(picard_stats), file(picard_gc_stats), file(goleft_ped), file(goleft_roc), file(samblaster), file(fastp) from joined_for_multiqc

        output:
            file "*report.html"

        shell:
        '''
        for file in $(cat input.* | sed -e 's/\\[//g' | sed -e 's/, \\|\\]/\\n/g'); do ln -s ${file} ./; done
        cat <<CONFIG > multiqc_config.yaml 
    title: Bwameth Alignment Summary - !{flowcell}
    extra_fn_clean_exts:
        - '.md'
        - '_combined_fastp'
    custom_plot_config:
        picard_insert_size:
            xmax: 1000
    table_columns_placement:
        Samtools Stats:
            raw_total_sequences: 10
            reads_mapped_percent: 20
            reads_properly_paired_percent: 30
            reads_MQ0_percent: 35
        Samblaster:
            pct_dups: 40
        Picard:
            summed_median: 50
    table_columns_visible:
        Picard:
            PCT_PF_READS_ALIGNED: False
            summed_mean: False
        Samtools Stats:
            reads_mapped: False
            mapped_passed: False
            non-primary_alignments: False
            reads_MQ0_percent: True
        Samtools Flagstat:
            mapped_passed: False
        samtools_idxstats_always:
            - plasmid_puc19c
            - phage_lambda
        FastQC:
            percent_duplicates: False
            total_sequences: False
            avg_sequence_length: False
            percent_fails: False
            total_sequences: False
    CONFIG

        multiqc -ip  .
        '''
    }

    process combine_mbias_tsv {
    publishDir "${outputPath}", mode: 'copy', pattern: 'combined*'

    input:
        tuple email, file(tsv) from mbias_output_tsv.groupTuple(by: 0)

    output:
        file "combined-mbias.tsv"

    shell:
    '''
        echo -ne 'flowcell\tlibrary\t' > combined-mbias.tsv
        ls *_mbias.tsv | head -n 1 | xargs head -n 1 >> combined-mbias.tsv
        for f in *_mbias.tsv; do 
            filebase=`basename "${f}" _combined_mbias.tsv`
            paste <( yes "!{flowcell}	${filebase}" | head -n `nl "$f" | tail -n 1 | cut -f 1` ) "$f" | tail -n +2  >> combined-mbias.tsv
        done
    '''
    }

    process combine_mbias_svg {
    publishDir "${outputPath}", mode: 'copy', pattern: 'combined*'
    conda 'cairosvg=2.4.2 ghostscript=9.22'

    input:
        tuple email, file(svg) from mbias_output_svg.transpose().groupTuple(by: 0)

    output:
        file "combined-mbias.pdf"

    shell:
    '''
    for f in *.svg; do
        cairosvg <(sed s/-nan/-1/ $f) -o $f.pdf
    done
    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=combined-mbias.pdf *.pdf
    '''
    }

