
library = params.library

process methylDackel_mbias {
        label 'cpus_8'
        errorStrategy 'retry'
        tag {library}
        conda "methyldackel=0.4.0 samtools=1.9"
        publishDir "${library}/methylDackelExtracts/mbias"

        input:
            tuple val(library), path(md_file), path(md_bai) 
            // from md_files_for_mbias.groupTuple()

        output:
            path('*.svg'), emit: mbias_output_svg
            path('*.tsv'), emit: mbias_output_tsv
            tuple val(library), path('*.tsv'), emit: mbias_for_aggregate

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
                    MethylDackel mbias --noSVG $arg -@ !{task.cpus} -r $chr !{params.genome} !{md_file} | \
                    tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
                ) \
                <( \
                    MethylDackel mbias --noSVG --keepDupes -F 2816 $arg -@ !{task.cpus} -r $chr !{params.genome} !{md_file} | \
                    tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
                ) \
                | sed "s/^/${chr}\t${context}\t/" \
                >> !{library}_combined_mbias.tsv
            done
        done
        # makes the svg files for trimming checks
        MethylDackel mbias -@ !{task.cpus} --noCpG --CHH --CHG -r ${chrs[0]} !{params.genome} !{md_file} !{library}_chn
        for f in *chn*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CHN <\\/text>/" $f; done;

        MethylDackel mbias -@ !{task.cpus} -r ${chrs[0]} !{params.genome} !{md_file} !{library}_cpg
        for f in *cpg*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CpG<\\/text>/" $f; done;

        '''

    }


process methylDackel_extract {
        label 'cpus_8'
        tag library
        publishDir "${library}/methylDackelExtracts", mode: 'copy'
        conda "methyldackel=0.4.0 pigz=2.4"

        input:
            tuple val(library), path(md_file), path(md_bai) // from md_files_for_extract.groupTuple()

        output:
            tuple val(library), file('*.methylKit.gz'), emit: extract_output //into extract_output

        shell:
        '''
        MethylDackel extract --methylKit --nOT 0,0,0,5 --nOB 0,0,5,0 -@ !{task.cpus} --CHH --CHG -o !{library} !{params.genome} !{md_file}
        pigz -p !{task.cpus} *.methylKit
        '''

    }