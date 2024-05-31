

process methylDackel_mbias {
        errorStrategy 'retry'
        tag "${library}"
        conda "bioconda::methyldackel=0.6.1 bioconda::samtools=1.19.2 conda-forge::pigz=2.8"
        publishDir params.outputDir, mode: 'copy'

        input:
            tuple val(library), path(md_bam), path(md_bai), val(barcodes)

        output:
            path('*.svg'), emit: mbias_output_svg
            path('*.tsv'), emit: mbias_output_tsv
            tuple val(library), path('*.tsv'), emit: for_agg

        shell:
        '''
        summarize_combined_mbias.sh !{md_bam} > !{library}_!{barcodes}_combined_mbias_summary.tsv
        
        chrs=( $(samtools view -H !{md_bam} | grep "^@SQ" | sed 's/^.*SN://g' | sed 's/\t.*//g') )

        # makes the svg files for trimming checks with just the first contig
        MethylDackel mbias -@ !{task.cpus} --noCpG --CHH --CHG -r ${chrs[0]} !{params.genome} !{md_bam} !{library}_chn
        for f in *chn*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CHN <\\/text>/" $f; done;

        MethylDackel mbias -@ !{task.cpus} -r ${chrs[0]} !{params.genome} !{md_bam} !{library}_cpg
        for f in *cpg*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CpG<\\/text>/" $f; done;

        '''
    }


process methylDackel_extract {
        tag "${library}"
        publishDir params.outputDir mode: 'copy'
        conda "bioconda::methyldackel=0.6.1 bioconda::samtools=1.19.2 conda-forge::pigz=2.8"
        input:
            tuple val(library), path(md_bam), path(md_bai), val(barcodes) 

        output:
            tuple val(library), file('*.methylKit.gz'), emit: extract_output 

        shell:
        '''
        MethylDackel extract --methylKit -q 20 --nOT 0,0,0,5 --nOB 0,0,5,0 -@ !{task.cpus} --CHH --CHG -o !{library}.!{barcodes} !{params.genome} !{md_bam} 
        pigz -p !{task.cpus} *.methylKit 
        '''
    }