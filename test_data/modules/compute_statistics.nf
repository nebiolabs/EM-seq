
    process picard_gc_bias {
        cpus 1
        tag { library }
        conda "picard=2.20.7"

        input:
            tuple val(library), file(bam), file(bai), val(barcodes)

        output:
            tuple val(library), path('*gc_metrics'), emit: gc_metrics

        shell:
        '''
        prefix=$(echo $bam | sed 's/bam//')
        picard -Xmx4g CollectGcBiasMetrics IS_BISULFITE_SEQUENCED=true VALIDATION_STRINGENCY=SILENT I=!{bam} O=${prefix}.gc_metrics S=${prefix}.gc_summary_metrics CHART=${prefix}.gc.pdf R=!{params.genome}
        '''
    }