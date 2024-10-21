default_dest_path = '/mnt/galaxy/tmp/users'

/* Are the followint lines to keep a single repo for both projects? 
 * Don't we use different deployments for these purposes?
 * e.g. development and production.
 */ 
// if ("$workflow.projectDir".contains('dev')) {
//     dest_modifier = '-dev'
//    seq_shepherd_dir = '/mnt/bioinfo/prg/seq-shepherd-dev/current'
// } else {
//     dest_modifier = ''
    seq_shepherd_dir = '/mnt/bioinfo/prg/seq-shepherd/current'
//}


tiles = Channel.fromPath(params.path_to_tiles + '/L_*_*', type: 'dir')
read_counts = Channel.fromPath(params.path_to_tiles + '/L*_metrics.txt')
sample_sheet_path = Channel.fromPath(params.sample_sheet)
downsampling_params = Channel.fromPath(params.path_to_tiles + '/../*_downsampling.txt')

process parallel_yamls { 
    cpus 1

    conda 'samtools=1.9'

    input:
        val(tile)

    output:
        file '*.yaml'

    shell:
    '''
    !{seq_shepherd_dir}/post_run_scripts/generate_tile_yaml.rb --path_to_bams !{tile} --email !{params.email} 
    '''
}
process merge_tiles {
    cpus 1

    conda 'conda-forge::pyyaml=5.4.1'

    input:
        path yaml
        path reads
        path sample_sheet
        path downsampling
    
    output:
        path('merged.yaml')            

    shell:
    '''
    cat *.yaml | python3 !{seq_shepherd_dir}/post_run_scripts/merge_yaml.py --email !{params.email} > merged.yaml
    '''
}

process call_again {
    cpus 1
    input:
        path(yaml)

    shell:
    '''
    flowcell=$(grep flowcell !{yaml} | cut -d " " -f 2)
    nextflow run !{seq_shepherd_dir}/post_run_scripts/em-seq.nf -resume \
    -w /mnt/hpc_scratch/novaseq_processing/${flowcell} \
    -params-file !{yaml} \
    -with-report  "!{params.run_path}/emseq_workflow_report.html" \
    -with-timeline "!{params.run_path}/emseq_workflow_timeline.html" \
    -with-dag "!{params.run_path}/emseq_workflow_dag.html" \
    -with-weblog !{params.run_uri} \
    -N !{params.email},galaxyadmin@nebiolabs.onmicrosoft.com
    '''
}


workflow path_to_tiles_known {
    take: tiles,
          read_counts,
          sample_sheet_path,
          downsampling_params
    main:
        parallel_yamls_channel = parallel_yamls(tiles.tile)
        for_recalling_channel = merge_tiles(all_yamls.collect.yaml,
                                            read_counts.collect.reads,
                                            sample_sheet_path.sample_sheet,
                                            downsampling_params.downsampling

        )
        call_again(for_recalling_channel.yaml)
}