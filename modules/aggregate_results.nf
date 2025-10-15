def format_ngs_agg_opts(optlist) {
    '''
    optlist: list of lists of format:
        [
            [ opt1_name, Channel( val(library), results1 ) ],
            [ opt2_name, Channel( val(library), results2 ) ],
        ]

    Returns a channel of format: [lib, [opt1, opt2], [results1, results2]]
    '''

    // Initialize agg channel with first opt/value pair
    agg_tuple = Channel.of(optlist[0][0]).combine(optlist[0][1]).map{opt, lib, val -> [lib, [opt, val]]}
    // Add remaining opt/value pairs
    for(item in optlist[1..-1]) {
        lib_opt_val = Channel.of(item[0]).combine(item[1]).map{opt, lib, val -> [lib, [opt, val]]}
        agg_tuple = agg_tuple.join(lib_opt_val)
    }
    // Transpose to change from (lib, [opt1, val1], [opt2, val2]) to (lib, [opt1, opt2], [val1, val2])
    agg_tuple = agg_tuple.map{[it[0], *it[1..-1].transpose()]}

    return agg_tuple
}


process aggregate_results {
    tag { library }
    conda "bioconda::samtools=1.9 git=2.40.1"
    publishDir "${params.outputDir}/aggregate_results", mode: 'copy'
    label 'process_single'

    input:
	tuple val(library), val(ngs_agg_opts), path(ngs_agg_paths)
    val(workflow_name)

    output:
        path("*arguments.txt")

    script:
    def opt_val = [ngs_agg_opts, ngs_agg_paths].transpose().collect{ opt, fp -> "${opt} ${fp}" }.join(' ')
    opt_val = opt_val.replaceFirst(/fastqc.zip/, "fastqc/fastqc_data.txt")
    """
    echo "${opt_val}" | sed 's/--/\\n--/g' > ${library}.arguments.txt
    unzip -o *fastqc.zip

    export GIT_HASH=\$(git -C "${workflow.projectDir}" log -1 --pretty=format:"%H")

    export RBENV_VERSION=\$(cat ${params.path_to_ngs_agg}/.ruby-version)
    RAILS_ENV=production \
    ${params.path_to_ngs_agg}/bin/bundle exec \
    ${params.path_to_ngs_agg}/aggregate_results.rb \\
        --workflow "${workflow_name}" \\
        --commit_hash \$GIT_HASH \\
        --contact_email "${params.email}" \\
        ${opt_val} \\

    """
    stub:
    def opt_val = [ngs_agg_opts, ngs_agg_paths].transpose().collect{ opt, fp -> "${opt} ${fp}" }.join(' ')
    """
    echo "${opt_val}" | sed 's/--/\\n--/g' > ${library}.arguments.txt
    echo "--workflow ${workflow_name}" >> ${library}.arguments.txt
    """
}
