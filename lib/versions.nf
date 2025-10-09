def createVersionsFile(versions_topic) {
    versions_topic
      | unique()
      | groupTuple()
      | map{process, names, versions ->
        def pairs = [names, versions].transpose()
        """${process.tokenize(':').last()}:\n${pairs.collect { name, version -> "  ${name}: \"${version}\"" }.join('\n')}\n""".stripIndent()
        }
      | collectFile(name: 'emseq_mqc_versions.yml')
}
