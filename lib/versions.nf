def createVersionsFile(versions_topic) {
    // ESC character (0x1B) — computed to avoid embedding control chars in source.
    // Some tools (e.g. picard, fgbio) emit ANSI colour codes even without a TTY;
    // strip them here so emseq_mqc_versions.yml stays valid YAML.
    def ESC = String.valueOf((char) 27)

    versions_topic
      | unique()
      | groupTuple()
      | map{process, names, versions ->
        def pairs = [names, versions].transpose()
        """${process.tokenize(':').last()}:\n${pairs.collect { name, version ->
            def clean = version
                .replaceAll(ESC + "\\[[0-9;]*[A-Za-z]", "")  // strip CSI sequences: ESC [ params letter
                .replaceAll(ESC, "")                           // strip any remaining bare ESC chars
                .trim()
            "  ${name}: \"${clean}\""
        }.join('\n')}\n""".stripIndent()
        }
      | collectFile(name: 'emseq_mqc_versions.yml')
}
