process test_flagstats {
  input:
    tuple val(library), path (flagstats_file)

  output:
    stdout

  script:
    """
    #!/bin/bash
    set -e

    if grep -q "1972 + 0 properly paired" ${flagstats_file}; then
        echo "flagstats OK"
    else
        echo "flagstats not OK"
        return 1
    fi
    """
}


process test_alignment_metrics {
  input:
    tuple val(library), path (picard_alignment_summary_metrics_file)

  output:
    stdout

  script:
    """
    #!/bin/bash
    set -e

    alignment_result=\$(tail -n2 ${picard_alignment_summary_metrics_file} | \
        awk 'BEGIN{result="alignment metrics not OK"}{if (\$1==150 && \$3>2200) {result="alignment metrics OK"}}END{print result}')
    echo "\${alignment_result}" 

    if [[ "\${alignment_result}" != "alignment metrics OK" ]]; then
        return 1
    fi
    """
}
