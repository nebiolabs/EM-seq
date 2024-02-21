def inferReadLength(inputGlob) {
    def readLengthScript
    if (inputGlob.contains('.bam')) {
        readLengthScript = "samtools view \$(ls ${inputGlob} | head -n 1) | head -n 10000 | awk -F'\t' '{if (length(\$10) > longest_read) {longest_read = length(\$10)}} END {print longest_read}'"
    } else { // assumes fastq if not bam
        readLengthScript = "zcat -f \$(ls ${inputGlob} | head -n 1) | head -n 40000 | awk \'NR%4==2 {if (length(\$1) > longest_read) {longest_read = length(\$1)}} END {print longest_read}'"
    }
    println "trying to infer read length using: ${readLengthScript}" 
    def readLength = ['/bin/sh', '-c', readLengthScript].execute().text.trim().toInteger()
    println "Detected read length: ${readLength}"
    return readLength
}