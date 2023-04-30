g='mnt/galaxy/data/genome/grch38_hs38d1_noalt+bs_controls/bwameth_index/grch38+bs_controls/grch38_core+bs_controls.fa'

nextflow run em-seq.nf --fastq_glob $1 \
                       --genome $g \
                       --flowcell 'AFLOWCELL'
                       --cpus 4
#'/mnt/galaxy/tmp/users/aerijman@neb.com/*bam' \
