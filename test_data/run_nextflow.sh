g='mnt/galaxy/data/genome/grch38_hs38d1_noalt+bs_controls/bwameth_index/grch38+bs_controls/grch38_core+bs_controls.fa'

nextflow run main.nf --input_files '/mnt/galaxy/tmp/users/aerijman@neb.com/*bam' \
                     --genome $g \
                     --flowcell 'AFLOWCELL'
