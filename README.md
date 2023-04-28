This repository contains tools and data related to [Enzymatic Methylation Sequencing](https://www.neb.com/products/e7120-nebnext-enzymatic-methyl-seq-kit).

There are 3 [nextflow](https://www.nextflow.io/) scripts:
 - em-seq.nf (to align reads, filter and call methylation)
 - bins.nf (to calculate binned coverage around the TSS
 - cov_vs_meth.nf (to generate the "coverage by feature type" figure from the paper)

A GRCh38 reference genome containing spike-in methylation controls is available via an amazon s3 bucket: s3://neb-em-seq-sra/grch38_core+bs_controls.fa 
or via https: https://neb-em-seq-sra.s3.amazonaws.com/grch38_core%2Bbs_controls.fa

See https://www.biorxiv.org/content/10.1101/2019.12.20.884692v2 for method details and use of these tools.



## Use example:
1. Download the genome (e.g. `wget https://neb-em-seq-sra.s3.amazonaws.com/grch38_core%2Bbs_controls.fa`)
2. If conda is not configured already, you will need [bioconda](https://bioconda.github.io/) for automatic installation of tool dependencies
3. run the nextflow script `nextflow run em-seq.nf --fastq_glob "*{1,2}.fastq" --genome "grch38_core+bs_controls.fa" --flowcell "HCVHLDMXX" --cpus 8`
