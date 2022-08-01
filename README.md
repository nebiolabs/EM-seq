This repository contains tools and data related to [Enzymatic Methylation Sequencing](https://www.neb.com/products/e7120-nebnext-enzymatic-methyl-seq-kit).

There are 3 [nextflow](https://www.nextflow.io/) scripts:
 - em-seq.nf (to align reads, filter and call methylation)
 - bins.nf (to calculate binned coverage around the TSS
 - cov_vs_meth.nf (to generate the "coverage by feature type" figure from the paper)

A GRCh38 reference genome containing spike-in methylation controls is available via an amazon s3 bucket: s3://neb-em-seq-sra/grch38_core+bs_controls.fa 
or via https: https://neb-em-seq-sra.s3.amazonaws.com/grch38_core%2Bbs_controls.fa

See https://www.biorxiv.org/content/10.1101/2019.12.20.884692v2 for method details and use of these tools.
