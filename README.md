# EM-seq Analysis Pipeline

This repository contains Nextflow-based analysis tools for [Enzymatic Methylation Sequencing (EM-seq)](https://www.neb.com/products/e7120-nebnext-enzymatic-methyl-seq-kit) and [Enzymatic 5hmC-seq (E5hmC-seq)](https://www.neb.com/en-us/products/e3350nebnext-enzymatic-methyl-seq-5hmc-kit) data processing.

### Main Analysis Pipeline (`main.nf`)
Complete EM-seq processing pipeline that accepts both FASTQ and BAM inputs:
- Adapter trimming and read alignment with (fastp, bwa-meth)
- Duplicate marking (Picard)
- Methylation calling (MethylDackel)
- Quality control metrics and statistics (Picard, Samtools, FastQC, MultiQC)
- Optional BED file intersection for targeted analysis (bedtools)

## Quick Start
1. Install [miniforge](https://conda-forge.org/download/) and [bioconda](https://bioconda.github.io/) (see [Requirements](<README#Requirements>))
1. Install [Nextflow](https://www.nextflow.io/) (e.g. conda install nextflow, or see [Nextflow installation guide](https://www.nextflow.io/docs/latest/getstarted.html#installation))
1. Clone this repository (git clone https://github.com/nebiolabs/EM-seq.git), Copy `nextflow.config.example` to `nextflow.config` and modify default settings as needed for your environment
1. Donwload or prepare a genome reference FASTA file (see [Reference Genomes](<README#Reference Genomes>))
1. Run the pipeline with appropriate parameters (see [Basic Usage](<README#Basic Usage>))
1. Examine results in the EM-seq_output directory
   - `EM-seq-Alignment-Summary-<FLOWCELL_ID>_multiqc_report.html` in em-seq_output for overall QC summary
   - Mbias files `em-seq_output/methylation/mbias` (to identify sample-dependent positional biases)
   - Methylation output files in `em-seq_output/methylation` (suitable for analysis with [methylKit](https://bioconductor.org/packages/release/bioc/html/methylKit.html))
   - Aligned reads in `em-seq_output/bams` (methylation coloring is recommended for visualization in [IGV](https://igv.org/doc/desktop/#UserGuide/tracks/alignments/bisulfite_sequencing/))

### Basic Usage
```bash
nextflow run main.nf \
  --path_to_genome_fasta /path/to/genome.fa \
  --input_glob "*_R1.fastq*" \
  --email your.email@example.com \
  --flowcell FLOWCELL_ID
```
`input_glob` should capture your bam files (aligned or unaliged) or `read1` fastq files (read2 will be inferred for paired-end libraries). Be sure to quote the glob pattern to avoid shell expansion. Test your glob pattern with `ls` (without the quotes) to confirm it matches the expected files.

### Key Parameters
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--path_to_genome_fasta` | Path to reference genome FASTA file | Required |
| `--input_glob` | Glob pattern for input files (FASTQ read 1/BAM) | `*_R1.fastq*`/`*.bam` |
| `--email` | Email for notifications | Required |
| `--output_dir` | Output directory | `em-seq_output` |
| `--max_input_reads` | Max reads to process (or "all_reads") | `all_reads` |
| `--min_mapq` | Minimum mapping quality | `20` |
| `--trim_bases` | Number of bases to trim (correcting end repair) | `5` |

### Advanced Options
- `--flowcell`- Flowcell identifier for reports (Optional)
- `--target_bed` - BED file for supplemental reporting of methylation in  specfic regions - (Optional)
- `--downsample_seed` - Random seed for downsampling (default: 42)
- `--tmp_dir` - Temporary directory (default: `/tmp`)
- `--workflow` - Workflow identifier (default: `EM-seq`)
- `--enable_neb_agg` - Enable NEB aggregation reporting (default: `False`)
- `--publish_intermediate_output` - Publish intermediate files to the output directory (default: `False`)
- `--publish_dir_mode` - Nextflow [publishDir mode](https://www.nextflow.io/docs/latest/reference/process.html#mode) for published files (default: `symlink`)


## Reference Genomes
Pre-built reference genomes with methylation spike-in controls:
- **T2T CHM13**: https://neb-em-seq-sra.s3.amazonaws.com/T2T_chm13v2.0%2Bbs_controls.fa
- **GRCh38**: https://neb-em-seq-sra.s3.amazonaws.com/grch38_core%2Bbs_controls.fa
- Create your own reference by appending the [control sequences](methylation_controls.fa) to your preferred genome fasta (e.g. `cat genome.fa methylation_controls.fa > genome+methylation_controls.fa`)

## Requirements
- [Nextflow](https://www.nextflow.io/)
- [Miniforge](https://conda-forge.org/download/), [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html), or [Conda](https://docs.conda.io/projects/conda/en/stable/) for dependency management
- [Bioconda](https://bioconda.github.io/) channel configured
- Sufficient computational resources (memory scales with input size)

### Historical Workflows
These legacy scripts are retained for reference and reproducibility but are not actively maintained and are not compatible with the latest Nextflow versions. Use `NXF_VER=22.10.4 nextflow run ...` to reproduce the reuslts in the [EM-seq paper](README#Citation).
- `em-seq.nf` - Original alignment and methylation calling workflow from
- `bins.nf` - TSS-centered binned coverage analysis
- `cov_vs_meth.nf` - Coverage vs methylation analysis for genomic features

## Citation
Analysis methods in this repository were used in the following publication:

Vaisvila R, Ponnaluri VKC, Sun Z, et al. **Enzymatic methyl sequencing detects DNA methylation at single-base resolution from picograms of DNA.** *Genome Res.* 2021;31(7):1280-1289. doi:[10.1101/gr.266551.120](https://doi.org/10.1101/gr.266551.120)
`
### Related Projects
You may also be interested in the [nf-core methylseq project](https://nf-co.re/methylseq/2.5.0)
