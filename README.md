# EM-seq Analysis Pipeline

[![Test Status](https://github.com/nebiolabs/EM-seq/actions/workflows/test.yml/badge.svg)](https://github.com/nebiolabs/EM-seq/actions)

This repository contains Nextflow-based analysis tools for [Enzymatic Methylation Sequencing (EM-seq)](https://www.neb.com/products/e7120-nebnext-enzymatic-methyl-seq-kit) and [Enzymatic 5hmC-seq (E5hmC-seq)](https://www.neb.com/en-us/products/e3350nebnext-enzymatic-methyl-seq-5hmc-kit) data processing.

### Main Analysis Pipeline (`main.nf`)
Complete EM-seq processing pipeline that accepts UBAM inputs:
- Adapter trimming and read alignment with (fastp, bwa-meth)
- Duplicate marking (Picard)
- Methylation calling (MethylDackel)
- Quality control metrics and statistics (Picard, Samtools, FastQC, MultiQC)
- Optional BED file intersection for targeted analysis (bedtools)

### Fastq to uBam pipeline (`fastq_to_ubam.nf`)
If your files are in fastq format you will need to convert them to uBams prior to running the main pipeline, e.g.:
```bash
nextflow run fastq_to_ubam.nf \
  --input_glob "tests/fixtures/fastq/emseq-test*{.ds.1,.ds.2}.fastq.gz" \
  --read_format 'paired-end'
```
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input_glob` | glob for your gzipped fastq files | `['*.{1,2}.fastq.gz']` |
| `--read_format` | 'paired-end' or 'single-end' | `'paired-end'` |

## Quick Start
1. Install [miniforge](https://conda-forge.org/download/) and [bioconda](https://bioconda.github.io/) (see [Requirements](<README#Requirements>))
2. Install [Nextflow](https://www.nextflow.io/) (e.g. conda install nextflow, or see [Nextflow installation guide](https://www.nextflow.io/docs/latest/getstarted.html#installation))
3. Clone this repository (`git clone https://github.com/nebiolabs/EM-seq.git`). Copy `nextflow.config.example` to `nextflow.config` and modify default settings as needed for your environment
4. Download or prepare a genome reference FASTA file (see [Reference Genomes](<README#Reference Genomes>))
5. Create a bwameth index for the fasta and configure your references in conf/references.config
6. Run the pipeline with appropriate parameters (see [Basic Usage](<README#Basic Usage>))
7. Examine results in the EM-seq_output directory
   - `EM-seq-Alignment-Summary-<FLOWCELL_ID>_multiqc_report.html` in em-seq_output for overall QC summary
   - Mbias files `em-seq_output/methylDackelExtracts/mbias` (to identify sample-dependent positional biases)
   - Methylation output files in `em-seq_output/methylDackelExtracts` (suitable for analysis with [methylKit](https://bioconductor.org/packages/release/bioc/html/methylKit.html))
   - Aligned reads in `em-seq_output/markduped_bams` (methylation coloring is recommended for visualization in [IGV](https://igv.org/doc/desktop/#UserGuide/tracks/alignments/bisulfite_sequencing/))

### Basic Usage
```bash
nextflow run main.nf \
  --genome 'test' \
  --ubam_dir './' \
  --email your.email@example.com \
  --flowcell FLOWCELL_ID
```
`ubam_dir` should be the folder where your ubam files are.

### Key Parameters
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--genome` | reference genome found in conf/references.config | Required |
| `--email` | Email for notifications | Required |
| `--flowcell` | Flowcell identifier | Optional |
| `--outputDir` | Output directory | `em-seq_output` |
| `--project` | Project name | `project_undefined` |
| `--enable_neb_agg` | Enable NEB aggregation reporting | `False` |

### References Config

Modify the conf/references.config file to specify your genome files
- `genome_fa` path to your genome fasta file 
- `genome_fai` path to your genome fasta fai file 
- `bwameth_index` path to your genome fasta file where bwameth indices exist 
- `target_bed` BED file for targeted analysis, Optional 

### Advanced Options
- `--tmp_dir` - Temporary directory (default: `/tmp`)
- `--workflow` - Workflow identifier (default: `EM-seq`)
- `--enable_neb_agg` - Enable NEB aggregation reporting (default: `False`)


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
These in the "legacy" folder are retained for reference and reproducibility but are not actively maintained and are not compatible with the latest Nextflow versions. Use `NXF_VER=22.10.4 nextflow run ...` to reproduce the results in the [EM-seq paper](README#Citation).
- `em-seq.nf` - Original alignment and methylation calling workflow
- `bins.nf` - TSS-centered binned coverage analysis
- `cov_vs_meth.nf` - Coverage vs methylation analysis for genomic features

## Citation
Analysis methods in this repository were used in the following publication:

Vaisvila R, Ponnaluri VKC, Sun Z, et al. **Enzymatic methyl sequencing detects DNA methylation at single-base resolution from picograms of DNA.** *Genome Res.* 2021;31(7):1280-1289. doi:[10.1101/gr.266551.120](https://doi.org/10.1101/gr.266551.120)

## Related Projects
You may also be interested in the [nf-core methylseq project](https://nf-co.re/methylseq/2.5.0)

## Developer documentation
### Production:
 - git tag -f current_production
 - git push -f origin current_production

 ### Development:
 - development workflow will run from master branch

 ### Testing:
 - Tests are run using nf-test and are integrated into github actions
 - install nf-test from bioconda using conda/mamba
 - To run all tests:
 ```bash
nf-test test
```
- When new tests are added or results change, to update the results snapshot:
```bash
nf-test test --updateSnapshot
```

