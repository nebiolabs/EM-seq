# EM-seq Analysis Pipeline

This repository contains Nextflow-based analysis tools for [Enzymatic Methylation Sequencing (EM-seq)](https://www.neb.com/products/e7120-nebnext-enzymatic-methyl-seq-kit) and [Enzymatic 5hmC-seq (E5hmC-seq)](https://www.neb.com/en-us/products/e3350nebnext-enzymatic-methyl-seq-5hmc-kit) data processing.

## Workflows

The repository includes multiple [Nextflow](https://www.nextflow.io/) workflows:

### Main Analysis Pipeline (`main.nf`)
Complete EM-seq processing pipeline that accepts both FASTQ and BAM inputs:
- Read alignment with BWA-meth
- Duplicate marking and filtering  
- Methylation calling with MethylDackel
- Quality control metrics and statistics
- Optional BED file intersection for targeted analysis

### Legacy Workflows
- `em-seq.nf` - Original alignment and methylation calling workflow
- `bins.nf` - TSS-centered binned coverage analysis
- `cov_vs_meth.nf` - Coverage vs methylation analysis for genomic features

## Quick Start

### Basic Usage
```bash
nextflow run main.nf \
  --path_to_genome_fasta /path/to/genome.fa \
  --input_glob "*_R1.fastq*" \
  --email your.email@example.com \
  --flowcell FLOWCELL_ID
```
`input_glob` should capture your bam files or `read1` fastq files  (even for paired-end libraries)

### Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--path_to_genome_fasta` | Path to reference genome FASTA file | Required |
| `--input_glob` | Glob pattern for input files (FASTQ/BAM) | `*_R1.fastq*`/`*.bam` |
| `--email` | Email for notifications | Required |
| `--flowcell` | Flowcell identifier | Optional |
| `--outputDir` | Output directory | `em-seq_output` |
| `--project` | Project name | `project_undefined` |
| `--max_input_reads` | Max reads to process (or "all_reads") | `all_reads` |
| `--target_bed` | BED file for targeted analysis | Optional |
| `--enable_neb_agg` | Enable NEB aggregation reporting | `False` |
| `--min_mapq` | Minimum mapping quality | `20` |

### Advanced Options
- `--downsample_seed` - Random seed for downsampling (default: 42)
- `--tmp_dir` - Temporary directory (default: `/tmp`)
- `--workflow` - Workflow identifier (default: `EM-seq`)

## Reference Genomes

Pre-built reference genomes with methylation spike-in controls:
- **GRCh38**: https://neb-em-seq-sra.s3.amazonaws.com/grch38_core%2Bbs_controls.fa
- **T2T CHM13**: https://neb-em-seq-sra.s3.amazonaws.com/T2T_chm13v2.0%2Bbs_controls.fa

## Requirements

- Nextflow (version 22.10.4 recommended for legacy scripts)
- Conda/Mamba for dependency management
- Sufficient computational resources (memory scales with input size)

## Configuration

Copy `nextflow.config.example` to `nextflow.config` and modify as needed for your environment.

## Related Projects

You may also be interested in the [nf-core methylseq project](https://nf-co.re/methylseq/2.5.0)
