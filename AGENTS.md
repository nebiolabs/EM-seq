# Copilot Instructions for EM-seq Repository

## Project Overview
This repository contains Nextflow-based analysis tools for Enzymatic Methylation Sequencing (EM-seq) and Enzymatic 5hmC-seq (E5hmC-seq) data processing. The pipeline handles adapter trimming, read alignment, duplicate marking, methylation calling, quality control metrics, and optional BED file intersection for targeted analysis.

## Repository Structure

### Main Pipelines
- `main.nf` - Main EM-seq processing pipeline (accepts UBAM inputs)
- `fastq_to_ubam.nf` - Converts fastq files to uBam format before running main pipeline

### Key Directories
- `modules/` - Modular Nextflow processes for various pipeline steps (fastp, alignment, methylation calling, QC, etc.)
- `conf/` - Configuration files including base settings and reference genomes
- `lib/` - Shared library code (versions tracking)
- `tests/` - nf-test test suite with fixtures and snapshots
- `assets/` - Supporting files like methylation control sequences
- `legacy/` - **IMPORTANT: DO NOT MODIFY** - Historical workflows retained for reference and reproducibility (see below)

### Legacy Folder
The `legacy/` folder contains historical workflows that are **preserved as-is** and should **never be modified**:
- `em-seq.nf` - Original alignment and methylation calling workflow
- `bins.nf` - TSS-centered binned coverage analysis
- `cov_vs_meth.nf` - Coverage vs methylation analysis
- Associated scripts and configuration files

These workflows are retained for reproducibility of results from the EM-seq paper and are not compatible with latest Nextflow versions. They require `NXF_VER=22.10.4` to run.

## Development Guidelines

### Code Style
- Follow existing Nextflow DSL2 patterns used in the repository
- Use modular process definitions in the `modules/` directory
- Maintain consistent naming conventions (snake_case for files, camelCase for Nextflow processes)
- Include appropriate comments for complex logic

### Testing
- Tests use nf-test framework (installed via bioconda)
- Test files: `tests/main.nf.test`, `tests/fastq_to_ubam.nf.test`
- Run tests: `nf-test test`
- Update snapshots when results change: `nf-test test --updateSnapshot`
- Tests run automatically in GitHub Actions (`.github/workflows/test.yml`)

### Configuration
- Genome references are configured in `conf/references.config`
- Base execution settings in `conf/base.config`
- Test-specific settings in `conf/test.config`
- Main configuration in `nextflow.config`

### Dependencies
- Use Conda/Mamba/Micromamba for dependency management
- Bioconda channel is required
- Nextflow is the primary workflow engine

### Key Parameters
When working with pipeline parameters, be aware of:
- `--genome` - Reference genome from conf/references.config (required)
- `--ubam_dir` - Input directory for uBAM files
- `--outputDir` - Output directory (default: `em-seq_output`)
- `--email` - Email for notifications (required)
- `--flowcell` - Flowcell identifier
- `--enable_neb_agg` - Enable NEB aggregation reporting

## Best Practices

### When Making Changes
1. Understand the modular structure - each process should be self-contained in `modules/`
2. Test changes using the nf-test framework
3. Update snapshots if expected outputs change
4. Ensure compatibility with both conda and other container/environment managers
5. Document any new parameters or significant changes in README.md

### Common Operations
- **Adding new module**: Create in `modules/` directory, follow existing pattern
- **Modifying workflow**: Edit `main.nf` or `fastq_to_ubam.nf`
- **Updating references**: Modify `conf/references.config`
- **Adding tests**: Update test files in `tests/` and regenerate snapshots

### Quality Control
- Pipeline generates MultiQC reports for overall QC
- Individual QC metrics from Picard, Samtools, FastQC
- Mbias plots for identifying positional biases
- Ensure new features integrate with existing QC reporting

## Citation and Context
Methods used in this repository were published in:
Vaisvila R, et al. "Enzymatic methyl sequencing detects DNA methylation at single-base resolution from picograms of DNA." *Genome Res.* 2021;31(7):1280-1289.

## Related Resources
- NEB EM-seq Kit: https://www.neb.com/products/e7120-nebnext-enzymatic-methyl-seq-kit
- Related project: nf-core methylseq (https://nf-co.re/methylseq)
- Pre-built references available for T2T CHM13 and GRCh38 with methylation controls
