# Benchmark results: bamadap + fgumi zipper + picodup vs the current pipeline

Performance changes in branch: `bamslice_zip`

**Bottom line**: ~2.4x faster end-to-end on a clean dedicated node (16m/6m45s), with indistinguishable output. Correctness verified in depth (§2);

## 1. What changed, and the measured impact

| Change | Impact |
| --- | --- |
| `fastp`+`bwa`+Picard-dedup pipeline → `bamadap` + `fgumi zipper` + `picodup`, fused into single-pipe `trimAndAlign` + streamed `mergeAndPicodup` | End-to-end: **16m → 6m45s (~2.4x)** |
| Dedup: Picard `MarkDuplicates` → `picodup`, streamed straight off `samtools merge` (no intermediate merged BAM) | **2m02s → 27s (~4.6x)**, peak RSS 8.5GB → <1GB |
| Fixed a real `bwameth.py` bug (naive interleave-detection failed on `bamadap`'s `read/1`,`read/2` mate names, silently mis-converting bisulfite reads); fixed via `bamadap --no-mate-suffix` rather than a downstream patch | Removes a 3x slowdown *and* a silent correctness bug in the align stage |
| `convert_methylkit_to_bed`: single-threaded `gawk` → `mawk` under `parallel --pipepart` | **7m53s → 35s (~14x)** — biggest single bottleneck in either pipeline; benefits `master` too |
| `tasmanian` 1.x (single-threaded, capped at 2M reads, errored on some real reads) → `tasmanian-mismatch` 2.x (parallel, full library) | **3m42s → 17.5s (~13x), on more data** |
| `fastqc` → `falco` | 57.8s → 30.8s (~2x) |
| Fixed `gc_bias` (silently broken: `CollectGcBiasMetrics` requires an `Rscript` on PATH even when no chart is produced; `picard-slim` excludes R) | No-op `Rscript` shim; restores a previously-silent-failing QC step |
| Right-sized cpu reservations for QC steps that don't read `task.cpus` (new `single_threaded_qc` label, fixed at 2 cpus instead of scaling with `--max_cpus`) | Frees queue slots on shared executors; no effect on single-task time |
| Kept `bwameth --read-group` (populated from the uBAM's first `@RG`) instead of dropping it as originally planned; `zipper` still fixes the per-read `RG:Z:` tag | Avoids `picodup` mislabeling every metric row "Unknown Library" |
| Fixed `bamadap`'s EL8 build (glibc-2.34-only symbols, `target-cpu=native` SIGILL'd on other CPUs) | Portability now solved at deploy time, not build time: Capistrano (`capistrano-rust-buildcache`) compiles a `-C target-cpu=native` binary per distinct machine type on first run and caches it, rather than shipping one generic binary |

Chunk-size sweep on the real SGE cluster (18/37/55/90MB): smaller chunks buy wall-clock speed at
the cost of more aggregate CPU-hours (7m14s/7.1 CPU-h at 18MB vs. 10m18s/4.4 CPU-h at 90MB); no
change made to the 37MB default — it's a reasonable middle ground, not a correctness question.

## 2. Correctness (final `.md.bam`, same test uBAM)

| | Baseline | New | Δ |
| --- | --- | --- | --- |
| Total reads | 10,017,473 | 10,017,757 | +0.003% |
| Duplicates | 1,362,823 | 1,362,959 | +0.01% |
| PERCENT_DUPLICATION | 0.136415 | 0.136426 | matches to 4th decimal |
| CpG methylation (Pearson r) | — | **0.999998** | 5 sites (0.0001%) differ by >1pp |

Residual differences trace to expected causes (`bamadap` vs `fastp` trimming, corrected
per-read `RG` moving optical-duplicate grouping, each dedup tool's own library-size extrapolation
formula) — no read loss, no regression. `ngs-aggregate_results` compatibility for
`tasmanian-mismatch` 2.x's new output format was verified against the real deployed parser
(PR #932, merged and deployed 2026-09-07) and an end-to-end real (non-stubbed) aggregation run
succeeded.
