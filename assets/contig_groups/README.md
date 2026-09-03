# Contig group definitions

One TSV per reference, assigning every contig to a group so the pipeline can report a GC bias
curve per organism instead of a single curve blended across the whole composite. Setting a
genome's `gc_groups_dir` in `conf/references.config` turns the feature on; see **GC bias curves**
in the top-level README.

These cover the references linked under [Reference
Genomes](../../README.md#reference-genomes), and are checked against the genome's `.fai` at
startup, so they only apply to those exact FASTAs.

| Column | Meaning |
|---|---|
| `contig_name` | must match the reference `.fai` exactly |
| `sequence_length` | must match the `.fai`; a mismatch means the TSV was built against a different reference |
| `frac_gc` | GC over ACGT positions only, so an N-heavy contig is not reported as low-GC. Informational — the pipeline does not read it |
| `contig_group` | the curve this contig contributes to. `[A-Za-z0-9_.-]+`. Blank means no curve of its own |
| `contig_type` | one of `autosomal`, `sex_chromosome`, `organelle`, `unplaced`, `control`, `unknown` |

Groups are named `<organism>_autosome`, `<organism>_sex`, `<organism>_mito` and
`<organism>_unplaced` for the host, and `lambda`, `pUC19`, `T4`, `Xp12` (plus `EBV` for the human
references) for the non-host contigs:

| Reference | Groups |
|---|---|
| `T2T_chm13v2.0+meth_controls` | human_autosome (22), human_sex (2), human_mito, EBV, lambda, pUC19, T4, Xp12 |
| `grch38_core+meth_controls` | human_autosome (63), human_sex (3), human_mito, human_unplaced (127), EBV, lambda, pUC19, T4, Xp12 |
| `grcm39+meth_controls` | mouse_autosome (58), mouse_sex (2), mouse_mito, lambda, pUC19, T4, Xp12 |

## Provenance

These are the group assignments NEB uses internally, and match published references.

One caveat carried over from how the references were built: **`grcm39+meth_controls` names its
chromosomes by GenBank accession** (`CM000994.3` for chr1, and so on) rather than `chr1`-style,
unlike the two human references. This is deliberate — it is why `feature_cov_meth.nf`'s mouse
preset sets `cpg_chr_lookup = '$10,$5'` to translate assembly-report chromosome names to
accessions, where the T2T preset uses `'$10,$10'`.

## Adding another reference

Write a TSV with the five columns above: one row per contig in the reference's `.fai`, in `.fai`
order, with `sequence_length` matching the `.fai` exactly. Give every contig that should get its own
curve a `contig_group`, and leave the column blank for contigs that should not.

Then build one subset reference per group and point that genome's `gc_groups_dir` at the directory
holding them:

```bash
TSV=contig_groups.tsv
FA=genome.fa            # must already have a .fai alongside it

for grp in $(awk -F'\t' 'NR>1 && $4!=""{print $4}' "$TSV" | sort -u); do
    # contig list and BED are taken in .fai order, which is the order picard requires
    awk -F'\t' -v g="$grp" 'NR==FNR{if(FNR>1&&$4==g)k[$1];next} $1 in k{print $1}' \
        "$TSV" "$FA.fai" > "$grp.contigs.txt"
    awk -F'\t' -v g="$grp" 'NR==FNR{if(FNR>1&&$4==g)k[$1];next} $1 in k{print $1"\t0\t"$2}' \
        "$TSV" "$FA.fai" > "$grp.bed"
    samtools faidx "$FA" -r "$grp.contigs.txt" > "$grp.fa"
    samtools faidx "$grp.fa"
    samtools dict  "$grp.fa" > "$grp.dict"
done
```

Each group's subset FASTA is a copy of those contigs' sequence, so the directory costs roughly as
much disk as the reference itself. The pipeline validates the TSV against `genome_fai` and checks
each group's contig order against its subset reference before running picard, so a mistake here
fails at startup rather than part-way through a run.
