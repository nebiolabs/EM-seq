// GC bias curves for one set of aligned BAMs: always the whole-reference curve, plus one curve per
// organism when the genome is a composite that ships a gc_groups_dir.
//
// Encapsulated as a subworkflow so main.nf sees one call regardless of whether the genome has
// contig groups, and so the conditional wiring and the per-library merge stay in one place.

include { hasContigGroups; readContigGroups; selectMultiqcGroup } from '../lib/contig_groups.nf'
include { gc_bias }                           from '../modules/gc_bias.nf'
include { gc_bias_by_contig_group }           from '../modules/gc_bias_by_contig_group.nf'
include { concatenate_gc_bias }               from '../modules/concatenate_gc_bias.nf'

workflow gc_bias_curves {
    take:
        md_bams  // tuple(library, bam, bai)

    main:
        reference_list = params.genomes[params.genome]

        gc_bias(md_bams.map { library, bam, bai ->
            tuple(library, bam, bai, reference_list.genome_fa)
        })
        gc_curves = gc_bias.out.for_agg
        // Falls through to the whole-reference curve for a genome with no contig groups; see the
        // warning below for why that is not the curve we would rather show.
        multiqc_curve = gc_bias.out.for_agg

        // Each organism needs a matching pair -- its own reads and its own reference -- because
        // picard normalizes against whatever reference it is handed.
        if (hasContigGroups(reference_list)) {
            contig_groups = readContigGroups(params.genome, reference_list)

            gc_bias_by_contig_group(
                md_bams.flatMap { library, bam, bai ->
                    contig_groups.collect { group ->
                        tuple(library, bam, bai, group.name, group.bed, group.contig_list,
                              group.fasta.toString())
                    }
                }
            )

            // One group's curve, under the filename the whole-reference curve has always had, so
            // the report stays comparable to historic runs.
            multiqc_group = selectMultiqcGroup(contig_groups, params.multiqc_gc_group)
            multiqc_curve = gc_bias_by_contig_group.out.for_multiqc
                .filter { library, group, curve -> group == multiqc_group }
                .map    { library, group, curve -> tuple(library, curve) }

            concatenate_gc_bias(
                gc_bias.out.for_agg
                    // size: lets a library's merge start as soon as its own groups are done;
                    // without it groupTuple waits for every library's groups first.
                    .join(gc_bias_by_contig_group.out.for_agg
                              .map { library, group, curve -> tuple(library, curve) }
                              .groupTuple(size: contig_groups.size()))
                    .map { library, whole_reference, group_curves ->
                        tuple(library, [whole_reference] + group_curves.sort(false) { it.name })
                    }
            )
            gc_curves = concatenate_gc_bias.out.for_agg
        }
        else {
            if (params.multiqc_gc_group) {
                log.warn "--multiqc_gc_group '${params.multiqc_gc_group}' was given but genome " +
                         "'${params.genome}' sets no gc_groups_dir, so there are no per-organism " +
                         "curves to choose from."
            }
            log.warn "Genome '${params.genome}' sets no gc_groups_dir, so the GC bias curve covers " +
                     "the whole reference, spike-in controls included. That curve sits higher at " +
                     "the GC-rich end than the autosome-only curve EM-seq reported historically, " +
                     "so do not compare it directly with older runs."
        }

    emit:
        for_agg = gc_curves
        // A single curve, because MultiQC's picard/gcbias module is configured with
        // use_filename_as_sample_name: every metrics block in a merged file resolves to one sample
        // name and the last silently wins, so a combined file would plot one arbitrary organism's
        // curve labelled as the library.
        for_multiqc = multiqc_curve
}
