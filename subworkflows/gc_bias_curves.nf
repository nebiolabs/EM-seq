// GC bias curves for one set of aligned BAMs: always the whole-reference curve, plus one curve per
// organism when the genome is a composite that ships a gc_groups_dir.
//
// Encapsulated as a subworkflow because a process may only be invoked once per workflow context:
// the whole-reference and per-group curves are the same process, so the second invocation must be
// an alias.

include { hasContigGroups; readContigGroups } from '../lib/contig_groups.nf'
include { gc_bias; gc_bias as gc_bias_by_group } from '../modules/gc_bias.nf'
include { subset_bam_to_contig_group }          from '../modules/subset_bam_to_contig_group.nf'
include { concatenate_gc_bias }                 from '../modules/concatenate_gc_bias.nf'

workflow gc_bias_curves {
    take:
        md_bams  // tuple(library, bam, bai)

    main:
        reference_list = params.genomes[params.genome]

        gc_bias(md_bams.map { library, bam, bai ->
            tuple(library, bam, bai, '', reference_list.genome_fa)
        })
        gc_curves = gc_bias.out.for_agg

        // Each organism needs a matching pair because picard normalizes against
        // whatever reference it is handed.
        if (hasContigGroups(reference_list)) {
            contig_groups = readContigGroups(params.genome, reference_list)
            group_by_name = contig_groups.collectEntries { [(it.name): it] }

            subset_bam_to_contig_group(
                md_bams.flatMap { library, bam, bai ->
                    contig_groups.collect { group ->
                        tuple(library, bam, bai, group.name, group.bed, group.contig_list)
                    }
                }
            )

            gc_bias_by_group(
                subset_bam_to_contig_group.out.bam.map { library, group, bam, bai ->
                    tuple(library, bam, bai, group, group_by_name[group].fasta.toString())
                }
            )


            concatenate_gc_bias(
                gc_bias.out.for_agg
                    // size: lets a library's merge start as soon as its own groups are done;
                    // without it groupTuple waits for every library's groups first.
                    .join(gc_bias_by_group.out.for_agg.groupTuple(size: contig_groups.size()))
                    .map { library, whole_reference, group_curves ->
                        tuple(library, [whole_reference] + group_curves.sort(false) { it.name })
                    }
            )
            gc_curves = concatenate_gc_bias.out.for_agg
        }

    emit:
        for_agg = gc_curves
}
