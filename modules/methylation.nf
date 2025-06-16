process methylDackel_mbias {
    label 'medium_cpu'
    errorStrategy 'retry'
    tag "${library}"
    conda "bioconda::methyldackel=0.6.1 bioconda::samtools=1.21 conda-forge::pigz=2.8 conda-forge::sed=4.9"
    publishDir "${params.outputDir}/methylDackelExtracts/mbias"

    input:
        tuple val(library), path(md_bam), path(md_bai), val(barcodes)
        tuple path(genome_fa), path(genome_fai)

    output:
        path('*.svg'), emit: mbias_output_svg
        path('*.tsv'), emit: mbias_output_tsv
        tuple val(params.email), val(library), path('*.tsv'), emit: for_agg

    script:
    """
    echo -e "chr\tcontext\tstrand\tRead\tPosition\tnMethylated\tnUnmethylated\tnMethylated(+dups)\tnUnmethylated(+dups)" > ${library}_${barcodes}_combined_mbias.tsv
    chrs=(`samtools view -H "${md_bam}" | grep @SQ | cut -f 2 | sed 's/SN://'| grep -v _random | grep -v chrUn | sed 's/|/\\|/'`)

    for chr in \${chrs[*]}; do
        for context in CHH CHG CpG; do
            arg=''
            if [ "\$context" = 'CHH' ]; then
            arg='--CHH --noCpG'
            elif [ "\$context" = 'CHG' ]; then
            arg='--CHG --noCpG'
            fi
            # need two calls to add columns containing the counts without filtering duplicate reads (for rrEM-seq where start/end is constrained)
            # not sure why we need both --keepDupes and -F, probably a bug in mbias
            join -t \$'\t' -j1 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -a 1 -e 0 \
            <( \
                MethylDackel mbias --noSVG \$arg -@ ${task.cpus} -r \$chr ${genome_fa} "${md_bam}" | \
                tail -n +2 | awk '{print \$1"-"\$2"-"\$3"\t"\$0}' | sort -k 1b,1
            ) \
            <( \
                MethylDackel mbias --noSVG --keepDupes -F 2816 \$arg -@ ${task.cpus} -r \$chr ${genome_fa} "${md_bam}" | \
                tail -n +2 | awk '{print \$1"-"\$2"-"\$3"\t"\$0}' | sort -k 1b,1
            ) \
            | sed "s/^/\${chr}\t\${context}\t/" \
            >> ${library}_${barcodes}_combined_mbias.tsv
        done
    done
    # makes the svg files for trimming checks
    MethylDackel mbias -@ ${task.cpus} --noCpG --CHH --CHG -r \${chrs[0]} ${genome_fa} "${md_bam}" ${library}_chn
    for f in *chn*.svg; do sed -i "s/Strand<\\/text>/Strand \$f \${chrs[0]} CHN <\\/text>/" \$f; done;

    MethylDackel mbias -@ ${task.cpus} -r \${chrs[0]} ${genome_fa} "${md_bam}" ${library}_cpg
    for f in *cpg*.svg; do sed -i "s/Strand<\\/text>/Strand \$f \${chrs[0]} CpG<\\/text>/" \$f; done;
    """
}


process methylDackel_extract {
    label 'high_cpu'
    tag "${library}"
    publishDir "${params.outputDir}/methylDackelExtracts", mode: 'copy'
    conda "bioconda::methyldackel=0.6.1 bioconda::samtools=1.21 conda-forge::pigz=2.8"

    input:
        tuple val(library), path(md_bam), path(md_bai), val(barcodes) 
        tuple path(genome_fa), path(genome_fai)

    output:
        tuple val(library), path('*CHG.methylKit.gz'), path('*CHH.methylKit.gz'),path('*CpG.methylKit.gz'), emit: extract_output 

    script:
    """
    MethylDackel extract --methylKit -q 20 --nOT 0,0,0,5 --nOB 0,0,5,0 -@ ${task.cpus} \
        --CHH --CHG -o ${library}.${barcodes} ${genome_fa} "${md_bam}" 
    pigz -p ${task.cpus} *.methylKit 
    """
}


process convert_methylkit_to_bed {
    label 'low_cpu'
    tag "${library}"
    conda "conda-forge::pigz=2.8 conda-forge::gawk=5.3.1 conda-forge::sed=4.9"

    input:
        tuple val(library), path(methylkit_CHH_gz),path(methylkit_CHG_gz),path(methylkit_CpG_gz),path(genome_fa),path(genome_fai)

    output:
        tuple val(library), path('*.bed'), emit: methylkit_bed

    script:
    """
    methylkit_basename=\$(basename "${methylkit_CpG_gz}" _CpG.methylKit.gz)
    
    # create sed scripts to interconvert chr names from genome index with line numbers for efficient sorting
    awk '{print "s/\\t"\$1"\\t/\\t"NR-1"\\t/"}' "${genome_fai}" > convert_chr_to_num.sed
    awk '{print "s/^"NR-1"\\t/"\$1"\\t/"}' "${genome_fai}" > revert_num_to_chr.sed 
    # merge initial methylkit files into a single BED format
    # Convert methylKit to BED format (0-based coordinates)
    # methylKit format: chr.base\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT
    # Output format: chr\tstart\tend\tname\tmethylation\tstrand
    sort -m -k2,2n -k3,3n \
            <(pigz -d -c "${methylkit_CHH_gz}" | sed 's/\$/\\tCHH/' | sed -f convert_chr_to_num.sed) \
            <(pigz -d -c "${methylkit_CHG_gz}" | sed 's/\$/\\tCHG/' | sed -f convert_chr_to_num.sed) \
            <(pigz -d -c "${methylkit_CpG_gz}" | sed 's/\$/\\tCpG/' | sed -f convert_chr_to_num.sed) \
    | awk -v OFS="\\t" 'NR>3{print \$2,\$3-1,\$3,\$8,\$6,(\$4 == "F") ? "+" : "-"}' \
    | sed -f revert_num_to_chr.sed \
    > "\${methylkit_basename}.bed"
    """
}

process prepare_target_bed {
    label 'low_cpu'
    tag "${target_bed.baseName}"
    conda "bioconda::bedtools=2.31.1"

    input:
        path(target_bed)
        tuple path(genome_fa), path(genome_fai)

    output:
        path('*_slop_sorted.bed'), emit: prepared_bed

    script:
    """
    # Apply slop to target BED file (50bp on each side)
    slop_len=50
    target_basename=\$(basename "${target_bed}" .bed)
    
    echo "Applying \${slop_len} bp slop to \${target_basename}..."
    bedtools slop -i <(grep -v '^#' "${target_bed}") -g "${genome_fai}" -b \${slop_len} \\
        | sort -k1,1 -k2,2n > \${target_basename}_slop_sorted.bed
    """
}

process intersect_beds {
    label 'low_cpu'
    tag "${library}"
    conda "bioconda::bedtools=2.31.1"

    input:
        tuple val(library), path(methylkit_bed)
        path(target_bed_prepared)
        tuple path(genome_fa), path(genome_fai)

    output:
        tuple val(library), path('*_intersect.tsv'), emit: intersections

    script:
    """
    methylkit_basename=\$(basename "${methylkit_bed}" .bed)
    target_basename=\$(basename "${target_bed_prepared}" _slop_sorted.bed)
    
    # Uses bedtools intersect with -wa -wb to keep both annotations
    # nonamecheck since control contigs don't conform
    bedtools intersect -nonamecheck \\
        -a "${target_bed_prepared}" \\
        -b "${methylkit_bed}" \\
        -wa -wb -sorted -g "${genome_fai}" > \${methylkit_basename}_vs_\${target_basename}_intersect.tsv
    """
}

process process_intersections {
    label 'low_cpu'
    tag "${library}"
    publishDir "${params.outputDir}/methylKit_intersections", mode: 'copy'
    conda "conda-forge::gawk=5.3.1"
    input:
        tuple val(library), path(intersect_file)

    output:
        path('*_intersections.tsv'), emit: intersection_results
        path('*_summary.tsv'), emit: intersection_summary
        tuple val(library), path('*_intersections.tsv'), path('*_summary.tsv'), emit: for_agg

    script:
    """
    # Extract base names from the intersect file
    intersect_basename=\$(basename "${intersect_file}" _intersect.tsv)
    
    output_file="\${intersect_basename}_intersections.tsv"
    summary_file="\${intersect_basename}_summary.tsv"
    
    # Add headers
    echo -e "methylkit_file\\tchr\\tstart\\tend\\tcontext\\tmethylation\\ttarget_locus\\ttarget_name" > \${output_file}
    echo -e "methylkit_file\\ttarget_length\\tposition\\tcontext\\tmean_methylation\\tn_loci\\tn_measurements" > \${summary_file}

    # Extract methylkit basename from the intersect filename
    methylkit_basename=\$(echo "\${intersect_basename}" | cut -d'_' -f1)

    # Process intersection results with awk
    if [ -s "${intersect_file}" ]; then
        awk -v methylkit_basename="\${methylkit_basename}" \\
            -v output_file="\${output_file}" \\
            -v summary_file="\${summary_file}" '
        BEGIN { 
            OFS="\\t"
        }
        !/^#/ {
            # Extract target info (first part from BED)
            target_chr = \$1
            target_start = \$2
            target_end = \$3
            target_name = \$4
            
            # Extract methylkit info (second part after intersection)
            # Assuming target BED has 4+ columns, methylkit data starts after that
            num_target_cols = 4
            if (NF > num_target_cols + 4) {
                # Try to detect actual number of target columns
                for (i = 5; i <= NF-4; i++) {
                    if (\$i ~ /^chr/ && \$(i+1) ~ /^[0-9]+\$/ && \$(i+2) ~ /^[0-9]+\$/) {
                        num_target_cols = i - 1
                        break
                    }
                }
            }
            
            meth_chr = \$(num_target_cols + 1)
            meth_start = \$(num_target_cols + 2)
            meth_end = \$(num_target_cols + 3)
            meth_context = \$(num_target_cols + 4)
            meth_value = \$(num_target_cols + 5)
            
            # Create locus string (convert 0-based BED to 1-based for display)
            locus = target_chr ":" (target_start + 1) "-" target_end
            
            # Calculate target length and position within target
            target_length = target_end - target_start
            # Position relative to start of target (1-based)
            position_in_target = meth_start - target_start + 1
            
            # Output detailed results
            print methylkit_basename, meth_chr, meth_start, meth_end, meth_context, meth_value, locus, target_name >> output_file
            
            # Collect summary data
            key = methylkit_basename "\\t" target_length "\\t" position_in_target "\\t" meth_context
            sum[key] += meth_value
            count[key]++
            loci[key][locus] = 1  # Track unique loci
        }
        END {
            # Write summary data
            for (k in sum) {
                split(k, parts, "\\t")
                methylkit_file = parts[1]
                target_length = parts[2] 
                position = parts[3]
                context = parts[4]
                mean_meth = sum[k] / count[k]
                n_measurements = count[k]
                n_loci = length(loci[k])
                print methylkit_file, target_length, position, context, mean_meth, n_loci, n_measurements >> summary_file
            }
        }' "${intersect_file}"
    else
        echo "No intersections found in ${intersect_file}"
        # Create empty files with headers only
        touch \${output_file}
        touch \${summary_file}
    fi
    """
}

process concatenate_intersections {
    label 'low_cpu'
    publishDir "${params.outputDir}/methylKit_intersections", mode: 'copy'

    input:
        path(intersection_files)
        path(summary_files)

    output:
        path('all_intersections_combined.tsv'), emit: combined_intersections
        path('all_summaries_combined.tsv'), emit: combined_summaries

    when:
        params.target_bed != 'undefined'

    script:
    """
    # Combine all intersection files
    echo -e "methylkit_file\\tchr\\tstart\\tend\\tcontext\\tmethylation\\ttarget_locus\\ttarget_name" > all_intersections_combined.tsv
    
    # Add all intersection data (skip headers from individual files)
    for file in ${intersection_files}; do
        if [ -s "\$file" ]; then
            tail -n +2 "\$file" >> all_intersections_combined.tsv
        fi
    done
    
    # Combine all summary files
    echo -e "methylkit_file\\ttarget_length\\tposition\\tcontext\\tmean_methylation\\tn_loci\\tn_measurements" > all_summaries_combined.tsv
    
    # Add all summary data (skip headers from individual files)
    for file in ${summary_files}; do
        if [ -s "\$file" ]; then
            tail -n +2 "\$file" >> all_summaries_combined.tsv
        fi
    done
    
    # Report statistics
    total_intersections=\$(tail -n +2 all_intersections_combined.tsv | wc -l)
    total_summaries=\$(tail -n +2 all_summaries_combined.tsv | wc -l)
    
    echo "Combined \${total_intersections} intersection records from \$(echo ${intersection_files} | wc -w) files"
    echo "Combined \${total_summaries} summary records from \$(echo ${summary_files} | wc -w) files"
    """
}
