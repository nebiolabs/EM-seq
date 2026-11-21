#!/usr/bin/env python
"""Count CpG dinucleotides per contig in a reference FASTA.

The count is the denominator for "fraction of reference Cs covered at depth X".
Each CG dinucleotide holds two CpG-context Cs, one on each strand, so n_c = 2 * n_cpg.
MethylDackel reports those two Cs separately, and the depth-dependent analyses count them
separately.

Contigs are chosen in one of two ways. With --contig-groups, keep contigs whose
contig_type is in --include-types. Without it, keep every contig whose name does
not match --exclude-regex. Excluding spike-in controls and organelles keeps their
extreme depth and fixed methylation out of genome-wide correlations.

Output columns: contig, length, n_cpg, n_c.
"""
import argparse
import re
import sys

DEFAULT_INCLUDE_TYPES = 'autosomal,sex_chromosome'
# Used only without a contig-groups file. Matches the spike-in controls and
# organelles in the bundled references.
DEFAULT_EXCLUDE_REGEX = (
    r'(?i)^(chrM|MT|M|chrEBV|EBV|phage_.*|plasmid_.*|lambda.*|puc19.*|pUC19.*|'
    r'T4|Xp12|chrC|Pt|chloroplast|mitochondri.*)$'
)


def read_contig_groups(path, include_types):
    keep = set()
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        try:
            name_i = header.index('contig_name')
            type_i = header.index('contig_type')
        except ValueError:
            sys.exit(f'{path}: needs contig_name and contig_type columns, got {header}')
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) > type_i and fields[type_i] in include_types:
                keep.add(fields[name_i])
    return keep


def iter_fasta(path):
    """Yield (name, sequence bytes) for each record."""
    name, chunks = None, []
    with open(path, 'rb') as fh:
        for line in fh:
            if line.startswith(b'>'):
                if name is not None:
                    yield name, b''.join(chunks)
                name = line[1:].split()[0].decode()
                chunks = []
            else:
                chunks.append(line.rstrip())
    if name is not None:
        yield name, b''.join(chunks)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fasta')
    ap.add_argument('-o', '--out', default='reference_cpgs.tsv')
    ap.add_argument('--contig-groups', help='contig_groups.tsv with a contig_type column')
    ap.add_argument('--include-types', default=DEFAULT_INCLUDE_TYPES,
                    help=f'comma-separated contig_type values to keep (default {DEFAULT_INCLUDE_TYPES})')
    ap.add_argument('--exclude-regex', default=DEFAULT_EXCLUDE_REGEX,
                    help='contig names to drop when no --contig-groups is given')
    args = ap.parse_args()

    if args.contig_groups:
        keep = read_contig_groups(args.contig_groups, set(args.include_types.split(',')))
        if not keep:
            sys.exit(f'no contigs in {args.contig_groups} have contig_type in {args.include_types}')

        def wanted(name):
            return name in keep
    else:
        exclude = re.compile(args.exclude_regex)

        def wanted(name):
            return not exclude.match(name)

    n_kept = 0
    with open(args.out, 'w') as out:
        out.write('contig\tlength\tn_cpg\tn_c\n')
        for name, seq in iter_fasta(args.fasta):
            if not wanted(name):
                continue
            n_cpg = seq.upper().count(b'CG')
            out.write(f'{name}\t{len(seq)}\t{n_cpg}\t{2 * n_cpg}\n')
            n_kept += 1
    if n_kept == 0:
        sys.exit('no contigs kept; check --contig-groups / --exclude-regex')
    print(f'[count_reference_cpgs] kept {n_kept} contigs', file=sys.stderr)


if __name__ == '__main__':
    main()
