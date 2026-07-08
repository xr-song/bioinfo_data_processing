#!/usr/bin/env python3

import gzip
import sys
import os
from collections import defaultdict
import argparse

try:
    import pybedtools
    PYBEDTOOLS_AVAILABLE = True
except ImportError:
    print("Warning: pybedtools not available.")


def load_cpg_bed_pybedtools(bed_file):
    cpg_positions = defaultdict(set)
    
    bed = pybedtools.BedTool(bed_file)
    
    for interval in bed:
        chrom = str(interval.chrom)
        start = int(interval.start)  # 0-based
        pos1=start
        pos2=start+1
        cpg_positions[chrom].add((pos1,pos2))
    
    total = sum(len(v) for v in cpg_positions.values())
    print(f"Loaded {total} CpG positions from {bed_file}")
    return cpg_positions

def merge_allc_to_cpg_bed(allc_file, cpg_positions, output_file):
    merged_data = defaultdict(lambda: {'mc': 0, 'cov': 0})
    
    pos_lookup = {}
    pos_lookup2 = {}
    for chrom, positions in cpg_positions.items():
        pos_lookup[chrom] = {pos[0]: pos[1] for pos in positions}
        pos_lookup2[chrom] = {pos[1]: pos[0] for pos in positions}
    
    if allc_file.endswith('.gz'):
        file_handle = gzip.open(allc_file, 'rt')
    else:
        file_handle = open(allc_file, 'r')
    
    try:
        for line in file_handle:
            line = line.strip()
            if not line:
                continue
            
            columns = line.split('\t')
            if len(columns) < 7:
                continue
            
            try:
                chrom = str(columns[0])
                pos_1based = int(columns[1])  # 1-based position from ALLC
                mc = int(columns[4])
                cov = int(columns[5])
            except (ValueError, IndexError):
                continue
            
            # Fast lookup using dictionary
            if chrom in pos_lookup:
                if pos_1based in pos_lookup[chrom]:
                    pos1 = pos_1based
                elif pos_1based in pos_lookup2[chrom]:
                    pos1 = pos_lookup2[chrom][pos_1based]
                else:
                    continue
                key = (chrom, pos1)
                merged_data[key]['mc'] += mc
                merged_data[key]['cov'] += cov
    
    finally:
        file_handle.close()
    
    output_handle = gzip.open(output_file, 'wt') if output_file.endswith('.gz') else open(output_file, 'w')
    
    try:
        for (chrom, pos1), values in sorted(merged_data.items()):
            mc = values['mc']
            cov = values['cov']
            # BED format: chrom, chromStart (0-based), chromEnd, mc, cov
            output_line = f"{chrom}\t{pos1-1}\t{pos1 + 1}\t{mc}\t{cov}\n"
            output_handle.write(output_line)
    
    finally:
        output_handle.close()
    
    print(f"Wrote {len(merged_data)} merged CpG sites to {output_file}")
    return len(merged_data)


def main(allc_file, bed_file, output_file, use_pybedtools=True):
    
    if use_pybedtools and PYBEDTOOLS_AVAILABLE:
        print("Using pybedtools for faster BED loading...")
        cpg_positions = load_cpg_bed_pybedtools(bed_file)
    
    num_sites = merge_allc_to_cpg_bed(allc_file, cpg_positions, output_file)
    
    print(f"Merged {num_sites} CpG dinucleotide sites")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate merged (CpG dinucleotide) level BED file from input ALLC using reference BED"
    )
    parser.add_argument(
        '--allc',
        required=True,
        help="Input ALLC file (gzipped or uncompressed, 1-based positions)"
    )
    parser.add_argument(
        '--bed',
        required=True,
        help="Reference BED file with CpG positions (1-based, standard BED format: chrom\\tstart\\tend)"
    )
    parser.add_argument(
        '--out',
        required=True,
        help="Output merged BED file (0-based, format: chrom\\tstart\\tend\\tmc\\tcov)"
    )
    parser.add_argument(
        '--no-pybedtools',
        action='store_true',
        help="Force manual BED loading instead of using pybedtools"
    )
    
    args = parser.parse_args()
    
    allc_file = args.allc
    bed_file = args.bed
    output_file = args.out
    use_pybedtools = not args.no_pybedtools
    
    if not os.path.isfile(allc_file):
        print(f"Error: {allc_file} not found")
        sys.exit(1)
    
    if not os.path.isfile(bed_file):
        print(f"Error: {bed_file} not found")
        sys.exit(1)
    
    main(allc_file, bed_file, output_file, use_pybedtools)
