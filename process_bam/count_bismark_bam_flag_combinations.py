#!/usr/bin/env python3
import pysam
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(
    description="Count all combinations of read flags, XR, and XG tags in a BAM file"
)
parser.add_argument("--input_bam", required=True, help="Input BAM file")
parser.add_argument("--output_tsv", required=True, help="Output TSV summary file")
args = parser.parse_args()

# Dictionary to store counts: key = (is_read1, FLAG, XR, XG, READ), value = count
combo_counts = defaultdict(int)

bam_in = pysam.AlignmentFile(args.input_bam, "rb")

# Track read names to determine order within pairs
read_order = {}  # key = read_name, value = order (1 or 2)

for read in bam_in.fetch(until_eof=True):
    if read.is_unmapped:
        continue

    is_read1 = read.is_read1
    flag = read.flag
    read_name = read.query_name

    # Determine READ order (1 = first appearance, 2 = second appearance)
    if read_name not in read_order:
        read_order[read_name] = 1
        read_pos = 1
    else:
        read_pos = 2
        # Clean up to save memory (assuming pairs are close together)
        del read_order[read_name]

    # Get XR and XG tags, handle missing tags
    try:
        xr = read.get_tag("XR")
    except KeyError:
        xr = "NA"

    try:
        xg = read.get_tag("XG")
    except KeyError:
        xg = "NA"

    # Store the combination
    combo_counts[(is_read1, flag, xr, xg, read_pos)] += 1

bam_in.close()

# Write output TSV
with open(args.output_tsv, 'w') as out:
    out.write("is_read1\tFLAG\tXR\tXG\tREAD\tcount\n")

    # Sort by count descending, then by other fields
    for (is_read1, flag, xr, xg, read_pos), count in sorted(
        combo_counts.items(),
        key=lambda x: (-x[1], x[0][1], x[0][0], x[0][4])  # sort by count desc, then flag, then is_read1, then READ
    ):
        out.write(f"{is_read1}\t{flag}\t{xr}\t{xg}\t{read_pos}\t{count}\n")

print(f"Found {len(combo_counts)} unique combinations")
print(f"Total reads analyzed: {sum(combo_counts.values())}")
print(f"Output written to: {args.output_tsv}")
