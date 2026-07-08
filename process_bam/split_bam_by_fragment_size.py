#!/usr/bin/env python3
"""
Split a BAM file into 4 BAM files based on fragment size:
  - small.bam:        <= 50 bp
  - medium.bam:       51-100 bp
  - long.bam:         101-150 bp
  - xlarge.bam:       > 150 bp
"""

import argparse
import pysam
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-bam", "-i", required=True, help="Input BAM file")
    ap.add_argument("--output-dir", "-o", default=".", help="Output directory (default: current dir)")
    args = ap.parse_args()

    input_bam = args.input_bam
    output_dir = args.output_dir

    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)

    # Generate output filenames based on input BAM name
    base_name = os.path.splitext(os.path.basename(input_bam))[0]
    out_small = os.path.join(output_dir, f"{base_name}.small_le50bp.bam")
    out_medium = os.path.join(output_dir, f"{base_name}.medium_51-100bp.bam")
    out_long = os.path.join(output_dir, f"{base_name}.long_101-150bp.bam")
    out_xlarge = os.path.join(output_dir, f"{base_name}.xlarge_gt150bp.bam")

    inbam = pysam.AlignmentFile(input_bam, "rb")

    # Open output BAM files with same header as input
    out_small_fh = pysam.AlignmentFile(out_small, "wb", template=inbam)
    out_medium_fh = pysam.AlignmentFile(out_medium, "wb", template=inbam)
    out_long_fh = pysam.AlignmentFile(out_long, "wb", template=inbam)
    out_xlarge_fh = pysam.AlignmentFile(out_xlarge, "wb", template=inbam)

    read_buffer = {}
    frag_counts = {"small": 0, "medium": 0, "long": 0, "xlarge": 0}

    for read in inbam.fetch(until_eof=True):
        # Skip unpaired or unmapped reads
        if not read.is_paired or read.is_unmapped or read.mate_is_unmapped:
            continue

        if read.query_name not in read_buffer:
            read_buffer[read.query_name] = read
        else:
            mate = read_buffer.pop(read.query_name)

            # Ensure read1 and read2 are in the correct order
            if read.is_read1:
                read1, read2 = read, mate
            else:
                read1, read2 = mate, read

            # Skip if on different chromosomes
            if read1.reference_id != read2.reference_id:
                continue

            # Calculate fragment size and orientation
            frag_start = min(read1.reference_start, read2.reference_start)
            frag_end = max(read1.reference_end, read2.reference_end)
            frag_size = frag_end - frag_start

            # Categorize by size and write to appropriate file
            if frag_size <= 50:
                out_small_fh.write(read1)
                out_small_fh.write(read2)
                frag_counts["small"] += 1
            elif 51 <= frag_size <= 100:
                out_medium_fh.write(read1)
                out_medium_fh.write(read2)
                frag_counts["medium"] += 1
            elif 101 <= frag_size <= 150:
                out_long_fh.write(read1)
                out_long_fh.write(read2)
                frag_counts["long"] += 1
            else:  # > 150
                out_xlarge_fh.write(read1)
                out_xlarge_fh.write(read2)
                frag_counts["xlarge"] += 1

    # Close all files
    inbam.close()
    out_small_fh.close()
    out_medium_fh.close()
    out_long_fh.close()
    out_xlarge_fh.close()

    # Index output files
    #pysam.index(out_small)
    #pysam.index(out_medium)
    #pysam.index(out_long)
    #pysam.index(out_xlarge)

    # Print summary
    total = sum(frag_counts.values())
    print(f"Fragment size distribution:")
    print(f"  ≤50 bp (small):       {frag_counts['small']:,} ({100*frag_counts['small']/total:.1f}%)")
    print(f"  51-100 bp (medium):   {frag_counts['medium']:,} ({100*frag_counts['medium']/total:.1f}%)")
    print(f"  101-150 bp (long):    {frag_counts['long']:,} ({100*frag_counts['long']/total:.1f}%)")
    print(f"  >150 bp (xlarge):     {frag_counts['xlarge']:,} ({100*frag_counts['xlarge']/total:.1f}%)")
    print(f"  Total:                {total:,}")
    print()
    print(f"Output BAM files:")
    print(f"  {out_small}")
    print(f"  {out_medium}")
    print(f"  {out_long}")
    print(f"  {out_xlarge}")


if __name__ == "__main__":
    main()
