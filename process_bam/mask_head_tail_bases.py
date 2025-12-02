#!/usr/bin/env python3

import pysam
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Set specified ends' base qualities to 0 for paired-end BAM.")
    parser.add_argument("--input_bam", required=True, help="Input BAM file (must be paired-end and sorted)")
    parser.add_argument("--output_bam", required=True, help="Output BAM file with modified qualities")
    parser.add_argument("--ignore_5prime_r1", type=int, default=None, help="Set R1 5' base qualities to 0 for first N bases")
    parser.add_argument("--ignore_3prime_r1", type=int, default=None, help="Set R1 3' base qualities to 0 for last N bases")
    parser.add_argument("--ignore_5prime_r2", type=int, default=None, help="Set R2 5' base qualities to 0 for first N bases")
    parser.add_argument("--ignore_3prime_r2", type=int, default=None, help="Set R2 3' base qualities to 0 for last N bases")
    return parser.parse_args()

def set_base_qualities_zero(qualities, head_trim=None, tail_trim=None):
    quals = list(qualities)
    length = len(quals)
    if head_trim and head_trim > 0:
        quals[:min(length, head_trim)] = [0] * min(length, head_trim)
    if tail_trim and tail_trim > 0:
        quals[-tail_trim:] = [0] * min(length, tail_trim)
    return quals

def main():
    args = parse_args()
    with pysam.AlignmentFile(args.input_bam, "rb") as infile, \
         pysam.AlignmentFile(args.output_bam, "wb", template=infile) as outfile:
        for read in infile:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                outfile.write(read)
                continue
            quals = list(read.query_qualities) if read.query_qualities is not None else None
            if quals is not None:
                if read.is_read1:
                    quals = set_base_qualities_zero(
                        quals,
                        args.ignore_5prime_r1,
                        args.ignore_3prime_r1
                    )
                elif read.is_read2:
                    quals = set_base_qualities_zero(
                        quals,
                        args.ignore_5prime_r2,
                        args.ignore_3prime_r2
                    )
                read.query_qualities = quals
            outfile.write(read)

if __name__ == "__main__":
    main()
