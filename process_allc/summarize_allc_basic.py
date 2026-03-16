#!/usr/bin/env python3
import gzip
import sys

def summarize_allc(input_file, output_file):
    total_cov = 0
    total_count = 0
    total_methylated = 0

    with gzip.open(input_file, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue  # skip malformed lines
            try:
                meth = int(fields[4])
                cov = int(fields[5])
            except ValueError:
                continue  # skip lines with non-integer meth/cov
            total_methylated += meth
            total_cov += cov
            total_count += 1  # each line counts as one site

    with open(output_file, 'w') as out:
        out.write("total_cov\ttotal_count\ttotal_methylated\n")
        out.write(f"{total_cov}\t{total_count}\t{total_methylated}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input.allc.tsv.gz output_summary.tsv")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    summarize_allc(input_file, output_file)

