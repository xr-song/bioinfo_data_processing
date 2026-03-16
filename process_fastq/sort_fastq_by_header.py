#!/usr/bin/env python3
import os
import glob
import gzip
import shutil
import tempfile
import argparse
import subprocess

def get_file_type(filename):
    if filename.endswith('.gz'):
        return 'gz'
    else:
        return 'plain'

def sort_fastq(input_file, output_file):
    file_type = get_file_type(input_file)
    with tempfile.NamedTemporaryFile(delete=False, mode='w') as temp_tab:
        # Convert FASTQ to tab-delimited records: header<TAB>seq<TAB>plus<TAB>qual
        if file_type == 'gz':
            opener = gzip.open
            mode = 'rt'
        else:
            opener = open
            mode = 'r'
        with opener(input_file, mode) as fin:
            while True:
                header = fin.readline()
                if not header:
                    break
                seq = fin.readline()
                plus = fin.readline()
                qual = fin.readline()
                temp_tab.write(f"{header.strip()}\t{seq.strip()}\t{plus.strip()}\t{qual.strip()}\n")

    # Sort by header (first column)
    temp_sorted = tempfile.NamedTemporaryFile(delete=False, mode='w')
    # Use system 'sort' for efficiency
    subprocess.run(f"sort -k1,1 {temp_tab.name} > {temp_sorted.name}", shell=True, check=True)
    os.remove(temp_tab.name)

    # Convert back tab-delimited records to FASTQ, compress if needed
    if output_file.endswith('.gz'):
        fout = gzip.open(output_file, 'wt')
    else:
        fout = open(output_file, 'w')
    with open(temp_sorted.name, 'r') as fin:
        for line in fin:
            header, seq, plus, qual = line.rstrip('\n').split('\t')
            fout.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
    fout.close()
    os.remove(temp_sorted.name)

def main(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    fq_exts = ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz']
    fq_files = []
    for ext in fq_exts:
        fq_files.extend(glob.glob(os.path.join(input_dir, ext)))
    if not fq_files:
        print("No FASTQ files found.", file=sys.stderr)
        return

    for input_file in fq_files:
        base = os.path.basename(input_file)
        output_file = os.path.join(output_dir, base)
        print(f"Sorting {base} ...")
        sort_fastq(input_file, output_file)
        print(f"✓ Wrote {output_file}")

if __name__ == "__main__":
    import sys
    parser = argparse.ArgumentParser(description="Sort all FASTQ files in a directory by read name.")
    parser.add_argument("-i", "--input_dir", required=True, help="Input directory")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory for sorted FASTQ files")
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
