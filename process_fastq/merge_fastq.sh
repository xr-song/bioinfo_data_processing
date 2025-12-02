#!/bin/bash

# merge all R1 R2 without considering sample names

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "usage: $0 <input_dir> <output_dir> <prefix>"
    exit 1
fi

input_dir="$1"
output_dir="$2"
prefix="$3"

mkdir -p "$output_dir"

r1_out="$output_dir/${prefix}.merged.R1.fastq.gz"
r2_out="$output_dir/${prefix}.merged.R2.fastq.gz"

echo "searching for R1 files..."
find "$input_dir" -type f -name "*R1*.fastq.gz" | sort | xargs cat > "$r1_out"

echo "searching for R2 files..."
find "$input_dir" -type f -name "*R2*.fastq.gz" | sort | xargs cat > "$r2_out"

echo "done:"
ls -lh "$r1_out" "$r2_out"

