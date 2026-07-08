#!/bin/bash

input_dir1=$1
input_dir2=$2
output_dir=$3

mkdir -p "$output_dir"

# Find all fastq.gz files in dir1 and concatenate with matching file in dir2
find "$input_dir1" -name "*.fastq.gz" -type f | \
    parallel -j 8 '
        filename=$(basename {})
        cat {} '"$input_dir2"'/"$filename" > '"$output_dir"'/"$filename"
        echo "Done: $filename"
    '

echo "All files concatenated to $output_dir"
