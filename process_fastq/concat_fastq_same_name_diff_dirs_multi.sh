#!/bin/bash
# Concatenating FASTQ files from multiple directories
# Usage: ./scripts.sh dir1 dir2 [dir3 ...]

if [ $# -lt 2 ]; then
    echo "Usage: $0 dir1 dir2 [dir3 ...]"
    exit 1
fi

input_dirs=("$@")
output_dir="$PWD/rawdata_combined"
mkdir -p "$output_dir"

echo "Concatenating FASTQ files from ${#input_dirs[@]} directories"
echo "Output: $output_dir"

# get all unique sample IDs
all_samples=$(
    for dir in "${input_dirs[@]}"; do
        find "$dir" -name '*R1*.fastq.gz' | sed 's/.*\///g; s/_S[0-9]*.*//g; s/_L[0-9]*.*//g; s/_R1.*//g'
    done | sort -u
)

export output_dir
export -A input_dirs_array
for i in "${!input_dirs[@]}"; do
    input_dirs_array[$i]="${input_dirs[$i]}"
done

# process each sample
export log_dir="${output_dir}/logs"
mkdir -p "$log_dir"

echo "$all_samples" | parallel -j 6 '
    sample={}
    echo "Processing: $sample"
    log_file=${log_dir}/${sample}.log

    # Collect all files for each read type (and sort by filename for consistency)
    r1_files=$(for dir in '"${input_dirs[*]}"'; do find "$dir" -name "${sample}*R1*.fastq.gz" 2>/dev/null; done | sort)
    r2_files=$(for dir in '"${input_dirs[*]}"'; do find "$dir" -name "${sample}*R2*.fastq.gz" 2>/dev/null; done | sort)
    i1_files=$(for dir in '"${input_dirs[*]}"'; do find "$dir" -name "${sample}*I1*.fastq.gz" 2>/dev/null; done | sort)
    i2_files=$(for dir in '"${input_dirs[*]}"'; do find "$dir" -name "${sample}*I2*.fastq.gz" 2>/dev/null; done | sort)

    # Merge, log file order
    if [ -n "$r1_files" ]; then
        echo "Merging R1 to ${output_dir}/${sample}_R1_001.fastq.gz" > "$log_file"
        echo "R1 files:" >> "$log_file"
        echo "$r1_files" >> "$log_file"
        cat $r1_files > "${output_dir}/${sample}_R1_001.fastq.gz"
    fi

    if [ -n "$r2_files" ]; then
        echo "Merging R2 to ${output_dir}/${sample}_R2_001.fastq.gz" >> "$log_file"
        echo "R2 files:" >> "$log_file"
        echo "$r2_files" >> "$log_file"
        cat $r2_files > "${output_dir}/${sample}_R2_001.fastq.gz"
    fi

    if [ -n "$i1_files" ]; then
        echo "Merging I1 to ${output_dir}/${sample}_I1_001.fastq.gz" >> "$log_file"
        echo "I1 files:" >> "$log_file"
        echo "$i1_files" >> "$log_file"
        cat $i1_files > "${output_dir}/${sample}_I1_001.fastq.gz"
    fi

    if [ -n "$i2_files" ]; then
        echo "Merging I2 to ${output_dir}/${sample}_I2_001.fastq.gz" >> "$log_file"
        echo "I2 files:" >> "$log_file"
        echo "$i2_files" >> "$log_file"
        cat $i2_files > "${output_dir}/${sample}_I2_001.fastq.gz"
    fi
'
echo "Done! Check results in: $output_dir"
