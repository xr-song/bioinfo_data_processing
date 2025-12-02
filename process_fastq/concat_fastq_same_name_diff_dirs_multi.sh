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
echo "$all_samples" | parallel -j 6 '
    sample={}
    echo "Processing: $sample"
    
    # collect all R1 files for this sample
    r1_files=""
    r2_files=""
    i1_files=""
    i2_files=""
    
    for dir in '"${input_dirs[*]}"'; do
        r1_found=$(find "$dir" -name "${sample}*R1*.fastq.gz" 2>/dev/null)
        r2_found=$(find "$dir" -name "${sample}*R2*.fastq.gz" 2>/dev/null)
	i1_found=$(find "$dir" -name "${sample}*I1*.fastq.gz" 2>/dev/null)
	i2_found=$(find "$dir" -name "${sample}*I2*.fastq.gz" 2>/dev/null)
        
        [ -n "$r1_found" ] && r1_files="$r1_files $r1_found"
        [ -n "$r2_found" ] && r2_files="$r2_files $r2_found"
	[ -n "$i1_found" ] && i1_files="$i1_files $i1_found"
	[ -n "$i2_found" ] && i2_files="$i2_files $i2_found"
    done
    
    # concatenate files
    [ -n "$r1_files" ] && cat $r1_files > "${output_dir}/${sample}_R1_001.fastq.gz"
    [ -n "$r2_files" ] && cat $r2_files > "${output_dir}/${sample}_R2_001.fastq.gz"
    [ -n "$i1_files" ] && cat $i1_files > "${output_dir}/${sample}_I1_001.fastq.gz"
    [ -n "$i2_files" ] && cat $i2_files > "${output_dir}/${sample}_I2_001.fastq.gz"

'

echo "Done! Check results in: $output_dir"
