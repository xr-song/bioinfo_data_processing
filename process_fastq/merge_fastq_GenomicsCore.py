#!/usr/bin/env python3

# Merge fastq file from GenomicsCore on i7. For easier file processing for sciMET

import os
import glob
import subprocess
from collections import defaultdict

def merge_fastq_files(input_dir, output_dir):
    """
    Merge FASTQ files by the part before '-' and assign sequential S numbers.
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get all fastq files
    fastq_files = glob.glob(os.path.join(input_dir, "*.fastq.gz"))
    
    if not fastq_files:
        print("No FASTQ files found in the directory")
        return
    
    # Group files by prefix (part before "-") and read type (I1, I2, R1, R2)
    file_groups = defaultdict(lambda: defaultdict(list))
    
    for file_path in fastq_files:
        filename = os.path.basename(file_path)
        
        # Extract prefix (part before "-")
        if "-" not in filename:
            print(f"Warning: No '-' found in {filename}, skipping...")
            continue
            
        prefix = filename.split("-")[0]
        
        # Determine read type (I1, I2, R1, R2)
        if "_I1_" in filename:
            read_type = "I1"
        elif "_I2_" in filename:
            read_type = "I2"
        elif "_R1_" in filename:
            read_type = "R1"
        elif "_R2_" in filename:
            read_type = "R2"
        else:
            print(f"Warning: Cannot determine read type for {filename}, skipping...")
            continue
        
        file_groups[prefix][read_type].append(file_path)
    
    # Sort prefixes to ensure consistent S number assignment
    sorted_prefixes = sorted(file_groups.keys())
    
    # Process each prefix
    for s_num, prefix in enumerate(sorted_prefixes, 1):
        s_string = f"S{s_num:02d}"  # S01, S02, S03, etc.
        
        print(f"Processing {prefix} -> {s_string}")
        
        # Process each read type for this prefix
        for read_type in ["I1", "I2", "R1", "R2"]:
            if read_type in file_groups[prefix]:
                input_files = sorted(file_groups[prefix][read_type])
                output_filename = f"{prefix}_{s_string}_{read_type}_001.fastq.gz"
                output_path = os.path.join(output_dir, output_filename)
                
                print(f"  Merging {len(input_files)} {read_type} files -> {output_filename}")
                
                # Merge files using cat (for gzipped files)
                try:
                    with open(output_path, 'wb') as outfile:
                        for input_file in input_files:
                            print(f"    Adding: {os.path.basename(input_file)}")
                            with open(input_file, 'rb') as infile:
                                outfile.write(infile.read())
                    
                    print(f"  ✓ Created: {output_filename}")
                    
                except Exception as e:
                    print(f"  ✗ Error merging {read_type} files for {prefix}: {e}")

def main():
    import sys
    
    if len(sys.argv) != 3:
        print("Usage: python merge_fastq.py <input_directory> <output_directory>")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} is not a valid directory")
        sys.exit(1)
    
    print(f"Processing FASTQ files in: {input_dir}")
    print(f"Output directory: {output_dir}")
    merge_fastq_files(input_dir, output_dir)
    print("Done!")

if __name__ == "__main__":
    main()
