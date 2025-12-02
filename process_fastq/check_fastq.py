#!/usr/bin/env python3
import os
import sys
import gzip
import argparse
import multiprocessing as mp
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

def check_fastq_file(filepath):
    """Check FASTQ file integrity and count reads."""
    try:
        read_count = 0
        line_count = 0
        
        with gzip.open(filepath, 'rt') as f:
            for line in f:
                line_count += 1
                line_type = line_count % 4
                
                if line_type == 1:  # Header line
                    if not line.startswith('@'):
                        return filepath.name, "INVALID_HEADER", 0
                    read_count += 1
                elif line_type == 3:  # Quality header
                    if not line.startswith('+'):
                        return filepath.name, "INVALID_QUALITY_HEADER", 0
        
        # Check if complete FASTQ records
        if line_count % 4 != 0:
            return filepath.name, "INCOMPLETE_RECORD", 0
            
        return filepath.name, "PASS", read_count
        
    except Exception as e:
        return filepath.name, f"ERROR: {str(e)}", 0

def main():
    parser = argparse.ArgumentParser(description='FASTQ integrity and read count check')
    parser.add_argument('-i', '--input', required=True, help='Input directory')
    parser.add_argument('-p', '--cores', type=int, default=4, help='Number of cores')
    
    args = parser.parse_args()
    
    input_dir = Path(args.input)
    if not input_dir.exists():
        print(f"Error: {input_dir} does not exist")
        sys.exit(1)
    
    print(f"=== FASTQ Integrity and Read Count Check ===")
    print(f"Input: {input_dir}")
    print(f"Cores: {args.cores}")
    print()
    
    # Find FASTQ files
    fastq_files = list(input_dir.glob("*.fastq.gz"))
    if not fastq_files:
        print("No FASTQ.gz files found")
        sys.exit(1)
    
    print(f"Found {len(fastq_files)} FASTQ files")
    print()
    
    # Process files in parallel
    print("=== Processing Files ===")
    with ProcessPoolExecutor(max_workers=args.cores) as executor:
        results = list(executor.map(check_fastq_file, fastq_files))
    
    print()
    
    # Analyze results
    passed = []
    failed = []
    
    for filename, status, read_count in results:
        if status == "PASS":
            passed.append((filename, read_count))
            print(f"✅ {filename}: {read_count:,} reads")
        else:
            failed.append((filename, status))
            print(f"❌ {filename}: {status}")
    
    print()
    print("=== Summary ===")
    print(f"Total files: {len(fastq_files)}")
    print(f"Passed: {len(passed)}")
    print(f"Failed: {len(failed)}")
    
    if passed:
        total_reads = sum(count for _, count in passed)
        print(f"Total reads: {total_reads:,}")
        print(f"Average reads per file: {total_reads//len(passed):,}")
    
    if failed:
        print(f"\nFailed files:")
        for filename, status in failed:
            print(f"  {filename}: {status}")
        sys.exit(1)
    else:
        print("\n✅ All files passed!")

if __name__ == "__main__":
    main()
