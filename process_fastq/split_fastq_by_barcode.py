#!/usr/bin/env python3
"""
Split FASTQ files by barcode extracted from the header.
Barcodes are extracted from the last field of the header (e.g., AATTGAGA+GTAGCTCC)
"""

import gzip
import sys
import os
from collections import defaultdict
from pathlib import Path

def open_gz(filename, mode="rt"):
    """Open gzipped or regular file"""
    if str(filename).endswith(".gz"):
        return gzip.open(filename, mode=mode)
    return open(filename, mode=mode)

def extract_barcode(header):
    """
    Extract barcode from FASTQ header.
    Format: @LH00478:266:22YHLKLT3:5:1101:28840:1048 1:N:0:AATTGAGA+GTAGCTCC
    Barcode is the last field after the last colon: AATTGAGA+GTAGCTCC
    """
    # Get the last colon-delimited field
    parts = header.strip().split(':')
    if len(parts) > 0:
        barcode = parts[-1]
        return barcode
    return "unknown"

def fastq_parser(fastq_file):
    """Parse FASTQ file and yield (header, seq, plus, qual, barcode)"""
    with open_gz(fastq_file, 'rt') as f:
        while True:
            name = f.readline().strip()
            if not name:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            
            barcode = extract_barcode(name)
            yield name, seq, plus, qual, barcode

def split_fastq_by_barcode(fastq_file, output_dir):
    """
    Split FASTQ file by barcode.
    Creates separate .fastq.gz files for each barcode found.
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get input filename (without .gz)
    input_filename = Path(fastq_file).stem
    if input_filename.endswith('.fastq'):
        input_filename = input_filename  # Already has .fastq
    else:
        input_filename = input_filename + '.fastq'
    
    # Extract file type (R1, R2, I1, I2, etc) from filename
    parts = input_filename.split('_')
    file_type = None
    for i, part in enumerate(parts):
        if part in ['R1', 'R2', 'I1', 'I2']:
            file_type = part
            break
    
    print(f"Splitting {fastq_file}...")
    print(f"Output directory: {output_dir}")
    
    # Dictionary to store file handles for each barcode
    file_handles = {}
    barcode_counts = defaultdict(int)
    
    try:
        # First pass: collect all data by barcode (streaming to disk, not memory)
        for header, seq, plus, qual, barcode in fastq_parser(fastq_file):
            barcode_counts[barcode] += 1
            
            # Create output filename
            if file_type:
                output_filename = f"{barcode}_{file_type}_001.fastq.gz"
            else:
                output_filename = f"{barcode}.fastq.gz"
            
            output_path = os.path.join(output_dir, output_filename)
            
            # Open file handle if not already open (lazy opening)
            if barcode not in file_handles:
                file_handles[barcode] = gzip.open(output_path, 'wt')
                print(f"  Creating output: {output_filename}")
            
            # Write record to appropriate file
            fh = file_handles[barcode]
            fh.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
    
    finally:
        # Close all file handles
        for fh in file_handles.values():
            fh.close()
    
    # Print summary
    print(f"\n✅ Splitting complete!")
    print(f"\nBarcode summary:")
    total_reads = 0
    for barcode in sorted(barcode_counts.keys()):
        count = barcode_counts[barcode]
        total_reads += count
        print(f"  {barcode}: {count:,} reads")
    
    print(f"\nTotal reads: {total_reads:,}")
    print(f"Total barcodes: {len(barcode_counts)}")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Split FASTQ file by barcode"
    )
    parser.add_argument(
        "--fastq",
        required=True,
        help="Input FASTQ file (gzipped or uncompressed)"
    )
    parser.add_argument(
        "--out",
        default="split_output",
        help="Output directory (default: split_output)"
    )
    
    args = parser.parse_args()
    
    # Verify input file exists
    if not Path(args.fastq).exists():
        print(f"Error: {args.fastq} not found")
        sys.exit(1)
    
    split_fastq_by_barcode(args.fastq, args.out)
