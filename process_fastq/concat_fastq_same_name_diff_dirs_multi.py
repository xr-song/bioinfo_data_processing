#!/usr/bin/env python3
"""
Concatenate FASTQ files from multiple directories with guaranteed order matching.
Usage: python concatenate_fastq. py dir1 dir2 [dir3 ...]
"""

import sys
import os
import re
from pathlib import Path
from collections import defaultdict
import subprocess
from concurrent.futures import ProcessPoolExecutor
import logging

def setup_logging(log_file):
    """Setup logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging. StreamHandler(sys.stdout)
        ]
    )

def extract_sample_id(filename):
    """
    Extract sample ID from filename by removing lane, sample number, and read info.
    Examples: 
      Sample1_S1_L001_R1_001.fastq.gz -> Sample1
      Sample1_L001_R1.fastq.gz -> Sample1
      Sample1_R1.fastq.gz -> Sample1
    """
    name = Path(filename).name
    # Remove common suffixes
    name = re.sub(r'_S\d+', '', name)  # Remove _S1, _S2, etc.
    name = re.sub(r'_L\d+', '', name)  # Remove _L001, _L002, etc. 
    name = re.sub(r'_R[12].*', '', name)  # Remove _R1, _R2 and everything after
    name = re.sub(r'_I[12].*', '', name)  # Remove _I1, _I2 and everything after
    return name

def get_read_type(filename):
    """Determine if file is R1, R2, I1, or I2"""
    name = Path(filename).name
    if '_R1' in name or '_R1_' in name:
        return 'R1'
    elif '_R2' in name or '_R2_' in name:
        return 'R2'
    elif '_I1' in name or '_I1_' in name: 
        return 'I1'
    elif '_I2' in name or '_I2_' in name:
        return 'I2'
    return None

def collect_files(input_dirs):
    """
    Collect all FASTQ files and organize by sample and read type.
    Returns: dict[sample_id][read_type] = [sorted list of file paths]
    """
    files_by_sample = defaultdict(lambda: defaultdict(list))
    
    for input_dir in input_dirs: 
        input_path = Path(input_dir)
        if not input_path.exists():
            logging.warning(f"Directory does not exist: {input_dir}")
            continue
            
        # Find all fastq. gz files
        for fastq_file in input_path.rglob('*.fastq.gz'):
            sample_id = extract_sample_id(fastq_file. name)
            read_type = get_read_type(fastq_file.name)
            
            if read_type: 
                files_by_sample[sample_id][read_type]. append(str(fastq_file))
    
    # Sort files within each category to ensure consistent order
    for sample_id in files_by_sample: 
        for read_type in files_by_sample[sample_id]:
            files_by_sample[sample_id][read_type].sort()
    
    return files_by_sample

def verify_file_pairing(files_by_sample):
    """
    Verify that files are properly paired across read types.
    Returns list of warnings if issues found.
    """
    warnings = []
    
    for sample_id, read_types in files_by_sample.items():
        counts = {rt: len(files) for rt, files in read_types.items()}
        
        # Check if R1 and R2 have same number of files
        if 'R1' in counts and 'R2' in counts:
            if counts['R1'] != counts['R2']: 
                warnings.append(
                    f"Sample {sample_id}:  R1 has {counts['R1']} files, "
                    f"R2 has {counts['R2']} files - MISMATCH!"
                )
        
        # Check if I1 and I2 have same number of files (if both present)
        if 'I1' in counts and 'I2' in counts:
            if counts['I1'] != counts['I2']:
                warnings.append(
                    f"Sample {sample_id}: I1 has {counts['I1']} files, "
                    f"I2 has {counts['I2']} files - MISMATCH!"
                )
    
    return warnings

def concatenate_files(sample_id, read_type, input_files, output_file, log_file):
    """Concatenate files for a single sample/read_type combination"""
    try:
        # Write detailed log
        with open(log_file, 'a') as log:
            log.write(f"\n{'='*80}\n")
            log.write(f"Sample: {sample_id} | Read Type: {read_type}\n")
            log.write(f"Output: {output_file}\n")
            log.write(f"Number of input files: {len(input_files)}\n")
            log.write(f"Concatenation order:\n")
            for i, f in enumerate(input_files, 1):
                log.write(f"  {i}.  {f}\n")
            log.write(f"{'='*80}\n")
        
        # Concatenate using cat (faster than Python for this)
        with open(output_file, 'wb') as outf:
            for input_file in input_files: 
                with open(input_file, 'rb') as inf:
                    # Copy in chunks for memory efficiency
                    while True:
                        chunk = inf.read(1024 * 1024)  # 1MB chunks
                        if not chunk:
                            break
                        outf.write(chunk)
        
        logging.info(f"✓ Completed {sample_id}_{read_type}")
        return True
        
    except Exception as e:
        logging.error(f"✗ Failed {sample_id}_{read_type}: {str(e)}")
        return False

def process_sample(args):
    """Process a single sample (all read types) - wrapper for parallel execution"""
    sample_id, read_types_dict, output_dir, log_dir = args
    
    results = []
    log_file = log_dir / f"{sample_id}. log"
    
    for read_type in ['R1', 'R2', 'I1', 'I2']:
        if read_type in read_types_dict: 
            input_files = read_types_dict[read_type]
            output_file = output_dir / f"{sample_id}_{read_type}_001.fastq.gz"
            
            success = concatenate_files(
                sample_id, read_type, input_files, 
                str(output_file), str(log_file)
            )
            results.append((sample_id, read_type, success))
    
    return results

def main():
    if len(sys.argv) < 3:
        print("Usage: python concatenate_fastq.py dir1 dir2 [dir3 ...]")
        sys.exit(1)
    
    input_dirs = sys.argv[1:]
    output_dir = Path. cwd() / "rawdata_combined"
    log_dir = output_dir / "logs"
    
    # Create output directories
    output_dir.mkdir(exist_ok=True)
    log_dir.mkdir(exist_ok=True)
    
    # Setup logging
    setup_logging(log_dir / "concatenation_master.log")
    
    logging.info(f"Concatenating FASTQ files from {len(input_dirs)} directories")
    logging.info(f"Input directories: {input_dirs}")
    logging.info(f"Output directory: {output_dir}")
    
    # Collect all files
    logging.info("Collecting files...")
    files_by_sample = collect_files(input_dirs)
    
    if not files_by_sample: 
        logging.error("No FASTQ files found!")
        sys.exit(1)
    
    logging.info(f"Found {len(files_by_sample)} unique samples")
    
    # Verify pairing
    warnings = verify_file_pairing(files_by_sample)
    if warnings:
        logging.warning("File pairing issues detected:")
        for w in warnings: 
            logging.warning(f"  {w}")
    
    # Print summary
    logging.info("\nSample summary:")
    for sample_id in sorted(files_by_sample.keys()):
        read_counts = {rt: len(files) for rt, files in files_by_sample[sample_id].items()}
        logging.info(f"  {sample_id}: {read_counts}")
    
    # Process samples in parallel
    logging.info(f"\nProcessing samples with 6 parallel workers...")
    
    process_args = [
        (sample_id, read_types_dict, output_dir, log_dir)
        for sample_id, read_types_dict in files_by_sample. items()
    ]
    
    with ProcessPoolExecutor(max_workers=6) as executor:
        all_results = executor.map(process_sample, process_args)
    
    # Collect and report results
    success_count = 0
    fail_count = 0
    for sample_results in all_results:
        for sample_id, read_type, success in sample_results:
            if success:
                success_count += 1
            else:
                fail_count += 1
    
    logging.info(f"\n{'='*80}")
    logging.info(f"SUMMARY:  {success_count} successful, {fail_count} failed")
    logging.info(f"Output directory: {output_dir}")
    logging.info(f"Logs directory: {log_dir}")
    logging.info(f"{'='*80}")
    
    if fail_count > 0:
        sys.exit(1)

if __name__ == "__main__": 
    main()
