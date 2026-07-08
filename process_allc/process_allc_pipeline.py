#!/usr/bin/env python3
"""
Main pipeline for computing methylation across genomic regions.
Handles parallel processing of ALLC files and separate processing for forward/reverse strands.
"""

import pandas as pd
import numpy as np
import os
import sys
import argparse
import gzip
import logging
import yaml
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
from functools import partial

def setup_logging(output_dir, name='pipeline.log'):
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, name)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

def load_config(config_file):
    """Load configuration from YAML file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def load_allc_file(allc_file):
    """Load ALLC file into memory (chr, pos) -> (methylated, coverage)."""
    allc_data = {}
    line_count = 0
    
    open_func = gzip.open if allc_file.endswith('.gz') else open
    with open_func(allc_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            line_count += 1
            if line_count % 1000000 == 0:
                print(f"  Loaded {line_count:,} positions", file=sys.stderr)
            
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                chr_name = str(parts[0])
                pos = int(parts[1])
                methylated = int(parts[4])
                covered = int(parts[5])
                allc_data[(chr_name, pos)] = (methylated, covered)
    
    return allc_data

def load_bed_file(bed_file, strand_column, strand_values={'top':'+', 'bottom':'-'}):
    """
    Load BED file and split by strand (top/bottom).
    
    Args:
        bed_file: Path to BED file
        strand_column: 1-based column index for strand info
        strand_values: Dict with 'top' and 'bottom' values
    
    Returns:
        Tuple of (top_bed, bottom_bed) DataFrames
    """
    bed = pd.read_csv(bed_file, sep='\t', header=None)
    strand_column = strand_column-1
    
    if strand_column >= bed.shape[1]:
        raise ValueError(f"Strand column {strand_column} not found. BED has {bed.shape[1]} columns")
    
    strand_col = bed[strand_column]
    
    top_strand = strand_values.get('top', '+')
    bottom_strand = strand_values.get('bottom', '-')
    
    # Split by strand
    bed_top = bed[strand_col == top_strand].copy()
    bed_bottom = bed[strand_col == bottom_strand].copy()
    
    return bed_top, bed_bottom

def compute_region_methylation(allc_data, bed_df):
    """
    Compute methylation per region/bin from a BED dataframe.
    Assumes: col0=chr, col1=start, col2=end, col3=bin_id
    """
    results = []
    
    for _, row in bed_df.iterrows():
        chr_name = str(row[0])
        start = int(row[1])
        end = int(row[2])
        bin_id = row[3]
        
        total_coverage = 0
        total_methylated = 0
        sites_with_data = 0
        
        for pos in range(start, end + 1):
            if (chr_name, pos) in allc_data:
                methylated, covered = allc_data[(chr_name, pos)]
                total_coverage += covered
                total_methylated += methylated
                sites_with_data += 1
        
        results.append({
            'region_id': f"{chr_name}:{start}-{end}",
            'bin': bin_id,
            'total_coverage': total_coverage,
            'total_methylated': total_methylated,
            'sites_with_data': sites_with_data
        })
    
    return pd.DataFrame(results)

def process_single_file(allc_file, config, output_dir, logger):
    """Process a single ALLC file for all regions and contexts (top/bottom strands)."""
    barcode = Path(allc_file).stem
    # Remove context suffix to get clean barcode
    for context in config['contexts']:
        barcode = barcode.replace(f".{context}", "")
    
    logger.info(f"Processing: {allc_file} (barcode: {barcode})")
    
    # Load ALLC data once
    allc_data = load_allc_file(allc_file)
    logger.info(f"  Loaded {len(allc_data):,} methylation positions")
    
    # Extract context from filename
    filename = Path(allc_file).name
    context = None
    for ctx in config['contexts']:
        if ctx in filename:
            context = ctx.split('-')[0]  # Get 'WCG' or 'GCH' prefix
            break
    
    if context is None:
        logger.warning(f"  Could not detect context in {filename}")
        return
    
    # Process each region
    for region_name, region_config in config['regions'].items():
        logger.info(f"  Processing region: {region_name}")
        
        bed_file = region_config['bed_file']
        strand_column = region_config['strand_column']
        strand_values = strand_values={'top':'+', 'bottom':'-'} #region_config['strand_values']
        
        # Split BED by strand
        bed_top, bed_bottom = load_bed_file(bed_file, strand_column, strand_values)
        
        # Compute methylation for top strand
        df_top = compute_region_methylation(allc_data, bed_top)
        output_file_top = os.path.join(output_dir, region_name,
                                       f"{barcode}.{context}.top.tsv")
        os.makedirs(os.path.dirname(output_file_top), exist_ok=True)
        df_top.to_csv(output_file_top, sep='\t', index=False)
        logger.info(f"    Saved top: {output_file_top}")
        
        # Compute methylation for bottom strand
        df_bottom = compute_region_methylation(allc_data, bed_bottom)
        output_file_bottom = os.path.join(output_dir, region_name,
                                          f"{barcode}.{context}.bottom.tsv")
        df_bottom.to_csv(output_file_bottom, sep='\t', index=False)
        logger.info(f"    Saved bottom: {output_file_bottom}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute methylation across genomic regions from ALLC files"
    )
    parser.add_argument('--allc_list', required=True, 
                       help="Text file with ALLC file paths (one per line)")
    parser.add_argument('--config', required=True, 
                       help="YAML configuration file")
    parser.add_argument('--output_dir', required=True, 
                       help="Output directory for results")
    parser.add_argument('--threads', type=int, default=8, 
                       help="Number of parallel threads")
    
    args = parser.parse_args()
    
    # Setup
    logger = setup_logging(args.output_dir)
    logger.info(f"Starting pipeline with {args.threads} threads")
    
    config = load_config(args.config)
    logger.info(f"Loaded configuration from {args.config}")
    logger.info(f"Regions: {', '.join(config['regions'].keys())}")
    logger.info(f"Contexts: {', '.join(config['contexts'])}")
    
    # Read ALLC file list
    with open(args.allc_list, 'r') as f:
        allc_files = [line.strip() for line in f if line.strip()]
    
    logger.info(f"Found {len(allc_files)} ALLC files")
    
    # Process in parallel
    process_func = partial(process_single_file, config=config, 
                          output_dir=args.output_dir, logger=logger)
    
    with Pool(args.threads) as pool:
        pool.map(process_func, allc_files)
    
    logger.info("Pipeline complete!")

if __name__ == "__main__":
    main()
