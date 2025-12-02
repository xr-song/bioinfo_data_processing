import pandas as pd
import numpy as np
import os
import argparse
import gzip
import logging
import sys
from collections import defaultdict

def setup_logging(log_file='tss_methylation.log'):
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

def load_pseudobulk_allc(allc_file):
    logging.info(f"Loading pseudobulk ALLC file: {allc_file}")
    
    allc_data = {}
    line_count = 0
    
    # Handle both compressed and uncompressed files
    if allc_file.endswith('.gz'):
        file_handle = gzip.open(allc_file, 'rt')
    else:
        file_handle = open(allc_file, 'r')
    
    with file_handle:
        for line in file_handle:
            if line.startswith('#'):
                continue
            
            line_count += 1
            if line_count % 1000000 == 0:
                logging.info(f"  Loaded {line_count:,} positions")
            
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                chr_name = str(parts[0])
                pos = int(parts[1])
                methylated = int(parts[4])
                covered = int(parts[5])
                
                allc_data[(chr_name, pos)] = (methylated, covered)
    
    logging.info(f"Loaded {len(allc_data):,} positions from pseudobulk ALLC")
    return allc_data

def compute_tss_methylation(tss_bed_file, allc_data, range_around_tss=None):
    if range_around_tss is None:
        range_around_tss = list(range(-1000, 1001))
    
    logging.info(f"Computing methylation for {len(range_around_tss)} positions around TSS")
    
    tss_bed = pd.read_csv(tss_bed_file, sep='\t', header=None)
    tss_bed = tss_bed[[0, 1]]
    tss_bed.columns = ['chr', 'tss_start']
    tss_bed['chr'] = tss_bed['chr'].astype(str)
    
    logging.info(f"Loaded {len(tss_bed)} TSS sites")
    
    n_positions = len(range_around_tss)
    total_coverage = np.zeros(n_positions)
    total_methylated = np.zeros(n_positions)
    position_counts = np.zeros(n_positions) 
    
    for tss_idx, (_, row) in enumerate(tss_bed.iterrows()):
        if tss_idx % 1000 == 0 and tss_idx > 0:
            logging.info(f"  Processed {tss_idx:,}/{len(tss_bed)} TSS sites")
        
        chr_name = row['chr']
        tss_start = row['tss_start']
        
        for offset_idx, offset in enumerate(range_around_tss):
            pos = tss_start + offset
            chr_pos = (chr_name, pos)
            
            if chr_pos in allc_data:
                methylated, covered = allc_data[chr_pos]
                total_coverage[offset_idx] += covered
                total_methylated[offset_idx] += methylated
                position_counts[offset_idx] += 1
    
    methylation_percentages = np.zeros(n_positions)
    for i in range(n_positions):
        if total_coverage[i] > 0:
            methylation_percentages[i] = (total_methylated[i] / total_coverage[i]) * 100
        else:
            methylation_percentages[i] = np.nan
    
    logging.info("TSS methylation computation complete")
    
    return {
        'positions': range_around_tss,
        'coverage': total_coverage.astype(int),
        'methylated': total_methylated.astype(int),
        'methylation_percentage': methylation_percentages,
        'position_counts': position_counts.astype(int)
    }

def write_results(results, output_file):
    logging.info(f"Writing results to {output_file}")
    
    df = pd.DataFrame({
        'position_relative_to_tss': results['positions'],
        'total_coverage': results['coverage'],
        'total_methylated': results['methylated'],
        'methylation_percentage': results['methylation_percentage'],
        'sites_with_data': results['position_counts']
    })
    
    df.to_csv(output_file, sep='\t', index=False, float_format='%.3f')
    logging.info(f"Results written to {output_file}")

def main(pseudobulk_allc_file, tss_bed_file, output_file, range_start=-1000, range_end=1000):
    log_dir = os.path.dirname(output_file) if os.path.dirname(output_file) else '.'
    setup_logging(os.path.join(log_dir, 'tss_methylation.log'))
    
    range_around_tss = list(range(range_start, range_end + 1))
    logging.info(f"Analyzing range {range_start} to {range_end} around TSS ({len(range_around_tss)} positions)")
    
    allc_data = load_pseudobulk_allc(pseudobulk_allc_file)
    
    results = compute_tss_methylation(tss_bed_file, allc_data, range_around_tss)
    
    write_results(results, output_file)
    
    total_positions_with_data = np.sum(results['position_counts'] > 0)
    avg_coverage = np.nanmean(results['coverage'][results['coverage'] > 0])
    avg_methylation = np.nanmean(results['methylation_percentage'])
    
    logging.info(f"Summary:")
    logging.info(f"  Positions with data: {total_positions_with_data}/{len(range_around_tss)}")
    logging.info(f"  Average coverage per position: {avg_coverage:.1f}")
    logging.info(f"  Average methylation percentage: {avg_methylation:.2f}%")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute methylation around TSS from ALLC file")
    parser.add_argument('--allc', required=True, help="ALLC file (.allc.tsv.gz)")
    parser.add_argument('--tss_bed_file', required=True, help="BED file with TSS coordinates")
    parser.add_argument('--output_file', required=True, help="Output TSV file for results")
    parser.add_argument('--range_start', type=int, default=-1000, help='Start of range around TSS (default: -1000)')
    parser.add_argument('--range_end', type=int, default=1000, help='End of range around TSS (default: 1000)')

    args = parser.parse_args()
    main(args.pseudobulk_allc, args.tss_bed_file, args.output_file, args.range_start, args.range_end)
