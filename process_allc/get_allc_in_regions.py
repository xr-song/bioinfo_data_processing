#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
import argparse
import gzip
import logging
import sys
from collections import defaultdict

def setup_logging(log_file='region_methylation.log'):
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

def load_pseudobulk_allc(allc_file):
    logging.info(f"Loading ALLC file: {allc_file}")

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

    logging.info(f"Loaded {len(allc_data):,} positions from ALLC")
    return allc_data

def compute_region_methylation(regions_bed_file, allc_data):
    logging.info(f"Loading regions from: {regions_bed_file}")

    regions = pd.read_csv(regions_bed_file, sep='\t', header=None)
    regions = regions[[0, 1, 2]]
    regions.columns = ['chr', 'start', 'end']
    regions['chr'] = regions['chr'].astype(str)
    regions['start'] = regions['start'].astype(int)
    regions['end'] = regions['end'].astype(int)

    logging.info(f"Loaded {len(regions)} regions")

    results = []

    for region_idx, (_, row) in enumerate(regions.iterrows()):
        if region_idx % 1000 == 0 and region_idx > 0:
            logging.info(f"  Processed {region_idx:,}/{len(regions)} regions")

        chr_name = row['chr']
        start = row['start']
        end = row['end']

        total_coverage = 0
        total_methylated = 0
        sites_with_data = 0

        # Check all positions within the region
        for pos in range(start, end + 1):
            chr_pos = (chr_name, pos)
            if chr_pos in allc_data:
                methylated, covered = allc_data[chr_pos]
                total_coverage += covered
                total_methylated += methylated
                sites_with_data += 1

        # Calculate methylation percentage
        if total_coverage > 0:
            methylation_percentage = (total_methylated / total_coverage) * 100
        else:
            methylation_percentage = np.nan

        results.append({
            'chr': chr_name,
            'start': start,
            'end': end,
            'region_length': end - start + 1,
            'total_coverage': total_coverage,
            'total_methylated': total_methylated,
            'methylation_percentage': methylation_percentage,
            'sites_with_data': sites_with_data
        })

    logging.info("Region methylation computation complete")

    return results

def write_results(results, output_file):
    logging.info(f"Writing results to {output_file}")

    df = pd.DataFrame(results)

    df.to_csv(output_file, sep='\t', index=False, float_format='%.3f')
    logging.info(f"Results written to {output_file}")

def main(pseudobulk_allc_file, regions_bed_file, output_file):
    log_dir = os.path.dirname(output_file) if os.path.dirname(output_file) else '.'
    setup_logging(os.path.join(log_dir, 'region_methylation.log'))

    allc_data = load_pseudobulk_allc(pseudobulk_allc_file)

    results = compute_region_methylation(regions_bed_file, allc_data)

    write_results(results, output_file)

    # Summary statistics
    df_results = pd.DataFrame(results)
    regions_with_data = (df_results['sites_with_data'] > 0).sum()
    avg_coverage = df_results[df_results['total_coverage'] > 0]['total_coverage'].mean()
    avg_methylation = df_results['methylation_percentage'].mean()

    logging.info(f"Summary:")
    logging.info(f"  Regions with data: {regions_with_data}/{len(results)}")
    logging.info(f"  Average coverage per region: {avg_coverage:.1f}")
    logging.info(f"  Average methylation percentage: {avg_methylation:.2f}%")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute methylation within genomic regions from ALLC file")
    parser.add_argument('--allc', required=True, help="ALLC file (.allc.tsv.gz)")
    parser.add_argument('--regions', required=True, help="BED file with region coordinates (chr, start, end)")
    parser.add_argument('--output', required=True, help="Output TSV file for results")

    args = parser.parse_args()
    main(args.allc, args.regions, args.output)
