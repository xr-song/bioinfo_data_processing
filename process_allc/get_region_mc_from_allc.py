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

def load_allc(allc_file):
    logging.info(f"Loading ALLC file: {allc_file}")

    allc_data = {}
    line_count = 0

    if allc_file. endswith('.gz'):
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

    logging. info(f"Loaded {len(allc_data):,} positions from ALLC")
    return allc_data

def load_regions(bed_file, bed_is_zero_based=True):
    logging.info(f"Loading regions from: {bed_file}")
    logging.info(f"BED coordinate system: {'0-based [start, end)' if bed_is_zero_based else '1-based [start, end]'}")
    
    regions = pd.read_csv(bed_file, sep='\t', header=None)
    
    if len(regions. columns) < 3:
        raise ValueError("BED file must have at least 3 columns: chr, start, end")
    
    regions = regions.iloc[:, :3]
    regions.columns = ['chr', 'start', 'end']
    regions['chr'] = regions['chr'].astype(str)
    regions['start'] = regions['start'].astype(int)
    regions['end'] = regions['end'].astype(int)
    
    regions['zero_based'] = bed_is_zero_based
    
    regions['region_id'] = regions['chr'] + '_' + regions['start'].astype(str) + '_' + regions['end'].astype(str)
    
    logging. info(f"Loaded {len(regions)} regions")
    return regions

def compute_region_methylation(regions, allc_data):
    """Compute methylation statistics for each region

    Handles coordinate conversion between BED (0-based or 1-based) and ALLC (1-based)
    """
    logging.info(f"Computing methylation for {len(regions)} regions")

    results = []
    bed_is_zero_based = regions['zero_based'].iloc[0]

    for region_idx, (_, row) in enumerate(regions.iterrows()):
        if region_idx % 1000 == 0 and region_idx > 0:
            logging.info(f"  Processed {region_idx:,}/{len(regions)} regions")

        chr_name = row['chr']
        start = row['start']
        end = row['end']
        region_id = row['region_id']

        total_coverage = 0
        total_methylated = 0
        positions_with_data = 0

        # Convert BED coordinates to 1-based ALLC coordinates
        if bed_is_zero_based:
            # BED is 0-based [start, end), convert to 1-based [start+1, end]
            # Example: BED [100, 200) covers 0-based positions 100-199
            #          which correspond to 1-based positions 101-200
            allc_start = start + 1
            allc_end = end  # end is exclusive in BED, becomes inclusive in 1-based
        else:
            # BED is already 1-based [start, end]
            allc_start = start
            allc_end = end

        # Query ALLC data using 1-based coordinates
        for pos in range(allc_start, allc_end + 1):
            chr_pos = (chr_name, pos)

            if chr_pos in allc_data:
                methylated, covered = allc_data[chr_pos]
                total_coverage += covered
                total_methylated += methylated
                positions_with_data += 1

        # Calculate methylation percentage
        if total_coverage > 0:
            methylation_percentage = (total_methylated / total_coverage) * 100
        else:
            methylation_percentage = np.nan

        results.append({
            'region_id': region_id,
            'chr': chr_name,
            'start': start,
            'end': end,
            'total_coverage': total_coverage,
            'total_methylated': total_methylated,
            'methylation_percentage': methylation_percentage,
            'sites_with_data': positions_with_data
        })

    logging.info("Region methylation computation complete")

    return pd.DataFrame(results)

def write_results(results_df, output_file):
    logging.info(f"Writing results to {output_file}")

    # Set region_id as index
    results_df = results_df.set_index('region_id')
    
    results_df.to_csv(output_file, sep='\t', float_format='%.3f')
    logging.info(f"Results written to {output_file}")

def main(allc_file, bed_file, output_file, bed_is_zero_based=True):
    log_dir = os.path. dirname(output_file) if os.path.dirname(output_file) else '.'
    setup_logging(os.path.join(log_dir, 'region_methylation.log'))

    allc_data = load_allc(allc_file)

    regions = load_regions(bed_file, bed_is_zero_based)

    results_df = compute_region_methylation(regions, allc_data)

    write_results(results_df, output_file)

    regions_with_data = (results_df['sites_with_data'] > 0).sum()
    avg_coverage = results_df[results_df['total_coverage'] > 0]['total_coverage'].mean()
    avg_methylation = results_df['methylation_percentage'].mean()

    logging.info(f"Summary:")
    logging.info(f"  Regions with data: {regions_with_data}/{len(results_df)}")
    logging. info(f"  Average coverage per region: {avg_coverage:.1f}")
    logging.info(f"  Average methylation percentage: {avg_methylation:.2f}%")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute methylation statistics for genomic regions from ALLC file")
    parser.add_argument('--allc', required=True, help="ALLC file (.allc.tsv.gz or .allc.tsv)")
    parser.add_argument('--bed_file', required=True, help="BED file with genomic regions (chr, start, end)")
    parser. add_argument('--output_file', required=True, help="Output TSV file for results")
    parser. add_argument('--one-based', action='store_true', help="BED file uses 1-based coordinates (default: 0-based)")

    args = parser.parse_args()
    main(args.allc, args.bed_file, args.output_file, bed_is_zero_based=not args.one_based)
