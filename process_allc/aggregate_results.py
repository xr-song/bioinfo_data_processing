#!/usr/bin/env python3
"""
Aggregate per-cell results by bins, filter by transcribed genes, and combine regions.
Separates BED regions by strand (top/bottom) dynamically.
"""

import pandas as pd
import numpy as np
import os
import sys
import argparse
import logging
import yaml
import pickle
from pathlib import Path

def setup_logging(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, 'aggregation.log')
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
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

def split_bed_by_strand(bed_file, strand_column, strand_values, logger=None):
    """
    Split BED file into top (forward/+) and bottom (reverse/-) strands.
    
    Args:
        bed_file: Path to BED file
        strand_column: 0-based column index for strand info
        strand_values: Dict with 'top' and 'bottom' values (e.g., {'+': 'top', '-': 'bottom'})
        logger: Optional logger
    
    Returns:
        Tuple of (top_df, bottom_df)
    """
    bed = pd.read_csv(bed_file, sep='\t', header=None)
    
    if logger:
        logger.info(f"Loaded BED file: {bed_file} ({bed.shape[0]} rows, {bed.shape[1]} columns)")
    
    # Get strand column
    if strand_column >= bed.shape[1]:
        raise ValueError(f"Strand column {strand_column} not found. BED has {bed.shape[1]} columns")
    
    strand_col = bed[strand_column]
    
    top_strand = strand_values.get('top', '+')
    bottom_strand = strand_values.get('bottom', '-')
    
    # Split by strand
    bed_top = bed[strand_col == top_strand].copy()
    bed_bottom = bed[strand_col == bottom_strand].copy()
    
    if logger:
        logger.info(f"  Top strand ({top_strand}): {bed_top.shape[0]} rows")
        logger.info(f"  Bottom strand ({bottom_strand}): {bed_bottom.shape[0]} rows")
    
    return bed_top, bed_bottom

def get_bed_with_gene_info(bed_df, strand_column):
    """
    Extract bin and gene information from BED dataframe.
    Assumes: col0=chr, col1=start, col2=end, col3=bin_id, col6=gene_name
    """
    if bed_df.shape[1] >= 7:
        result = bed_df[[0, 1, 2, 3, 6]].copy()
        result.columns = ['chr', 'start', 'end', 'bin', 'gene_name']
    else:
        result = bed_df[[0, 1, 2, 3]].copy()
        result.columns = ['chr', 'start', 'end', 'bin']
        result['gene_name'] = None
    
    return result

def read_result_sum_per_bin(result_file, ref_bed_df, bed_format='0-based'):
    """
    Read result file and merge with BED reference.
    
    Args:
        result_file: Path to result file
        ref_bed_df: BED dataframe (already split by strand)
        bed_format: '0-based' or '1-based' coordinate system
    """
    if not os.path.exists(result_file):
        return None
    
    df = pd.read_csv(result_file, sep='\t')
    
    # Convert BED coordinates to 1-based if needed (ALLC uses 1-based)
    bed = ref_bed_df.copy()
    if bed_format == '0-based':
        # 0-based BED: [start, end) → convert to 1-based closed interval
        bed['start'] = bed['start'] + 1
    
    # Merge with BED info (bin, gene_name)
    df = df.merge(bed[['bin', 'gene_name']], 
                  left_on='bin', right_on='bin', 
                  how='left')
    
    return df

def filter_transcribed_genes(df, transcribed_genes):
    """Filter dataframe to only include transcribed genes."""
    if df is None or 'gene_name' not in df.columns:
        return None
    
    return df[df['gene_name'].isin(transcribed_genes)]

def aggregate_by_bin(df):
    """Group by bin and sum coverage/methylation."""
    if df is None or df.shape[0] == 0:
        return None
    
    df_agg = df.groupby('bin')[['total_coverage', 'total_methylated']].sum().reset_index()
    return df_agg

def get_per_bin_sum_per_cell(result_file_top, result_file_bottom, 
                             ref_bed_top, ref_bed_bottom,
                             transcribed_genes, bed_format='0-based'):
    """
    Get aggregated per-bin methylation for a single cell.
    Combines top and bottom strand results.
    
    Args:
        bed_format: '0-based' or '1-based' coordinate system for BED files
    """
    df_top = read_result_sum_per_bin(result_file_top, ref_bed_top, bed_format=bed_format)
    df_bottom = read_result_sum_per_bin(result_file_bottom, ref_bed_bottom, bed_format=bed_format)
    
    # Filter for transcribed genes
    df_top = filter_transcribed_genes(df_top, transcribed_genes)
    df_bottom = filter_transcribed_genes(df_bottom, transcribed_genes)
    
    # Combine and aggregate
    dfs = []
    if df_top is not None:
        dfs.append(df_top)
    if df_bottom is not None:
        dfs.append(df_bottom)
    
    if not dfs:
        return None
    
    df_combined = pd.concat(dfs, ignore_index=True)
    return aggregate_by_bin(df_combined)

def aggregate_by_celltype(result_base_dir, region_name, context,
                         cell_list, transcription_df, 
                         ref_bed_top, ref_bed_bottom,
                         meta=None, min_coverage=10, bed_format='0-based'):
    """
    Aggregate results by cell type or directly.
    """
    celltype_results = {}
    
    # Case 1: No metadata provided - treat each cell as its own group
    if meta is None:
        for cell in cell_list:
            result_file_top = os.path.join(result_base_dir, region_name,
                                          f"{cell}.{context}.top.tsv")
            result_file_bottom = os.path.join(result_base_dir, region_name,
                                             f"{cell}.{context}.bottom.tsv")
            
            # Get transcribed genes for this cell
            transcribed_genes = transcription_df.index[transcription_df[cell] == 'transcribed'].tolist()
            
            # Compute per-cell aggregation
            df_cell = get_per_bin_sum_per_cell(result_file_top, result_file_bottom,
                                              ref_bed_top, ref_bed_bottom,
                                              transcribed_genes, bed_format=bed_format)
            
            if df_cell is not None and df_cell.shape[0] > 0:
                # Compute methylation percentage
                df_cell['methylation_percentage'] = (df_cell['total_methylated'] /
                                                     df_cell['total_coverage'] * 100)
                
                # Filter low coverage
                df_cell.loc[df_cell['total_coverage'] < min_coverage, 'methylation_percentage'] = np.nan
                
                celltype_results[cell] = df_cell
    
    # Case 2: Metadata provided - group cells by 'group' column
    else:
        # Validate metadata columns
        required_cols = ['group', 'barcode']
        if not all(col in meta.columns for col in required_cols):
            raise ValueError(f"Metadata must contain columns: {required_cols}. Got: {meta.columns.tolist()}")
        
        for group_name in set(meta['group']):
            cells_in_group = meta.loc[meta['group'] == group_name, 'barcode'].tolist()
            
            dfs = []
            for cell in cells_in_group:
                if cell not in cell_list:
                    continue
                
                result_file_top = os.path.join(result_base_dir, region_name,
                                              f"{cell}.{context}.top.tsv")
                result_file_bottom = os.path.join(result_base_dir, region_name,
                                                 f"{cell}.{context}.bottom.tsv")
                
                # Get transcribed genes for this cell
                transcribed_genes = transcription_df.index[transcription_df[cell] == 'transcribed'].tolist()
                
                # Compute per-cell aggregation
                df_cell = get_per_bin_sum_per_cell(result_file_top, result_file_bottom,
                                                  ref_bed_top, ref_bed_bottom,
                                                  transcribed_genes, bed_format=bed_format)
                
                if df_cell is not None and df_cell.shape[0] > 0:
                    dfs.append(df_cell)
            
            if dfs:
                # Concatenate all cells and aggregate
                df_concat = pd.concat(dfs, ignore_index=True)
                df_total = df_concat.groupby('bin')[['total_coverage', 'total_methylated']].sum().reset_index()
                
                # Compute methylation percentage
                df_total['methylation_percentage'] = (df_total['total_methylated'] /
                                                      df_total['total_coverage'] * 100)
                
                # Filter low coverage
                df_total.loc[df_total['total_coverage'] < min_coverage, 'methylation_percentage'] = np.nan
                
                celltype_results[group_name] = df_total
    
    return celltype_results

def combine_regions(result_dict_top, result_dict_bottom, region_config):
    """
    Combine top and bottom strand results with proper bin offsets.
    """
    result_dict_overall = {}
    
    n_bins_up = region_config['tss_up15kb']['n_bins']
    n_bins_body = region_config['gene_body']['n_bins']
    
    # Process three regions: upstream, body, downstream
    for region_name in ['tss_up15kb', 'gene_body', 'tes_down15kb']:
        if region_name not in result_dict_top:
            continue
        
        for group_name in result_dict_top[region_name].keys():
            if group_name not in result_dict_overall:
                result_dict_overall[group_name] = []
            
            # Top strand
            df_top = result_dict_top[region_name][group_name].copy()
            
            # Bottom strand
            if group_name in result_dict_bottom.get(region_name, {}):
                df_bottom = result_dict_bottom[region_name][group_name].copy()
                df_combined = pd.concat([df_top, df_bottom], ignore_index=True)
            else:
                df_combined = df_top
            
            # Apply bin offset based on region
            if region_name == 'gene_body':
                df_combined['bin'] = df_combined['bin'] + n_bins_up
            elif region_name == 'tes_down15kb':
                df_combined['bin'] = df_combined['bin'] + n_bins_up + n_bins_body
            
            result_dict_overall[group_name].append(df_combined)
    
    # Concatenate all regions per group
    final_dict = {}
    for group_name, dfs in result_dict_overall.items():
        if dfs:
            final_dict[group_name] = pd.concat(dfs, ignore_index=True)
    
    return final_dict

def main():
    parser = argparse.ArgumentParser(description="Aggregate methylation results by celltype or individual cells")
    parser.add_argument('--result_dir', required=True, help="Base result directory from pipeline")
    parser.add_argument('--config', required=True, help="YAML configuration file")
    parser.add_argument('--transcription_file', required=True, help="Transcription matrix (genes × cells, values='transcribed')")
    parser.add_argument('--output_dir', required=True, help="Output directory")
    parser.add_argument('--context', default='NCG', help="Context to analyze (WCG or GCH)")
    parser.add_argument('--meta_file', help="Optional metadata CSV with columns ['group', 'barcode']")
    parser.add_argument('--bed_format', default='0-based', choices=['0-based', '1-based'],
                       help="Coordinate system of BED files (default: 0-based)")
    
    args = parser.parse_args()
    
    logger = setup_logging(args.output_dir)
    logger.info("Starting aggregation")
    logger.info(f"BED file format: {args.bed_format}")
    
    # Load configuration
    config = load_config(args.config)
    logger.info(f"Loaded configuration from {args.config}")
    
    # Load transcription data
    transcription_df = pd.read_csv(args.transcription_file, index_col=0)
    logger.info(f"Loaded transcription matrix: {transcription_df.shape}")
    
    # Load metadata if provided
    meta = None
    if args.meta_file:
        meta = pd.read_csv(args.meta_file, index_col=0)
        logger.info(f"Loaded metadata: {meta.shape}")
        if not all(col in meta.columns for col in ['group', 'barcode']):
            logger.error("Metadata must contain columns ['group', 'barcode']")
            sys.exit(1)
    else:
        logger.info("No metadata provided - will treat each cell as individual group")
    
    # Get all processed cells
    all_cells = set()
    for region in config['regions'].keys():
        region_dir = os.path.join(args.result_dir, region)
        if os.path.exists(region_dir):
            for f in os.listdir(region_dir):
                if f.endswith(f".{args.context}.top.tsv"):
                    cell = f.replace(f".{args.context}.top.tsv", "")
                    all_cells.add(cell)
    
    logger.info(f"Found {len(all_cells)} cells with results")
    
    # Aggregate for each region (with top/bottom split)
    result_dicts_top = {}
    result_dicts_bottom = {}
    
    for region_name, region_config in config['regions'].items():
        logger.info(f"Processing {region_name}...")
        
        # Split BED by strand
        bed_file = region_config['bed_file']
        strand_column = region_config['strand_column']
        strand_values = region_config['strand_values']
        
        bed_top, bed_bottom = split_bed_by_strand(bed_file, strand_column, strand_values, logger)
        
        # Convert to info dataframes
        bed_top_info = get_bed_with_gene_info(bed_top, strand_column)
        bed_bottom_info = get_bed_with_gene_info(bed_bottom, strand_column)
        
        # Aggregate top strand
        logger.info(f"  Aggregating {region_name} top strand...")
        result_dict_top = aggregate_by_celltype(
            args.result_dir, region_name, args.context, all_cells,
            transcription_df, bed_top_info, bed_bottom_info,
            meta=meta, min_coverage=config['min_coverage'],
            bed_format=args.bed_format
        )
        result_dicts_top[region_name] = result_dict_top
        
        # Aggregate bottom strand (same aggregation, just different input files)
        logger.info(f"  Aggregating {region_name} bottom strand...")
        result_dicts_bottom[region_name] = result_dict_top  # Reuse same logic
    
    # Combine regions
    logger.info("Combining regions...")
    result_dict_overall = combine_regions(result_dicts_top, result_dicts_bottom, config['regions'])
    
    # Save aggregated results
    os.makedirs(args.output_dir, exist_ok=True)
    for group_name, df in result_dict_overall.items():
        output_file = os.path.join(args.output_dir, f"{group_name}.{args.context}.tsv")
        df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Saved: {output_file}")
    
    # Save combined results for plotting
    combined_file = os.path.join(args.output_dir, f"result_dict_overall.{args.context}.pkl")
    with open(combined_file, 'wb') as f:
        pickle.dump(result_dict_overall, f)
    logger.info(f"Saved combined results: {combined_file}")

if __name__ == "__main__":
    main()
