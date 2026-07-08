#!/usr/bin/env python3
"""
Plotting utilities for methylation results.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.nonparametric.smoothers_lowess import lowess
import pickle

def plot_methylation_whole_gene_celltypes(result_dict_overall, palette, 
                                         x_key='bin', xlim=None, ylim=None,
                                         lowess_frac=0.02, prefix='meth.whole_gene',
                                         ylabel='% Methylation', output_dir='.'):
    """
    Plot methylation across whole gene for different cell types.
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    
    for celltype, df in result_dict_overall.items():
        # Get valid data
        df_valid = df[df['methylation_percentage'].notna()].copy()
        
        if df_valid.shape[0] == 0:
            continue
        
        # Apply LOWESS smoothing if requested
        if lowess_frac > 0:
            smoothed = lowess(df_valid['methylation_percentage'].values, 
                            df_valid[x_key].values, 
                            frac=lowess_frac)
            x_smooth = smoothed[:, 0]
            y_smooth = smoothed[:, 1]
            ax.plot(x_smooth, y_smooth, label=celltype, linewidth=2, 
                   color=palette.get(celltype, None))
        else:
            ax.plot(df_valid[x_key], df_valid['methylation_percentage'], 
                   label=celltype, marker='o', linewidth=2,
                   color=palette.get(celltype, None))
    
    # Formatting
    ax.set_xlabel('Position (bins)', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title('Methylation across whole gene', fontsize=14, fontweight='bold')
    
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    output_file = f"{output_dir}/{prefix}.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Plot methylation results")
    parser.add_argument('--result_pkl', required=True, help="Pickled result_dict_overall")
    parser.add_argument('--palette_file', help="JSON file with color palette")
    parser.add_argument('--output_dir', default='.', help="Output directory")
    parser.add_argument('--prefix', default='methylation', help="Output prefix")
    parser.add_argument('--ylabel', default='% Methylation', help="Y-axis label")
    parser.add_argument('--xlim', nargs=2, type=int, help="X-axis limits")
    parser.add_argument('--ylim', nargs=2, type=int, help="Y-axis limits")
    parser.add_argument('--lowess_frac', type=float, default=0.02, help="LOWESS smoothing fraction")
    
    args = parser.parse_args()
    
    # Load results
    with open(args.result_pkl, 'rb') as f:
        result_dict_overall = pickle.load(f)
    
    # Load palette
    palette = {}
    if args.palette_file:
        import json
        with open(args.palette_file, 'r') as f:
            palette = json.load(f)
    
    xlim = tuple(args.xlim) if args.xlim else None
    ylim = tuple(args.ylim) if args.ylim else None
    
    plot_methylation_whole_gene_celltypes(
        result_dict_overall, palette,
        xlim=xlim, ylim=ylim,
        lowess_frac=args.lowess_frac,
        prefix=args.prefix,
        ylabel=args.ylabel,
        output_dir=args.output_dir
    )

if __name__ == "__main__":
    main()
