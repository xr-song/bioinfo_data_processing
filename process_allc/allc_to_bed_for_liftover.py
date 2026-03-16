import pandas as pd
import gzip
import argparse

def allc_to_bed(allc_file, bed_file):
    """Convert ALLC to BED format for LiftOver, packing all info in name field with _"""
    
    # Read ALLC file
    if allc_file.endswith('.gz'):
        allc_df = pd.read_csv(allc_file, sep='\t', header=None, 
                              names=['chrom', 'pos', 'strand', 'context', 'mc', 'cov', 'flag'],
                              compression='gzip')
    else:
        allc_df = pd.read_csv(allc_file, sep='\t', header=None,
                              names=['chrom', 'pos', 'strand', 'context', 'mc', 'cov', 'flag'])
    
    # Pack all info into name field with _ separator
    bed_df = pd.DataFrame({
        'chrom': allc_df['chrom'],
        'chromStart': allc_df['pos'] - 1,  # Convert to 0-based
        'chromEnd': allc_df['pos'],         # Half-open interval
        'strand': allc_df['strand'],
        'name': allc_df.apply(lambda row: f"{row['context']}_{row['mc']}_{row['cov']}", axis=1),
    })
    
    # Write BED file
    bed_df.to_csv(bed_file, sep='\t', header=False, index=False)
    print(f"Converted {allc_file} to {bed_file}")
    
    return bed_df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--allc', required=True)
    parser.add_argument('--outbed', required=True)
    args = parser.parse_args()
    allc_to_bed(args.allc, args.outbed)

if __name__ == '__main__':
    main()
