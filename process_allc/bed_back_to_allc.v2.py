import argparse
import pandas as pd
import pybedtools


def bed_to_allc_with_cpg_filter(bed_file, cpg_bed, allc_output):
    # Read BED file with the new format: chrom, start, end, mc, cov
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, 
                         names=['chrom', 'start', 'end', 'mc', 'cov'])
    
    # Read CpG BED file
    cpg = pybedtools.BedTool(cpg_bed)
    cpg_df = pd.read_csv(cpg_bed, sep='\t', header=None,
                         names=['cpg_chrom', 'cpg_start', 'cpg_end'])
    
    records = []
    
    for _, row in bed_df.iterrows():
        chrom = row['chrom']
        start = int(row['start'])
        end = int(row['end'])
        mc = int(row['mc'])
        cov = int(row['cov'])
        
        # Find exact match in CpG regions (both start and end must match exactly)
        matching_cpg = cpg_df[
            (cpg_df['cpg_chrom'] == chrom) & 
            (cpg_df['cpg_start'] == start) & 
            (cpg_df['cpg_end'] == end)
        ]
        
        # Only keep if there's an exact match
        if not matching_cpg.empty:
            # Convert 0-based BED to 1-based position (use start + 1)
            records.append({
                'chrom': chrom,
                'pos': start + 1,
                'strand': '+',  # Default strand if not specified
                'context': 'CG',  # Assuming CpG context
                'mc': mc,
                'cov': cov,
                'flag': 1
            })
    
    if not records:
        raise RuntimeError("no records left after CpG exact match filtering")
    
    df = pd.DataFrame(records)
    
    # Remove duplicated genomic positions entirely
    dup_mask = df.duplicated(subset=['chrom', 'pos'], keep=False)
    df = df.loc[~dup_mask]
    n_duplicates = dup_mask.sum()
    print(f"Number of duplicated lines to remove: {n_duplicates}", flush=True)
    
    df = df.sort_values(['chrom', 'pos']).reset_index(drop=True)
    
    output_path = allc_output if allc_output.endswith('.gz') else f"{allc_output}.gz"
    df.to_csv(output_path, sep='\t', header=False, index=False, compression='gzip')
    
    print(f"Converted {bed_file} to ALLC format: {output_path}")
    print(f"Total records after filtering: {len(df)}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert BED (chrom, start, end, mc, cov) to ALLC, subset to exact CpG matches"
    )
    parser.add_argument('--bed', required=True, help='Input BED file (chrom, start, end, mc, cov)')
    parser.add_argument('--cpg_bed', required=True, help='CpG BED file for exact coordinate matching')
    parser.add_argument('--outallc', required=True, help='Output ALLC file')
    
    args = parser.parse_args()
    
    bed_to_allc_with_cpg_filter(args.bed, args.cpg_bed, args.outallc)


if __name__ == '__main__':
    main()
