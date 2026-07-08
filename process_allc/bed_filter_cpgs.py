import argparse
import pandas as pd
import pybedtools


def bed_to_filtered_bed_with_cpg_filter(bed_file, cpg_bed, bed_output, summary_output):
    # Read BED file with the new format: chrom, start, end, mc, cov
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, 
                         names=['chrom', 'start', 'end', 'mc', 'cov'])
    
    initial_count = len(bed_df)
    
    # Remove rows with NaN values
    bed_df_clean = bed_df.dropna()
    n_nan_removed = initial_count - len(bed_df_clean)
    
    if n_nan_removed > 0:
        print(f"Removed {n_nan_removed} entries with NaN values")
    
    if len(bed_df_clean) == 0:
        raise RuntimeError("No valid records left after removing NaN entries")
    
    # Convert to integer to ensure proper types
    bed_df_clean = bed_df_clean.copy()
    bed_df_clean['start'] = bed_df_clean['start'].astype(int)
    bed_df_clean['end'] = bed_df_clean['end'].astype(int)
    bed_df_clean['mc'] = bed_df_clean['mc'].astype(int)
    bed_df_clean['cov'] = bed_df_clean['cov'].astype(int)
    
    # Convert BED dataframe to pybedtools format
    bed_for_tool = bed_df_clean[['chrom', 'start', 'end']].copy()
    bed_for_tool['mc'] = bed_df_clean['mc'].astype(str)
    bed_for_tool['cov'] = bed_df_clean['cov'].astype(str)
    
    bed_tool = pybedtools.BedTool.from_dataframe(bed_for_tool)
    cpg_tool = pybedtools.BedTool(cpg_bed)
    
    # Intersect with exact match: -f 1.0 means 100% overlap on both query and subject
    intersected = bed_tool.intersect(cpg_tool, f=1.0, F=1.0, wa=True)
    
    # Parse results back to dataframe
    records = []
    for interval in intersected:
        records.append({
            'chrom': interval.chrom,
            'start': int(interval.start),
            'end': int(interval.end),
            'mc': int(interval.fields[3]),
            'cov': int(interval.fields[4])
        })
    
    if not records:
        print("Warning: No records left after CpG exact match filtering")
        filtered_df = pd.DataFrame(columns=['chrom', 'start', 'end', 'mc', 'cov'])
    else:
        filtered_df = pd.DataFrame(records)
    
    final_count = len(filtered_df)
    n_filtered_out = len(bed_df_clean) - final_count
    
    # Sort by chrom and position
    filtered_df = filtered_df.sort_values(['chrom', 'start']).reset_index(drop=True)
    
    # Write filtered BED file
    filtered_df.to_csv(bed_output, sep='\t', header=False, index=False)
    print(f"Wrote filtered BED file: {bed_output}")
    
    # Write summary file
    with open(summary_output, 'w') as f:
        f.write(f"Input CpGs: {initial_count}\n")
        f.write(f"CpGs with NaN values removed: {n_nan_removed}\n")
        f.write(f"CpGs before CpG region filtering: {len(bed_df_clean)}\n")
        f.write(f"CpGs after CpG region filtering: {final_count}\n")
        f.write(f"CpGs filtered out by CpG regions: {n_filtered_out}\n")
        retention_rate = (final_count / initial_count * 100) if initial_count > 0 else 0
        f.write(f"Overall retention rate: {retention_rate:.2f}%\n")
    
    print(f"Wrote summary file: {summary_output}")
    print(f"Input CpGs: {initial_count}")
    print(f"CpGs with NaN values removed: {n_nan_removed}")
    print(f"CpGs after CpG region filtering: {final_count}")
    print(f"CpGs filtered out by CpG regions: {n_filtered_out}")


def main():
    parser = argparse.ArgumentParser(
        description="Filter BED (chrom, start, end, mc, cov) to exact CpG matches using pybedtools"
    )
    parser.add_argument('--bed', required=True, help='Input BED file (chrom, start, end, mc, cov)')
    parser.add_argument('--cpg_bed', required=True, help='CpG BED file for exact coordinate matching')
    parser.add_argument('--outbed', required=True, help='Output filtered BED file')
    parser.add_argument('--summary', required=True, help='Output summary file')
    
    args = parser.parse_args()
    
    bed_to_filtered_bed_with_cpg_filter(args.bed, args.cpg_bed, args.outbed, args.summary)


if __name__ == '__main__':
    main()
