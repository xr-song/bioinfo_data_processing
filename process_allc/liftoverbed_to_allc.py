import pandas as pd
import argparse

def bed_to_allc(bed_file, allc_output):
    """Convert LiftOver BED back to ALLC format"""
    
    # Read BED file
    bed_df = pd.read_csv(bed_file, sep='\t', header=None,
                         names=['chrom', 'chromStart', 'chromEnd', 'strand', 'name'])
    
    # Extract packed info from name field
    allc_data = []
    for idx, row in bed_df.iterrows():
        try:
            # Parse name: context_mc_cov
            parts = row['name'].split('_')
            context = parts[0]
            mc = int(parts[1])
            cov = int(parts[2])
            
            allc_data.append({
                'chrom': row['chrom'],
                'pos': int(row['chromStart']) + 1,  # Convert back to 1-based
                'strand': row['strand'],
                'context': context,
                'mc': mc,
                'cov': cov,
                'flag': 1  # Default flag value
            })
        except (ValueError, IndexError) as e:
            print(f"Warning: Could not parse line {idx}: {row['name']}")
            continue
    
    allc_df = pd.DataFrame(allc_data)
    
    # Sort by position for consistency
    allc_df = allc_df.sort_values(['chrom', 'pos']).reset_index(drop=True)
    
    # Write ALLC file (gzipped)
    output_path = allc_output if allc_output.endswith('.gz') else f"{allc_output}.gz"
    allc_df.to_csv(output_path, sep='\t', header=False, index=False, compression='gzip')
    print(f"Converted {bed_file} back to ALLC format: {output_path}")
    
    return allc_df

def main():
    parser = argparse.ArgumentParser(description='Convert BED format back to ALLC format after liftOver')
    parser.add_argument('--bed', required=True, help='Input BED file path')
    parser.add_argument('--outallc', required=True, help='Output ALLC file path')
    
    args = parser.parse_args()
    
    bed_to_allc(args.bed, args.outallc)

if __name__ == '__main__':
    main()
