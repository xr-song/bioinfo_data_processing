import argparse
import gzip
import os

def detect_gzip(filename):
    """Check if file is gzipped"""
    return filename.endswith('.gz')

def open_file(filename, mode='rt'):
    """Open file, handling gzip automatically"""
    if detect_gzip(filename):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def convert_beta_to_allc(input_beta, output_allc):
    """Convert beta format to allc format"""
    
    print(f"Reading beta file: {input_beta}")
    input_is_gzipped = detect_gzip(input_beta)
    output_is_gzipped = detect_gzip(output_allc)
    
    line_count = 0
    
    # Open output file
    with open_file(output_allc, 'wt') as out_f:
        # Open input file
        with open_file(input_beta, 'rt') as in_f:
            for line in in_f:
                line = line.strip()
                if not line:
                    continue
                
                cols = line.split('\t')
                
                if len(cols) < 5:
                    print(f"Warning: Skipping line with insufficient columns: {line}")
                    continue
                
                try:
                    chrom = cols[0]
                    pos_end = int(cols[2])  # 3rd col (end position)
                    mc = int(cols[3])        # 4th col (methylated cytosines)
                    cov = int(cols[4])       # 5th col (coverage)
                    
                    # Convert to allc format
                    pos = pos_end - 1        # 2nd col in allc = 3rd col in beta - 1
                    strand = "+"             # Always +
                    context = "CG"           # Always CG
                    flag = 1                 # Always 1
                    
                    # Write allc format: chr, pos, strand, context, mc, cov, flag
                    out_line = f"{chrom}\t{pos}\t{strand}\t{context}\t{mc}\t{cov}\t{flag}\n"
                    out_f.write(out_line)
                    
                    line_count += 1
                    
                except (ValueError, IndexError) as e:
                    print(f"Warning: Error processing line: {line}")
                    print(f"  Error: {e}")
                    continue
    
    print(f"✓ Successfully converted {line_count} lines")
    print(f"✓ Output written to: {output_allc}")
    print(f"  Input gzipped: {input_is_gzipped}")
    print(f"  Output gzipped: {output_is_gzipped}")

def main():
    parser = argparse.ArgumentParser(description='Convert beta format to allc format')
    
    parser.add_argument('--input_beta', required=True, help='Input beta file (.tsv or .tsv.gz)')
    parser.add_argument('--output_allc', required=True, help='Output allc file (.tsv or .tsv.gz)')
    
    args = parser.parse_args()
    
    # Validate input file exists
    if not os.path.exists(args.input_beta):
        print(f"Error: Input file not found: {args.input_beta}")
        exit(1)
    
    convert_beta_to_allc(args.input_beta, args.output_allc)

if __name__ == '__main__':
    main()
