import argparse
import gzip

def compare_fastq_reads(r1_file, r2_file, output_file='mismatched_reads.txt'):
    """
    Compare read names between R1 and R2 FASTQ files. 
    Report all mismatched read names with their line numbers.
    """
    mismatches = []
    
    with gzip.open(r1_file, 'rt') as f1, gzip.open(r2_file, 'rt') as f2:
        line_num = 0
        for line_r1, line_r2 in zip(f1, f2):
            line_num += 1
            
            if line_num % 4 == 1:
                read_name_r1 = line_r1.strip().split()[0]
                read_name_r2 = line_r2.strip().split()[0]
                
                if read_name_r1 != read_name_r2:
                    mismatches.append({
                        'line_num': line_num,
                        'r1_name': read_name_r1,
                        'r2_name': read_name_r2
                    })
    
    with open(output_file, 'w') as fout:
        fout.write("Line Number\tR1 Read Name\tR2 Read Name\n")
        for mismatch in mismatches:
            fout.write(f"{mismatch['line_num']}\t{mismatch['r1_name']}\t{mismatch['r2_name']}\n")
    print(f"Total mismatches found: {len(mismatches)}")
    print(f"Results saved to: {output_file}")
    if mismatches: 
        print("\nFirst 10 mismatches:")
        for mismatch in mismatches[:10]:
            print(f"Line {mismatch['line_num']}: {mismatch['r1_name']} vs {mismatch['r2_name']}")

def main():
    parser = argparse.ArgumentParser(description="compare paired FASTQ reads")
    parser.add_argument("r1", help="R1 FASTQ file")
    parser.add_argument("r2", help="R2 FASTQ file")
    parser.add_argument("-o", "--output-file", default="mismatched_reads.txt")

    args = parser.parse_args()

    compare_fastq_reads(
        args.r1,
        args.r2,
        output_file=args.output_file
    )

if __name__ == "__main__":
    main()
