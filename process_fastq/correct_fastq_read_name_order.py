import gzip
import argparse
from collections import defaultdict

def process_fastq_files(i1_file, i2_file, r1_file, r2_file, output_dir='sorted_fastqs'):
    """
    Process 4 FASTQ files: 
    1. Write matching reads directly (group 1)
    2. Sort and write unmatched reads (group 2)
    3. Concatenate both groups
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Read all files and track which reads are matched
    print("Reading files...")
    reads = defaultdict(lambda: {'I1': None, 'I2': None, 'R1': None, 'R2':  None})
    
    # Parse each file
    for file_path, file_type in [(i1_file, 'I1'), (i2_file, 'I2'), 
                                   (r1_file, 'R1'), (r2_file, 'R2')]:
        with gzip.open(file_path, 'rt') as f:
            lines = f.readlines()
            for i in range(0, len(lines), 4):
                read_name = lines[i].strip().split()[0]
                reads[read_name][file_type] = ''.join(lines[i:i+4])
    
    print(f"Total unique read names: {len(reads)}")
    
    # Step 2:  Identify matched and unmatched reads
    matched_reads = []
    unmatched_reads = []
    
    for read_name, data in reads.items():
        if all(data. values()):  # All 4 files present
            matched_reads.append(read_name)
        else:
            unmatched_reads.append(read_name)
    
    print(f"Matched reads (all 4 files): {len(matched_reads)}")
    print(f"Unmatched reads: {len(unmatched_reads)}")
    
    # Step 3: Write group 1 (matched reads, in order they appeared)
    print("Writing group 1 (matched reads)...")
    temp_i1_g1 = f'{output_dir}/.temp_I1_g1.fastq.gz'
    temp_i2_g1 = f'{output_dir}/.temp_I2_g1.fastq.gz'
    temp_r1_g1 = f'{output_dir}/.temp_R1_g1.fastq.gz'
    temp_r2_g1 = f'{output_dir}/.temp_R2_g1.fastq.gz'
    
    with gzip.open(temp_i1_g1, 'wt') as f_i1, \
         gzip.open(temp_i2_g1, 'wt') as f_i2, \
         gzip. open(temp_r1_g1, 'wt') as f_r1, \
         gzip.open(temp_r2_g1, 'wt') as f_r2:
        for read_name in matched_reads: 
            f_i1.write(reads[read_name]['I1'])
            f_i2.write(reads[read_name]['I2'])
            f_r1.write(reads[read_name]['R1'])
            f_r2.write(reads[read_name]['R2'])
    
    print(f"Group 1 written:  {len(matched_reads)} reads")
    
    # Step 4: Write group 2 (unmatched reads, sorted)
    print("Writing group 2 (unmatched reads, sorted)...")
    temp_i1_g2 = f'{output_dir}/.temp_I1_g2.fastq.gz'
    temp_i2_g2 = f'{output_dir}/.temp_I2_g2.fastq.gz'
    temp_r1_g2 = f'{output_dir}/.temp_R1_g2.fastq.gz'
    temp_r2_g2 = f'{output_dir}/.temp_R2_g2.fastq.gz'
    
    sorted_unmatched = sorted(unmatched_reads)
    
    with gzip. open(temp_i1_g2, 'wt') as f_i1, \
         gzip.open(temp_i2_g2, 'wt') as f_i2, \
         gzip.open(temp_r1_g2, 'wt') as f_r1, \
         gzip.open(temp_r2_g2, 'wt') as f_r2:
        for read_name in sorted_unmatched:
            if reads[read_name]['I1']: 
                f_i1.write(reads[read_name]['I1'])
            if reads[read_name]['I2']:
                f_i2.write(reads[read_name]['I2'])
            if reads[read_name]['R1']:
                f_r1.write(reads[read_name]['R1'])
            if reads[read_name]['R2']:
                f_r2.write(reads[read_name]['R2'])
    
    print(f"Group 2 written: {len(unmatched_reads)} reads (sorted)")
    
    # Step 5: Concatenate group 1 + group 2
    print("Concatenating group 1 + group 2...")
    
    output_files = {
        'I1': f'{output_dir}/I1_sorted.fastq.gz',
        'I2': f'{output_dir}/I2_sorted.fastq.gz',
        'R1': f'{output_dir}/R1_sorted.fastq.gz',
        'R2':  f'{output_dir}/R2_sorted.fastq.gz'
    }
    
    for file_type, output_path in output_files.items():
        temp_g1 = f'{output_dir}/.temp_{file_type}_g1.fastq.gz'
        temp_g2 = f'{output_dir}/.temp_{file_type}_g2.fastq.gz'
        
        with gzip.open(output_path, 'wt') as fout:
            # Write group 1
            with gzip.open(temp_g1, 'rt') as fin:
                fout.write(fin.read())
            # Write group 2
            with gzip.open(temp_g2, 'rt') as fin:
                fout.write(fin.read())
        
        # Remove temp files
        os.remove(temp_g1)
        os.remove(temp_g2)
    
    print(f"\nFinal sorted files written to: {output_dir}/")
    print(f"  - I1_sorted.fastq.gz ({len(matched_reads)} matched + {len(unmatched_reads)} unmatched)")
    print(f"  - I2_sorted.fastq.gz")
    print(f"  - R1_sorted.fastq. gz")
    print(f"  - R2_sorted.fastq.gz")

def main():
    parser = argparse.ArgumentParser(
        description='Sort 4 FASTQ files (I1, I2, R1, R2) by matching read names'
    )
    parser.add_argument('--i1', required=True, help='Path to I1 FASTQ file')
    parser.add_argument('--i2', required=True, help='Path to I2 FASTQ file')
    parser.add_argument('--r1', required=True, help='Path to R1 FASTQ file')
    parser.add_argument('--r2', required=True, help='Path to R2 FASTQ file')
    parser.add_argument('--output-dir', default='sorted_fastqs', 
                        help='Output directory (default: sorted_fastqs)')
    
    args = parser.parse_args()
    
    process_fastq_files(args.i1, args.i2, args.r1, args. r2, args.output_dir)

if __name__ == '__main__':
    main()
