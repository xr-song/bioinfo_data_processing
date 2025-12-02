import os
import glob
import subprocess
from collections import defaultdict
import argparse
from concurrent.futures import ProcessPoolExecutor

def merge_single_readtype(args):
    prefix, read_type, input_files, output_dir, s_num = args
    s_string = f"S{s_num:02d}"
    output_filename = f"{prefix}_{s_string}_{read_type}_001.fastq.gz"
    output_path = os.path.join(output_dir, output_filename)
    
    print(f"  Merging {len(input_files)} {read_type} files -> {output_filename}")
    
    try:
        with open(output_path, 'wb') as outfile:
            for input_file in sorted(input_files):
                print(f"    Adding: {os.path.basename(input_file)}")
                with open(input_file, 'rb') as infile:
                    outfile.write(infile.read())
        
        print(f"  ✓ Created: {output_filename}")
        return True
        
    except Exception as e:
        print(f"  ✗ Error merging {read_type} files for {prefix}: {e}")
        return False

def merge_fastq_files(input_dir, output_dir, num_cores):
    os.makedirs(output_dir, exist_ok=True)

    fastq_files = glob.glob(os.path.join(input_dir, "*.fastq.gz"))

    if not fastq_files:
        print("No FASTQ files found in the directory")
        return

    file_groups = defaultdict(lambda: defaultdict(list))

    for file_path in fastq_files:
        filename = os.path.basename(file_path)

        if "-" not in filename:
            print(f"Warning: No '-' found in {filename}, skipping...")
            continue

        prefix = filename.split("-")[0]

        if "_I1_" in filename:
            read_type = "I1"
        elif "_I2_" in filename:
            read_type = "I2"
        elif "_R1_" in filename:
            read_type = "R1"
        elif "_R2_" in filename:
            read_type = "R2"
        else:
            print(f"Warning: Cannot determine read type for {filename}, skipping...")
            continue

        file_groups[prefix][read_type].append(file_path)

    sorted_prefixes = sorted(file_groups.keys())

    tasks = []
    for s_num, prefix in enumerate(sorted_prefixes, 1):
        print(f"Processing {prefix} -> S{s_num:02d}")
        
        for read_type in ["I1", "I2", "R1", "R2"]:
            if read_type in file_groups[prefix]:
                input_files = file_groups[prefix][read_type]
                tasks.append((prefix, read_type, input_files, output_dir, s_num))

    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        results = list(executor.map(merge_single_readtype, tasks))

def main():
    parser = argparse.ArgumentParser(description="Merge FASTQ files by prefix")
    parser.add_argument("-i", "--input_dir", required=True, help="Input directory containing FASTQ files")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for merged files")
    parser.add_argument("-p", "--cores", type=int, default=4, help="Number of CPU cores to use")
    
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"Error: {args.input_dir} is not a valid directory")
        exit(1)

    print(f"Processing FASTQ files in: {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Using {args.cores} cores")
    merge_fastq_files(args.input_dir, args.output_dir, args.cores)
    print("Done!")

if __name__ == "__main__":
    main()
