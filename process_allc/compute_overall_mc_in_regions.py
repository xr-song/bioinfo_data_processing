import subprocess
import pandas as pd
import numpy as np
import os
import argparse
import logging
import sys
from pathlib import Path
from multiprocessing import Pool

def setup_logging(log_file='region_methylation.log'):
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )


def allc_to_bed(allc_file, tmpdir):
    """
    Convert ALLC to BED format (0-based)
    ALLC: chr, pos(1-based), strand, mc_class, methylated, covered, ...
    BED:  chr, start(0-based), end, mc_class, methylated, covered
    """
    sample_name = Path(allc_file).name.replace('.allc.tsv.gz', '').replace('.allc.tsv', '')
    bed_file = os.path.join(tmpdir, "{}.bed".format(sample_name))

    logging.info("Converting ALLC to BED: {}".format(sample_name))

    import gzip as gz
    opener = gz.open if allc_file.endswith('.gz') else open

    with opener(allc_file, 'rt') as fin, open(bed_file, 'w') as fout:
        for line in fin:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 6:
                chrom = fields[0]
                pos = int(fields[1])
                start = pos - 1  # convert 1-based to 0-based
                end = pos
                mc_class = fields[3]
                methylated = fields[4]
                covered = fields[5]
                fout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    chrom, start, end, mc_class, methylated, covered))

    return bed_file


def process_allc_file(args):
    allc_file, regions_bed_file, output_dir, tmpdir = args
    sample_name = Path(allc_file).name.replace('.allc.tsv.gz', '').replace('.allc.tsv', '')
    output_file = os.path.join(output_dir, "{}.methylation.tsv".format(sample_name))

    logging.info("Processing: {}".format(sample_name))

    bed_file = None

    try:
        # Convert ALLC to BED first
        bed_file = allc_to_bed(allc_file, tmpdir)

        # Run bedtools intersect
        cmd = [
            'bedtools', 'intersect',
            '-a', bed_file,
            '-b', regions_bed_file,
            '-wa'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        total_methylated = 0
        total_coverage = 0
        sites_with_data = 0

        for line in result.stdout.splitlines():
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 6:
                methylated = int(fields[4])
                covered = int(fields[5])
                if covered > 0:
                    total_methylated += methylated
                    total_coverage += covered
                    sites_with_data += 1

        methylation_fraction = total_methylated / total_coverage if total_coverage > 0 else np.nan

        df = pd.DataFrame([{
            'total_mC': total_methylated,
            'total_cov': total_coverage,
            'methylation_fraction': methylation_fraction,
            'sites_with_data': sites_with_data
        }])
        df.to_csv(output_file, sep='\t', index=False, float_format='%.6f')

        logging.info("✓ {}: mC={:,}, cov={:,}, frac={:.6f}".format(
            sample_name, total_methylated, total_coverage, methylation_fraction))

        return sample_name, total_methylated, total_coverage, methylation_fraction

    except subprocess.CalledProcessError as e:
        logging.error("✗ bedtools error for {}: {}".format(sample_name, e.stderr))
        return sample_name, None, None, None
    except Exception as e:
        logging.error("✗ Error processing {}: {}".format(sample_name, e))
        return sample_name, None, None, None

    finally:
        # Clean up tmp bed file
        if bed_file and os.path.exists(bed_file):
            os.remove(bed_file)
            logging.info("Cleaned up tmp BED: {}".format(bed_file))

def main(allc_dir, regions_bed_file, output_dir, num_cores, tmpdir):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(tmpdir, exist_ok=True)
    setup_logging(os.path.join(output_dir, 'region_methylation.log'))

    logging.info(f"Using tmpdir: {tmpdir}")

    # Check bedtools is available
    try:
        subprocess.run(['bedtools', '--version'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        logging.error("bedtools not found in PATH")
        sys.exit(1)

    allc_files = sorted(Path(allc_dir).glob('*.allc.tsv.gz'))
    logging.info(f"Found {len(allc_files)} ALLC files")

    if not allc_files:
        logging.error(f"No *.allc.tsv.gz files found in {allc_dir}")
        return

    args_list = [(str(f), regions_bed_file, output_dir, tmpdir) for f in allc_files]

    logging.info(f"Processing with {num_cores} cores")

    with Pool(num_cores) as pool:
        results = pool.map(process_allc_file, args_list)

    valid = [(s, mc, cov, frac) for s, mc, cov, frac in results if mc is not None]
    total_mC = sum(mc for _, mc, _, _ in valid)
    total_cov = sum(cov for _, _, cov, _ in valid)

    logging.info("=" * 60)
    logging.info("Summary Statistics:")
    logging.info("Files processed: {}/{}".format(len(valid), len(results)))
    logging.info("Aggregate mC: {:,}".format(total_mC))
    logging.info("Aggregate coverage: {:,}".format(total_cov))

    if total_cov > 0:
        logging.info("Overall methylation fraction: {:.6f}".format(total_mC / total_cov))
    else:
        logging.warning("Overall methylation fraction: N/A (zero total coverage — check bedtools intersect output)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute aggregate methylation within genomic regions from multiple ALLC files"
    )
    parser.add_argument('--allc-dir', required=True,
                       help="Directory containing *.allc.tsv.gz files")
    parser.add_argument('--regions', required=True,
                       help="BED file with region coordinates (chr, start, end)")
    parser.add_argument('--outdir', required=True,
                       help="Output directory for TSV files (one per sample)")
    parser.add_argument('--cores', type=int, default=4,
                       help="Number of cores for parallel processing (default: 4)")
    parser.add_argument('--tmpdir', default='./tmp',
                       help="Temporary directory (default: ./tmp)")

    args = parser.parse_args()
    main(args.allc_dir, args.regions, args.outdir, args.cores, args.tmpdir)
