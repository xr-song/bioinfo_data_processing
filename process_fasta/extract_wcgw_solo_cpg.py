#!/usr/bin/env python3
"""
Extract solo WCGW CpG sites from a genome FASTA file and output in BED format.
"""

import argparse
import sys
import re
import os
import time
from datetime import datetime
from multiprocessing import Pool
from bisect import bisect_left, bisect_right
from Bio import SeqIO


def log_message(message, level="INFO"):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("[{}] {:8s} | {}".format(timestamp, level, message), file=sys.stdout, flush=True)


def log_error(message):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("[{}] {:8s} | {}".format(timestamp, "ERROR", message), file=sys.stderr, flush=True)


def find_wcgw_solo_sites(sequence, chromosome, solo_distance=1000):
    """
    Find WCGW CpG sites that are solo using binary search — O(n log n) instead of O(n*m).
    """
    seq_upper = sequence.upper()

    all_cpg_positions = [m.start() for m in re.finditer(r'CG', seq_upper)]

    if not all_cpg_positions:
        return [], 0, 0

    wcgw_cpg_positions = []
    for match in re.finditer(r'[AT]CG[AT]', seq_upper):
        cpg_start = match.start() + 1  # CpG starts after W
        wcgw_cpg_positions.append(cpg_start)

    total_wcgw = len(wcgw_cpg_positions)

    solo_wcgw_sites = []
    for cpg_pos in wcgw_cpg_positions:
        lo = bisect_left(all_cpg_positions, cpg_pos - solo_distance)
        hi = bisect_right(all_cpg_positions, cpg_pos + solo_distance)

        neighbors = hi - lo - 1  # subtract 1 for self

        if neighbors == 0:
            solo_wcgw_sites.append((chromosome, cpg_pos, cpg_pos + 2))

    return solo_wcgw_sites, total_wcgw, len(solo_wcgw_sites)


def process_chromosome(args):
    """
    Worker function: process one chromosome, write to tmp file, return stats.
    """
    chromosome, sequence, tmpdir, solo_distance = args

    tmp_file = os.path.join(tmpdir, "{}.bed".format(chromosome))

    t0 = time.time()
    log_message("Starting: {} ({:,} bp)".format(chromosome, len(sequence)), level="DEBUG")

    try:
        solo_wcgw_sites, total_wcgw, solo_count = find_wcgw_solo_sites(
            sequence, chromosome, solo_distance=solo_distance
        )

        with open(tmp_file, 'w') as f:
            for chr_name, start, end in solo_wcgw_sites:
                f.write("{}\t{}\t{}\n".format(chr_name, start, end))

        elapsed = time.time() - t0
        log_message(
            "Done: {} | WCGW: {:,} | Solo: {:,} | {:.2f}s".format(
                chromosome, total_wcgw, solo_count, elapsed),
            level="INFO"
        )

        return tmp_file, chromosome, total_wcgw, solo_count

    except Exception as e:
        log_error("Error processing {}: {}".format(chromosome, e))
        raise


def process_fasta_file(input_file, output_file, solo_distance=1000, max_workers=None, tmpdir='./tmp'):
    """
    Process FASTA file using multiprocessing + binary search.
    """
    if max_workers is None:
        max_workers = os.cpu_count()

    os.makedirs(tmpdir, exist_ok=True)

    start_time = time.time()

    log_message("=" * 70)
    log_message("SOLO WCGW CpG SITE EXTRACTION STARTED")
    log_message("=" * 70)
    log_message("Input:          {}".format(input_file))
    log_message("Output:         {}".format(output_file))
    log_message("Solo distance:  {} bp".format(solo_distance))
    log_message("Workers:        {}".format(max_workers))
    log_message("Tmp dir:        {}".format(tmpdir))

    log_message("Loading FASTA records...")
    records = list(SeqIO.parse(input_file, "fasta"))
    log_message("Found {} chromosomes/contigs".format(len(records)))

    args_list = [
        (record.id, str(record.seq), tmpdir, solo_distance)
        for record in records
    ]

    log_message("Processing with {} cores...".format(max_workers))
    with Pool(max_workers) as pool:
        results = pool.map(process_chromosome, args_list)

    results.sort(key=lambda x: x[1])
    log_message("Merging results into: {}".format(output_file))
    total_wcgw = 0
    total_solo = 0
    chromosome_stats = []

    with open(output_file, 'w') as fout:
        for tmp_file, chromosome, wcgw_count, solo_count in results:
            total_wcgw += wcgw_count
            total_solo += solo_count
            chromosome_stats.append((chromosome, wcgw_count, solo_count))
            if os.path.exists(tmp_file):
                with open(tmp_file, 'r') as fin:
                    fout.write(fin.read())
                os.remove(tmp_file)

    total_time = time.time() - start_time
    log_message("=" * 70)
    log_message("SUMMARY")
    log_message("=" * 70)
    for chromosome, wcgw_count, solo_count in chromosome_stats:
        pct = (solo_count / wcgw_count * 100) if wcgw_count > 0 else 0
        log_message("{:>15} | WCGW: {:>10,} | Solo: {:>10,} ({:>5.1f}%)".format(
            chromosome, wcgw_count, solo_count, pct))
    log_message("-" * 70)
    solo_pct = (total_solo / total_wcgw * 100) if total_wcgw > 0 else 0
    log_message("{:>15} | WCGW: {:>10,} | Solo: {:>10,} ({:>5.1f}%)".format(
        "TOTAL", total_wcgw, total_solo, solo_pct))
    log_message("-" * 70)
    log_message("Output: {}".format(output_file))
    log_message("Total time: {:.2f}s".format(total_time))
    log_message("=" * 70)
    log_message("COMPLETED SUCCESSFULLY")
    log_message("=" * 70)


def main():
    parser = argparse.ArgumentParser(
        description="Extract solo WCGW CpG sites from genome FASTA file to BED format"
    )
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output BED file')
    parser.add_argument('-d', '--distance', type=int, default=35,
                        help='Solo distance threshold in bp (default: 35)')
    parser.add_argument('-t', '--threads', type=int, default=None,
                        help='Number of worker processes (default: all cores)')
    parser.add_argument('--tmpdir', default='./tmp',
                        help='Temporary directory (default: ./tmp)')

    args = parser.parse_args()

    if not os.path.exists(args.input):
        log_error("Input file '{}' not found.".format(args.input))
        sys.exit(1)

    process_fasta_file(
        args.input, args.output,
        solo_distance=args.distance,
        max_workers=args.threads,
        tmpdir=args.tmpdir
    )


if __name__ == "__main__":
    main()
