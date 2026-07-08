#!/usr/bin/env python3

import pysam
import argparse
import os

parser = argparse.ArgumentParser(
    description="Split Bismark BAM into strand-specific BAMs and sort outputs"
)

parser.add_argument("--input_bam", required=True)
parser.add_argument("--output_dir", required=True)
parser.add_argument("--output4bams", action="store_true",
                    help="Output OT, OB, CTOT, CTOB as four BAMs")
parser.add_argument("--threads", type=int, default=1,
                    help="Threads for sorting")

args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

basename = os.path.basename(args.input_bam)
prefix = os.path.splitext(basename)[0]

bam_in = pysam.AlignmentFile(args.input_bam, "rb")

flag_to_strand = {
    99: "OT",
    83: "OB",
    147: "CTOT",
    163: "CTOB"
}

writers = {}
tmp_paths = {}
final_paths = {}

if args.output4bams:

    strands = ["OT", "OB", "CTOT", "CTOB"]

    for s in strands:
        tmp = os.path.join(args.output_dir, f"{prefix}.{s}.tmp.bam")
        final = os.path.join(args.output_dir, f"{prefix}.{s}.bam")

        tmp_paths[s] = tmp
        final_paths[s] = final
        writers[s] = pysam.AlignmentFile(tmp, "wb", template=bam_in)

else:

    groups = ["OT_OB", "CTOT_CTOB"]

    for g in groups:
        tmp = os.path.join(args.output_dir, f"{prefix}.{g}.tmp.bam")
        final = os.path.join(args.output_dir, f"{prefix}.{g}.bam")

        tmp_paths[g] = tmp
        final_paths[g] = final
        writers[g] = pysam.AlignmentFile(tmp, "wb", template=bam_in)

buffer = {}

for read in bam_in.fetch(until_eof=True):

    if read.is_unmapped:
        continue

    name = read.query_name

    if name not in buffer:
        buffer[name] = read
        continue

    first = buffer.pop(name)
    second = read

    strand = flag_to_strand.get(first.flag)

    if strand is None:
        continue

    if args.output4bams:

        writers[strand].write(first)
        writers[strand].write(second)

    else:

        if strand in ("OT", "OB"):
            writers["OT_OB"].write(first)
            writers["OT_OB"].write(second)
        else:
            writers["CTOT_CTOB"].write(first)
            writers["CTOT_CTOB"].write(second)

bam_in.close()

for w in writers.values():
    w.close()

print("Sorting outputs...")

for key in tmp_paths:

    pysam.sort(
        "-@", str(args.threads),
        "-o", final_paths[key],
        tmp_paths[key]
    )

    pysam.index(final_paths[key])

    os.remove(tmp_paths[key])

print("Done.")
print("Output files:")

for p in final_paths.values():
    print(p)
