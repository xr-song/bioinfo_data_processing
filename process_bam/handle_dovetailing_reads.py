import pysam
import argparse
import os

def mask_to_overlap(read, overlap_start, overlap_end):
    positions = read.get_reference_positions(full_length=True)
    quals = list(read.query_qualities)
    for i, pos in enumerate(positions):
        if pos is None or not (overlap_start <= pos < overlap_end):
            quals[i] = 0
    read.query_qualities = quals

def mask_dovetail_extruding_bases(input_bam, output_bam, write_dovetail_bam=False):
    inbam = pysam.AlignmentFile(input_bam, "rb")
    outbam = pysam.AlignmentFile(output_bam, "wb", template=inbam)

    dovetail_bam = None
    if write_dovetail_bam:
        dt_bam_name = input_bam.replace(".bam", ".dovetailing_reads.bam")
        dovetail_bam = pysam.AlignmentFile(dt_bam_name, "wb", template=inbam)

    read_buffer = {}

    for read in inbam.fetch(until_eof=True):
        if not read.is_paired or read.is_unmapped or read.mate_is_unmapped:
            outbam.write(read)
            continue
        qname = read.query_name
        if qname not in read_buffer:
            read_buffer[qname] = read
        else:
            mate = read_buffer.pop(qname)
            if read.is_read1:
                read1, read2 = read, mate
            else:
                read1, read2 = mate, read

            # If not same reference, ignore (won't have this case anyway)
            if read1.reference_id != read2.reference_id:
                continue

            r1_start, r1_end = read1.reference_start, read1.reference_end
            r2_start, r2_end = read2.reference_start, read2.reference_end

            dovetail = False

            if not read1.is_reverse and read2.is_reverse:
                if (r2_end < r1_end) or (r2_start < r1_start):
                    dovetail = True
            elif read1.is_reverse and not read2.is_reverse:
                if (r1_end < r2_end) or (r1_start < r2_start):
                    dovetail = True

            if dovetail:
                overlap_start = max(r1_start, r2_start)
                overlap_end = min(r1_end, r2_end)
                mask_to_overlap(read1, overlap_start, overlap_end)
                mask_to_overlap(read2, overlap_start, overlap_end)

            outbam.write(read1)
            outbam.write(read2)

            if dovetail and dovetail_bam:
                dovetail_bam.write(read1)
                dovetail_bam.write(read2)

    # Write unpaired leftovers
    for leftover in read_buffer.values():
        outbam.write(leftover)

    inbam.close()
    outbam.close()
    if dovetail_bam:
        dovetail_bam.close()

def _mask_overhang(read, ref_a, ref_b, is_left_overhang=False):
    """Set base quality to zero for portion of `read` aligning from ref_a (inclusive) to ref_b (exclusive).
    If is_left_overhang, mask the left side; else the right side.
    """
    if ref_a >= ref_b:
        return
    ref_positions = read.get_reference_positions(full_length=True)
    to_mask = [
        idx for idx, rp in enumerate(ref_positions)
        if rp is not None and ref_a <= rp < ref_b
    ]
    if not to_mask:
        return
    quals = list(read.query_qualities)
    for idx in to_mask:
        quals[idx] = 0
    read.query_qualities = quals

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mask base qualities in dovetail overhangs. Optionally write out those pairs to a separate BAM.")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bam", help="Output BAM file (qualities masked in dovetail overhangs)")
    parser.add_argument("--write_dovetail_bam", action="store_true",
                        help="If set, write modified dovetailing read pairs to a second BAM (input.bam -> input.dovetailing_reads.bam)")
    args = parser.parse_args()
    mask_dovetail_extruding_bases(args.input_bam, args.output_bam, args.write_dovetail_bam)
