import pysam
import gzip
import argparse

def main(input_bam, output_bed_gz, frag_summary_tsv):
    inbam = pysam.AlignmentFile(input_bam, "rb")
    read_buffer = {}
    out = gzip.open(output_bed_gz, "wt")
    frag_min, frag_max = 20, 500
    frag_counts = {size: 0 for size in range(frag_min, frag_max+1)}
    for read in inbam.fetch(until_eof=True):
        if not read.is_paired or read.is_unmapped or read.mate_is_unmapped:
            continue
        if read.query_name not in read_buffer:
            read_buffer[read.query_name] = read
        else:
            mate = read_buffer.pop(read.query_name)
            if read.is_read1:
                read1, read2 = read, mate
            else:
                read1, read2 = mate, read
            if read1.reference_id != read2.reference_id:
                continue
            chrom = read1.reference_name
            # forward strand pair: R1 +, R2 -
            if not read1.is_reverse and read2.is_reverse:
                start = read1.reference_start
                end = read2.reference_end
                strand = '+'
            # reverse strand pair: R2 +, R1 -
            elif read1.is_reverse and not read2.is_reverse:
                start = read2.reference_start
                end = read1.reference_end
                strand = '-'
            else:
                continue
            if start >= end:
                continue

            frag_size = end - start

            #dovetail = False
            #r1_start, r1_end = read1.reference_start, read1.reference_end
            #r2_start, r2_end = read2.reference_start, read2.reference_end
            #if not read1.is_reverse and read2.is_reverse:
            #    if (r2_end < r1_end) or (r2_start < r1_start):
            #        dovetail = True
            #elif read1.is_reverse and not read2.is_reverse:
            #    if (r1_end < r2_end) or (r1_start < r2_start):
            #        dovetail = True
            #if dovetail:
            #    overlap_start = max(r1_start, r2_start)
            #    overlap_end = min(r1_end, r2_end)
            #    frag_size = overlap_end - overlap_start

            # bed file
            out.write(f"{chrom}\t{start}\t{end}\t{strand}\t{frag_size}\n")
            # fragment size summary
            if frag_min <= frag_size <= frag_max:
                frag_counts[frag_size] += 1

    out.close()
    inbam.close()

    with open(frag_summary_tsv, "w") as fs:
        print("fragment_size\tcount", file=fs)
        for size in range(frag_min, frag_max + 1):
            print(f"{size}\t{frag_counts[size]}", file=fs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract fragment info from paired-end BAM and output as BED (0-based) and fragment size summary.")
    parser.add_argument("input_bam", help="BAM file")
    parser.add_argument("output_bed_gz", help="Output .bed.gz")
    parser.add_argument("--frag-summary", required=True, help="Fragment size summary tsv output")
    args = parser.parse_args()
    main(args.input_bam, args.output_bed_gz, args.frag_summary)
