import pysam

def is_dovetail(read1, read2):
    if read1.reference_id != read2.reference_id:
        return False

    r1_start, r1_end = read1.reference_start, read1.reference_end
    r2_start, r2_end = read2.reference_start, read2.reference_end

    # Dovetail: R1's end extends past R2's start AND R2's end extends past R1's start
    # need to account for strand orientation!
    dovetail=False

    if not read1.is_reverse and read2.is_reverse:
        if (r2_end < r1_end) or (r2_start < r1_start):
            dovetail = True
    elif read1.is_reverse and not read2.is_reverse:
        if (r1_end < r2_end) or (r1_start < r2_start):
            dovetail = True
    return dovetail

def extract_dovetail(input_bam, output_bam):
    inbam = pysam.AlignmentFile(input_bam, "rb")
    outbam = pysam.AlignmentFile(output_bam, "wb", template=inbam)
    read_buffer = {}
    for read in inbam.fetch(until_eof=True):
        if not read.is_paired or read.is_unmapped or read.mate_is_unmapped:
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
            if is_dovetail(read1, read2):
                outbam.write(read1)
                outbam.write(read2)
    inbam.close()
    outbam.close()

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python extract_dovetailing_reads.py input.bam output.bam")
        sys.exit(1)
    extract_dovetail(sys.argv[1], sys.argv[2])
