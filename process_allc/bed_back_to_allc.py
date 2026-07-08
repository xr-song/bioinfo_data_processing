import argparse
import pandas as pd
import pybedtools


def bed_to_allc_with_cpg_filter(bed_file, cpg_bed, allc_output):
    bed = pybedtools.BedTool(bed_file)
    cpg = pybedtools.BedTool(cpg_bed)

    # keep only BED entries overlapping CpG regions
    intersected = bed.intersect(cpg, u=True)

    records = []
    for interval in intersected:
        try:
            chrom, start, end, strand, name = interval.fields[:5]
            context, mc, cov = name.split('_')
            records.append({
                'chrom': chrom,
                'pos': int(start) + 1,
                'strand': strand,
                'context': context,
                'mc': int(mc),
                'cov': int(cov),
                'flag': 1
            })
        except (ValueError, IndexError):
            continue

    if not records:
        raise RuntimeError("no records left after CpG intersection")

    df = pd.DataFrame(records)

    # remove duplicated genomic positions entirely
    dup_mask = df.duplicated(subset=['chrom', 'pos'], keep=False)  # False : Mark all duplicates as True
    df = df.loc[~dup_mask]
    n_duplicates = dup_mask.sum()
    print(f"Number of duplicated lines to remove: {n_duplicates}", flush=True)

    df = df.sort_values(['chrom', 'pos']).reset_index(drop=True)

    output_path = allc_output if allc_output.endswith('.gz') else f"{allc_output}.gz"
    df.to_csv(output_path, sep='\t', header=False, index=False, compression='gzip')
    
    print(f"Converted {bed_file} to ALLC format: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert lifted BED to ALLC, subset to CpG regions, drop duplicated sites"
    )
    parser.add_argument('--bed', required=True, help='Input BED file')
    parser.add_argument('--cpg_bed', required=True, help='CpG BED file')
    parser.add_argument('--outallc', required=True, help='Output ALLC file')

    args = parser.parse_args()

    bed_to_allc_with_cpg_filter(args.bed, args.cpg_bed, args.outallc)


if __name__ == '__main__':
    main()

