#!/usr/bin/env python3
import argparse
import gzip
import os
import sys
import logging

import numpy as np
import pandas as pd


def setup_logging(log_file=None, verbose=True):
    handlers = [logging.StreamHandler(sys.stdout)] if verbose else []
    if log_file:
        os.makedirs(os.path.dirname(log_file) or ".", exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=handlers
    )


def read_allc_as_df(allc_path):
    """
    Read ALLC as DataFrame with columns: chr, pos, meth, cov.
    ALLC pos is assumed 1-based.
    """
    logging.info(f"Reading ALLC: {allc_path}")
    opener = gzip.open if allc_path.endswith(".gz") else open
    # allc: chr, pos, strand?, context?, meth, cov, ...
    # We only need 0,1,4,5
    df = pd.read_csv(
        opener(allc_path, "rt"),
        sep="\t",
        header=None,
        comment="#",
        usecols=[0, 1, 4, 5],
        names=["chr", "pos", "meth", "cov"],
        dtype={"chr": "string", "pos": np.int64, "meth": np.int64, "cov": np.int64},
    )
    logging.info(f"ALLC rows loaded: {len(df):,}")
    return df


def read_tss_bed_split_strand(bed_path, strand_col_1based=5, top_value="+", bottom_value="-"):
    """
    Read a BED-like file and split into top/bottom by strand column.
    We only need chr and tss_start (assumed to be column 2 in your old script).
    Here we support generic bed with at least:
      col1=chr, col2=tss_start, col{strand_col}=strand

    strand_col_1based: e.g., 5 means 5th column.
    """
    logging.info(f"Reading BED: {bed_path}")
    bed = pd.read_csv(bed_path, sep="\t", header=None, dtype={0: "string"})
    strand_idx0 = strand_col_1based - 1
    if strand_idx0 >= bed.shape[1]:
        raise ValueError(f"strand_col={strand_col_1based} out of range; bed has {bed.shape[1]} columns")

    # keep chr (0) and tss_start (1) for compatibility with your existing bed
    # If your bed has start in column 1 (0-based), adjust here.
    tss = bed[[0, 1, strand_idx0]].copy()
    tss.columns = ["chr", "tss_start", "strand"]
    tss["tss_start"] = tss["tss_start"].astype(np.int64)

    top = tss[tss["strand"] == top_value][["chr", "tss_start"]].reset_index(drop=True)
    bottom = tss[tss["strand"] == bottom_value][["chr", "tss_start"]].reset_index(drop=True)

    logging.info(f"BED rows: total={len(tss):,} top={len(top):,} bottom={len(bottom):,}")
    return top, bottom


def compute_methylation_around_tss(allc_df, tss_df, range_start, range_end):
    """
    Vectorized approach:
      - create query table of all (chr, pos, offset_idx)
      - merge with allc_df on (chr, pos)
      - groupby offset to sum meth/cov
    """
    offsets = np.arange(range_start, range_end + 1, dtype=np.int64)
    n_pos = offsets.size

    if len(tss_df) == 0:
        # return empty summary with NaNs
        cov = np.zeros(n_pos, dtype=np.int64)
        meth = np.zeros(n_pos, dtype=np.int64)
        pct = np.full(n_pos, np.nan, dtype=float)
        counts = np.zeros(n_pos, dtype=np.int64)
        return pd.DataFrame({
            "position_relative_to_tss": offsets,
            "total_coverage": cov,
            "total_methylated": meth,
            "methylation_percentage": pct,
            "sites_with_data": counts,
        })

    # Expand: each TSS contributes n_pos query rows
    # Repeat chr and tss_start
    chr_rep = np.repeat(tss_df["chr"].to_numpy(), n_pos)
    tss_rep = np.repeat(tss_df["tss_start"].to_numpy(), n_pos)
    off_rep = np.tile(offsets, len(tss_df))
    pos_rep = tss_rep + off_rep

    query = pd.DataFrame({
        "chr": chr_rep,
        "pos": pos_rep,
        "offset": off_rep,
    })

    merged = query.merge(allc_df, on=["chr", "pos"], how="left", copy=False)

    # present rows = those with cov not null
    present = merged["cov"].notna()
    merged_present = merged.loc[present, ["offset", "meth", "cov"]]

    sums = merged_present.groupby("offset", sort=False)[["meth", "cov"]].sum()
    counts = merged_present.groupby("offset", sort=False).size().rename("sites_with_data")

    # reindex to all offsets
    sums = sums.reindex(offsets, fill_value=0)
    counts = counts.reindex(offsets, fill_value=0)

    cov = sums["cov"].to_numpy(dtype=np.int64)
    meth = sums["meth"].to_numpy(dtype=np.int64)
    pct = np.where(cov > 0, (meth / cov) * 100.0, np.nan)

    out = pd.DataFrame({
        "position_relative_to_tss": offsets,
        "total_coverage": cov,
        "total_methylated": meth,
        "methylation_percentage": pct,
        "sites_with_data": counts.to_numpy(dtype=np.int64),
    })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--allc", required=True)
    ap.add_argument("--bed", required=True, help="BED with strand column; will be split into top/bottom")
    ap.add_argument("--out_top", required=True)
    ap.add_argument("--out_bottom", required=True)
    ap.add_argument("--range_start", type=int, default=-3000)
    ap.add_argument("--range_end", type=int, default=3000)
    ap.add_argument("--strand_col", type=int, default=4, help="1-based strand column index in BED (default: 5)")
    ap.add_argument("--top_value", default="+")
    ap.add_argument("--bottom_value", default="-")
    ap.add_argument("--log_file", default=None)
    args = ap.parse_args()

    setup_logging(args.log_file)

    allc_df = read_allc_as_df(args.allc)

    bed_top, bed_bottom = read_tss_bed_split_strand(
        args.bed,
        strand_col_1based=args.strand_col,
        top_value=args.top_value,
        bottom_value=args.bottom_value,
    )

    logging.info("Computing TOP...")
    res_top = compute_methylation_around_tss(allc_df, bed_top, args.range_start, args.range_end)
    os.makedirs(os.path.dirname(args.out_top) or ".", exist_ok=True)
    res_top.to_csv(args.out_top, sep="\t", index=False, float_format="%.3f")

    logging.info("Computing BOTTOM...")
    res_bottom = compute_methylation_around_tss(allc_df, bed_bottom, args.range_start, args.range_end)
    os.makedirs(os.path.dirname(args.out_bottom) or ".", exist_ok=True)
    res_bottom.to_csv(args.out_bottom, sep="\t", index=False, float_format="%.3f")

    logging.info("Done.")


if __name__ == "__main__":
    main()
