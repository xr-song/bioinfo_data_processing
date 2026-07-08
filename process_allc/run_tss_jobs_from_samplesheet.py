#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import pandas as pd
from multiprocessing import Pool


def run_one(job):
    """
    job is a dict-like row with:
      allc_path, bed_path, out_dir, barcode, category, context, range_start, range_end, strand_col
    """
    allc_path = job["allc_path"]
    bed_path = job["bed_path"]
    out_dir = job["out_dir"]
    barcode = job.get("barcode", "sample")
    context = job.get("context", "NCG")
    category = job.get("category", "NA")
    range_start = int(job.get("range_start", -3000))
    range_end = int(job.get("range_end", 3000))
    strand_col = int(job.get("strand_col", 5))

    os.makedirs(out_dir, exist_ok=True)

    out_top = os.path.join(out_dir, f"{barcode}.{context}.{category}.mC_around_tss.top.tsv")
    out_bottom = os.path.join(out_dir, f"{barcode}.{context}.{category}.mC_around_tss.bottom.tsv")

    cmd = [
        sys.executable, job["script"],
        "--allc", allc_path,
        "--bed", bed_path,
        "--out_top", out_top,
        "--out_bottom", out_bottom,
        "--range_start", str(range_start),
        "--range_end", str(range_end),
        "--strand_col", str(strand_col),
    ]

    subprocess.run(cmd, check=True)
    return(out_top, out_bottom)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="samplesheet.csv")
    ap.add_argument("--script", required=True, help="path to get_allc_around_tss_modular.py")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    required = {"allc_path", "bed_path", "out_dir"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"samplesheet missing required columns: {sorted(missing)}")

    # attach script path to each job
    jobs = df.to_dict(orient="records")
    for j in jobs:
        j["script"] = args.script

    with Pool(args.threads) as pool:
        results = pool.map(run_one, jobs)

    # print outputs
    for out_top, out_bottom in results:
        print(f"{out_top}\t{out_bottom}")


if __name__ == "__main__":
    main()
