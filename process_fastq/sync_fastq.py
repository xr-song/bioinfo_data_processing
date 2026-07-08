#!/usr/bin/env python3
"""
Synchronize paired-end FASTQ files (R1, R2, optional I2) to match I1 read order.

Assumptions:
  - All files contain exactly the same set of read IDs (no filtering needed).
  - I1 is used as the ordering reference.

Memory model:
  - Phase 1 : id_to_rank dict   (~50 bytes/read)
  - Phase 2+: one file at a time (~350 bytes/read for typical 150bp reads)
  Peak ≈ id_to_rank + one file, vs the original's id_to_rank × 4 files simultaneously.
"""

import gzip
import os
import sys
import shutil
import argparse


# ---------------------------------------------------------------------------
# I/O helpers (binary throughout — avoids encode/decode overhead)
# ---------------------------------------------------------------------------

def _open(path: str, mode: str = "rb"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def iter_records(path: str):
    """Yield (header, seq, plus, qual) as raw byte strings, newlines included."""
    with _open(path, "rb") as fh:
        readline = fh.readline
        while True:
            h = readline()
            if not h:
                break
            yield h, readline(), readline(), readline()


def read_id_from_header(header: bytes) -> bytes:
    """
    Extract canonical read ID:
      b'@NS500_1234/1\\n'  -> b'NS500_1234'
      b'@NS500_1234 1:N:0:ATCG\\n' -> b'NS500_1234'
    Strips '@', splits on whitespace, removes /1 or /2 pair suffix.
    """
    rid = header.split()[0][1:]          # drop '@', take first token
    if rid[-2:] in (b"/1", b"/2"):
        rid = rid[:-2]
    return rid


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def build_rank_index(i1_path: str) -> dict:
    """
    Single pass over I1: returns {read_id_bytes: int_rank}.
    This is the only structure that lives in memory for the entire run.
    """
    print(f"[1/1] Indexing I1: {i1_path}")
    index = {}
    for rank, (h, s, p, q) in enumerate(iter_records(i1_path)):
        index[read_id_from_header(h)] = rank
    print(f"      {len(index):,} reads indexed.\n")
    return index


def reorder_file(src: str, dst: str, id_to_rank: dict, label: str) -> None:
    """
    Read `src`, place each record at its rank slot, write `dst` in I1 order.
    Uses a pre-allocated list → O(N) placement, no sort required.
    """
    n = len(id_to_rank)
    print(f"[reorder] {label}: reading {src} ...")

    # Pre-allocate: index = rank position, value = raw bytes blob for 4 lines
    slots: list = [None] * n

    for h, s, p, q in iter_records(src):
        rid = read_id_from_header(h)
        rank = id_to_rank[rid]          # O(1) lookup
        slots[rank] = h + s + p + q    # join 4 lines into one bytes object

    print(f"         writing {dst} ...")
    with _open(dst, "wb") as out:
        for blob in slots:
            out.write(blob)

    print(f"         done — {n:,} reads written.\n")


def copy_file(src: str, dst: str, label: str) -> None:
    """I1 is already in the reference order — just copy it."""
    print(f"[copy  ] {label}: {src} -> {dst}")
    with _open(src, "rb") as fi, _open(dst, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    print(f"         done.\n")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def synchronize(r1: str, r2: str, i1: str, i2: str | None, prefix: str) -> None:
    ext = ".fastq.gz" if r1.endswith(".gz") else ".fastq"

    # Phase 1: build rank index from I1 (tiny memory footprint)
    id_to_rank = build_rank_index(i1)

    # Phase 2: I1 is the reference — stream-copy, no RAM needed
    copy_file(i1, f"{prefix}.I1{ext}", "I1")

    # Phase 3: reorder each remaining file independently
    #          only ONE file's data is in RAM at a time
    reorder_file(r1, f"{prefix}.R1{ext}", id_to_rank, "R1")
    reorder_file(r2, f"{prefix}.R2{ext}", id_to_rank, "R2")

    if i2:
        reorder_file(i2, f"{prefix}.I2{ext}", id_to_rank, "I2")

    print(f"All synchronized files written to: {prefix}.*")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Synchronize R1/R2/I1/I2 FASTQ files to I1 read order. "
                    "Assumes all files contain identical read sets."
    )
    parser.add_argument("-r1", "--read1",  required=True, help="R1 FASTQ file")
    parser.add_argument("-r2", "--read2",  required=True, help="R2 FASTQ file")
    parser.add_argument("-i1", "--index1", required=True, help="I1 FASTQ file (ordering reference)")
    parser.add_argument("-i2", "--index2",                help="I2 FASTQ file (optional)")
    parser.add_argument("-o",  "--output", default="synchronized", help="Output prefix")

    args = parser.parse_args()
    synchronize(args.read1, args.read2, args.index1, args.index2, args.output)
