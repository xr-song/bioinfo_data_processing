#!/usr/bin/env bash
# Usage: bash submit_sync.sh
# Discovers all I1 files, groups into batches of 10, submits one SLURM job per batch.

set -euo pipefail

SCRIPT=$(realpath sync_fastq.py)
INPUT_DIR=fastq_outs
OUTPUT_DIR=fastq_synced
BATCH_SIZE=10
JOBS_PER_BATCH=5          # parallel -j N within each job
SAMPLE_LIST=sample_i1.txt  # intermediate list of I1 paths

mkdir -p "$OUTPUT_DIR"

# ── Collect all I1 paths ─────────────────────────────────────────────────────
echo "Scanning $INPUT_DIR for I1 files..."
find "$INPUT_DIR" -name "*_I1*.fastq.gz" | sort > "$SAMPLE_LIST"

TOTAL=$(wc -l < "$SAMPLE_LIST")
echo "Found $TOTAL I1 files."

if [[ $TOTAL -eq 0 ]]; then
    echo "ERROR: No I1 files found in $INPUT_DIR"; exit 1
fi

# ── Calculate array range (1-based, one task = BATCH_SIZE samples) ───────────
N_JOBS=$(( (TOTAL + BATCH_SIZE - 1) / BATCH_SIZE ))   # ceil division
echo "Submitting array of $N_JOBS jobs (batch size=$BATCH_SIZE, parallel -j $JOBS_PER_BATCH per job)"

sbatch \
    --array="1-${N_JOBS}" \
    --export=ALL,SCRIPT="$SCRIPT",SAMPLE_LIST="$SAMPLE_LIST",OUTPUT_DIR="$OUTPUT_DIR",BATCH_SIZE="$BATCH_SIZE",JOBS_PER_BATCH="$JOBS_PER_BATCH" \
    run_sync_batch.sh

echo "Done. Monitor with: squeue -u \$USER"
