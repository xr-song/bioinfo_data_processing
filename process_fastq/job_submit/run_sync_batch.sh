#!/usr/bin/bash
#SBATCH --job-name=sync_fastq
#SBATCH -t 10:00:00
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --cluster wice
#SBATCH -A account
#SBATCH -o logs/%j.o
#SBATCH -e logs/%j.e

export PATH=/staging/leuven/stg_00064/Xinran/sw/miniconda3/envs/myenv/bin:$PATH

set -euo pipefail
mkdir -p logs

# ── Slice: this task handles lines [START, END] of the sample list ───────────
START=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE + 1 ))
END=$(( START + BATCH_SIZE - 1 ))

echo "[Task ${SLURM_ARRAY_TASK_ID}] Processing samples ${START}–${END} from ${SAMPLE_LIST}"

# Extract this batch's I1 paths and feed to GNU parallel
sed -n "${START},${END}p" "$SAMPLE_LIST" | \
parallel -j "$JOBS_PER_BATCH" '
    i1={}
    i2=$(echo "$i1" | sed "s/_I1/_I2/")
    r1=$(echo "$i1" | sed "s/_I1/_R1/")
    r2=$(echo "$i1" | sed "s/_I1/_R2/")

    # Derive output prefix: strip input dir and _I1* suffix
    base=$(basename "$i1")
    sample="${base%%_I1*}"

    if [[ ! -f "$r1" ]]; then
        echo "WARNING: Missing R1 for $i1 — skipping"; exit 0
    fi
    if [[ ! -f "$r2" ]]; then
        echo "WARNING: Missing R2 for $i1 — skipping"; exit 0
    fi

    i2_arg=""
    [[ -f "$i2" ]] && i2_arg="-i2 $i2"

    echo "  -> $sample"
    python '"$SCRIPT"' \
        -i1 "$i1" \
        -r1 "$r1" \
        -r2 "$r2" \
        $i2_arg \
        -o '"$OUTPUT_DIR"'/"$sample"
'

echo "[Task ${SLURM_ARRAY_TASK_ID}] Finished."
