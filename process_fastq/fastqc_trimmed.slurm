#!/usr/bin/bash
#SBATCH -t 06:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --cluster wice
#SBATCH -A lp_big_wice_cpu
#SBATCH --partition dedicated_big_bigmem
#SBATCH -o %j.o
#SBATCH -e %j.e

export PATH=/staging/leuven/stg_00064/Xinran/sw/miniconda3/envs/myenv/bin:$PATH

set -euo pipefail

project_path=$(pwd)
data_path="$project_path/trim"

mkdir -p "$project_path/fastqc/trimmed"

export project_path data_path

find "$data_path" -type f -name "*_val_1.fq.gz" | \
parallel -j 8 '
    r1={}
    r2=$(echo "$r1" | sed "s/_val_1/_val_2/")
    echo "$r1"
    fastqc "$r1" "$r2" -o "$project_path/fastqc/trimmed"
'

cd "$project_path/fastqc/trimmed" && multiqc . -n multiqc_report.trimmed.html
