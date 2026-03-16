#!/usr/bin/bash
#SBATCH -t 05:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=8
#SBATCH --cluster genius
#SBATCH -A lp_labepi
#SBATCH -o %j.o
#SBATCH -e %j.e

tsv_folder=allc

export tsv_folder

check_gzip() {
    gzip -t "$1"
    if [ $? -ne 0 ]; then
        echo "Problematic file: $1"
    fi
}

export -f check_gzip

find "$tsv_folder" -name "*.tsv.gz" | parallel -j 8 check_gzip

echo "tsv.gz check done."
