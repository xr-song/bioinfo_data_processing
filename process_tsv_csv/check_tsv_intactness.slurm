#!/usr/bin/bash
#SBATCH -t 00:30:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --cluster genius
#SBATCH -A lp_labepi
#SBATCH -o %j.o
#SBATCH -e %j.e

tsv_folder=allc

for tsv_file in "$tsv_folder"/*.tsv.gz; do
  gzip -t "$tsv_file"

  if [ $? -ne 0 ]; then
    echo "Problematic file: $tsv_file"
  fi
done

echo "tsv.gz check done."

