#!/usr/bin/bash
#SBATCH -t 00:30:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --cluster genius
#SBATCH -A lp_labepi
#SBATCH -o BAM_check.o
#SBATCH -e BAM_check.e

export PATH=/staging/leuven/stg_00064/Xinran/sw/SAMTOOLS/bin:$PATH

for bam_file in bam_sc_outs/*.bam; do
  samtools quickcheck "$bam_file"

  if [ $? -ne 0 ]; then
    echo "Problematic BAM file: $bam_file"
  fi
done

echo "BAM check done."
