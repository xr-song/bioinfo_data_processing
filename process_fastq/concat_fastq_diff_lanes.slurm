#!/usr/bin/bash
#SBATCH -t 00:30:00
#SBATCH --mem=2G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --cluster wice
#SBATCH -A lp_big_wice_cpu
#SBATCH --partition dedicated_big_bigmem
#SBATCH -o %j.o
#SBATCH -e %j.e

sh /staging/leuven/stg_00064/Xinran/repo/bioinfo_data_processing/process_fastq/concat_fastq_diff_lanes.sh
