#!/usr/bin/bash
#SBATCH -t 02:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --cluster wice
#SBATCH -A lp_big_wice_cpu
#SBATCH --partition dedicated_big_bigmem
#SBATCH -o %j.o
#SBATCH -e %j.e

dir1=
dir2=
dir3=

sh /staging/leuven/stg_00064/Xinran/repo/bioinfo_data_processing/process_fastq/concat_fastq_same_name_diff_dirs_multi.sh $dir1 $dir2 $dir3
