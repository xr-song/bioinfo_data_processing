#!/usr/bin/bash
#SBATCH -t 01:00:00
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --cluster wice
#SBATCH -A lp_big_wice_cpu
#SBATCH --partition dedicated_big_bigmem
#SBATCH -o %j.o
#SBATCH -e %j.e

dir1=part1
dir2=part2
dir3=part3
export PATH=/staging/leuven/stg_00064/Xinran/sw/miniconda3/envs/myenv/bin:$PATH
python concat_fastq_same_name_diff_dirs_multi.py $dir1 $dir2 $dir3
