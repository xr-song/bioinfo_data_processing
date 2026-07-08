#!/usr/bin/bash
#SBATCH -t 20:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --cluster wice
#SBATCH -A account 
#SBATCH -o %j.o
#SBATCH -e %j.e

export PATH=/staging/leuven/stg_00064/Xinran/sw/miniconda3/envs/myenv/bin:$PATH
script=/staging/leuven/stg_00064/Xinran/repo/bioinfo_data_processing/process_fastq/split_fastq_by_barcode.py

for file in fastq/*S01*.fastq.gz; do
	echo $file
	python $script --fastq $file --out fastq_outs
done
