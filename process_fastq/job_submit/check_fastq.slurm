#!/usr/bin/bash
#SBATCH -t 08:00:00
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --cluster wice
#SBATCH -A account 
#SBATCH -o %j.o
#SBATCH -e %j.e

export PATH=/staging/leuven/stg_00064/Xinran/sw/miniconda3/envs/myenv/bin:$PATH

export script=/staging/leuven/stg_00064/Xinran/repo/bioinfo_data_processing/process_fastq/check_fastq_unmatched_read_names.py

input_dir=fastq_outs

find "$input_dir" -name "*_I1*.fastq.gz" | \
    parallel -j 4 '
        i1={}
        i2=$(echo {} | sed "s/_I1/_I2/")
        
        if [ -f "$i2" ]; then
            echo "Checking: $i1 vs $i2"
            python $script "$i1" "$i2"
        else
            echo "WARNING: No matching I2 file for $i1"
        fi
    '
