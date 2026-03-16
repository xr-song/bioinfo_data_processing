#!/usr/bin/bash

#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH -A lp_labepi
#SBATCH --cluster wice
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=qian.yu@kuleuven.be


cd $SLURM_SUBMIT_DIR
#sh demultiplex.sh
sh fastqc.sh
