#!/usr/bin/bash
#SBATCH -t 24:00:00
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --cluster wice
#SBATCH -A lp_big_wice_cpu
#SBATCH --partition dedicated_big_bigmem
#SBATCH -o %j.o
#SBATCH -e %j.e

export PATH=/staging/leuven/stg_00064/Xinran/sw/miniconda3/envs/myenv/bin:$PATH

bams=(bam_pseudobulk/*.bam)
fa=/staging/leuven/stg_00064/Xinran/db/T2T-CHM13v2.0/fasta/UCSC/genome.fa

for input_bam in "${bams[@]}"; do

tmp_vcf=$(basename "$input_bam" .pseudobulk.sorted.bam).unfiltered.vcf.gz
freebayes-parallel <(python /staging/leuven/stg_00064/Xinran/repo/sciMET/bin/fasta_generate_regions.py ${fa}.fai 100000) 20 \
    -f $fa --min-coverage 2 --min-base-quality 30 $input_bam | bgzip > $tmp_vcf

done
