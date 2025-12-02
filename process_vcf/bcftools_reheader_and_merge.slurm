#!/usr/bin/bash
#SBATCH -t 04:00:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --cluster wice
#SBATCH -A lp_big_wice_cpu
#SBATCH --partition dedicated_big_bigmem
#SBATCH -o %j.o
#SBATCH -e %j.e

export PATH=/staging/leuven/stg_00064/Xinran/sw/miniconda3/envs/tools/bin:$PATH

# Reheader VCF Files
for vcf in *.unfiltered.vcf.gz; do
    sample=$(echo "$vcf" | cut -d'.' -f1)
    echo "Processing $vcf -> sample name: $sample"
    bcftools reheader -s <(echo "$sample") "$vcf" -o "${sample}.unfiltered.reheadered.vcf.gz"
done

for vcf in *.reheadered.vcf.gz; do
    bcftools index "$vcf"
done

# Merge all reheadered VCF files
bcftools merge *.reheadered.vcf.gz -O z -o merged.unfiltered.vcf.gz

# Index the merged file
bcftools index merged.unfiltered.vcf.gz
