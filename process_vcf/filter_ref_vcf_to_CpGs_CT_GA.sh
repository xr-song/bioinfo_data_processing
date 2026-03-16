#!/usr/bin/bash

# script for filtering reference vcf file (human common SNPs)

export PATH=/staging/leuven/stg_00064/Xinran/sw/miniconda3/envs/myenv/bin:$PATH

ref_vcf=/staging/leuven/stg_00064/Xinran/db/hg38/SNP/1000G_phase3_cellsnp/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.chr.sorted.vcf.gz
bedfile=/staging/leuven/stg_00064/Xinran/db/hg38/CpG/Human_hg38.CpGsites.0based.sorted.bed
output_vcf=$(echo "$ref_vcf" | sed -E 's/\.vcf(\.gz)?$/.CpGs_only.vcf/')
output_vcf2=$(echo "$ref_vcf" | sed -E 's/\.vcf(\.gz)?$/.CpGs_only.CT_GA.vcf.gz/')
output_bed=$(echo "$ref_vcf" | sed -E 's/\.vcf(\.gz)?$/.CpGs_only.CT_GA.bed/')

# subset the vcf file to the ranges in the bed file
bcftools view -R $bedfile $ref_vcf > $output_vcf

# subset to only C->T and G->A conversions
bcftools view -i 'REF="C" && ALT="T" || REF="G" && ALT="A"' $output_vcf | bgzip -c > $output_vcf2
# index
tabix -p vcf $output_vcf2

# convert to bed file
bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' $output_vcf2 > $output_bed
