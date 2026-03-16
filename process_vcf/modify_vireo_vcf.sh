# 1. Create a "Strict" Header file (Vireo usually misses these)
cat <<EOF > fix_header.txt
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihoods">
##INFO=<ID=AD,Number=R,Type=Integer,Description="Total Allelic depths">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=10>
##contig=<ID=11>
##contig=<ID=12>
##contig=<ID=13>
##contig=<ID=14>
##contig=<ID=15>
##contig=<ID=16>
##contig=<ID=17>
##contig=<ID=18>
##contig=<ID=19>
##contig=<ID=20>
##contig=<ID=21>
##contig=<ID=22>
##contig=<ID=X>
##contig=<ID=Y>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
EOF

# 2. Extract and Sort the data lines (skipping existing broken headers)
# We use 'sort -V' for natural version sort (1, 2, 10 instead of 1, 10, 2)
zgrep -v "^#" GT_donors.vireo.vcf.gz | sort -V -k1,1 -k2,2n > sorted_body.txt

# 3. Merge Header and Data
cat fix_header.txt sorted_body.txt > GT_donors.final.vcf

# 4. Compress using BGZF (the format IGV requires)
bgzip -f GT_donors.final.vcf

# 5. Index the file
tabix -p vcf GT_donors.final.vcf.gz

# 6. Cleanup temporary files
rm fix_header.txt sorted_body.txt
