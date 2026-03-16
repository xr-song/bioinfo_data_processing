export PATH=/staging/leuven/stg_00064/Xinran/sw/miniconda3/envs/myenv/bin:$PATH

export outdir='allc_per_cell_type'
export fasta='genome.fa'

for celltype in celltype1 celltype2 celltype3; do
	echo $celltype
	outfile="${outdir}/${celltype}.NCGN-Both.allc.tsv.gz"
	bwfile=${outdir}/${celltype}.NCGN-Both.bw
	methylpy allc-to-bigwig --allc-file $outfile --output-file $bwfile --ref-fasta $fasta --mc-type NCGN --bin-size 100 --path-to-wigToBigWig /staging/leuven/stg_00064/Xinran/sw/miniconda3/bin 
done
