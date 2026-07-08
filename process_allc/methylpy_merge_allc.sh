export PATH=/staging/leuven/stg_00064/Xinran/sw/miniconda3/envs/myenv/bin:$PATH

export allc_list_dir=''
export outdir=''

for celltype in celltype1 celltype2 celltype3; do
	echo $celltype
	allc_list="${allc_list_dir}/allc_path.${celltype}.txt"
	outfile="${outdir}/${celltype}.NCGN-Both.allc.tsv.gz"
	methylpy merge-allc --allc-files $(cat $allc_list) --output-file $outfile --num-procs 4 --compress-output True
done | parallel -j 12

# OR
for f in ${allc_list_dir}/allc_path.*.txt; do
    celltype=$(basename "$f" .txt | sed 's/allc_path\.//')
    echo $celltype
    outfile="${outdir}/${celltype}.allc.tsv.gz"
    methylpy merge-allc --allc-files $(cat "$f") --output-file "$outfile" --num-procs 4 --compress-output True
done | parallel -j 12
