input_R1="${prefix}_downsampled_20_R1_001_val_1.fq.gz"
input_R2="${prefix}_downsampled_20_R2_001_val_2.fq.gz"

parallel "
    nreads={}
    millions=\$(echo \"scale=3; \$nreads/1000000\" | bc)
    output_R1=${prefix}_downsampled_\${millions}M_R1_001_val_1.fq.gz
    output_R2=${prefix}_downsampled_\${millions}M_R2_001_val_2.fq.gz
    #$seqtk sample -s100 $input_R1 \$nreads | gzip -c > \$output_R1
    $seqtk sample -s100 $input_R2 \$nreads | gzip -c > \$output_R2
    echo Downsampled to \${millions} million reads: \$output_R1 and \$output_R2
" ::: 50000 100000 500000 1000000 1500000 2000000 3000000 4000000
