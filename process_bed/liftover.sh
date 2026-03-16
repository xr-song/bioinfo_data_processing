export PATH=/lustre1/project/stg_00064/Xinran/sw/UCSC_liftOver:$PATH

chain=/staging/leuven/stg_00064/Xinran/db/T2T-CHM13v2.0/liftOver_chain/hg38ToHs1.over.chain.gz
input=chromHMM.A549.ori.bed
output=chromHMM.A549.T2T.bed
unmap=liftover.unmapped.bed

liftOver $input $chain $output $unmap
