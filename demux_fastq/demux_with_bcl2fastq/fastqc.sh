#!/bin/bash

fastqc_path="/data/leuven/334/vsc33470/miniconda3/bin/fastqc"
project_path="/staging/leuven/stg_00064/projects/cropseq/267genes/sequencing_June2024/demultiplex_Anna_sciMet"
data_path=$project_path/fastq/
sample_list="/staging/leuven/stg_00064/projects/cropseq/267genes/sequencing_June2024/demultiplex_Anna_sciMet/sample_list.txt"
PE=true

mkdir -p $project_path/fastqc/$sample

while read sample
do
mkdir -p $project_path/fastqc/$sample

if [ $PE == true ];
then
	$fastqc_path $data_path/${sample}*_R1*.fastq.gz $data_path/${sample}*_R2*.fastq.gz -o $project_path/fastqc/$sample
else
	$fastqc_path $data_path/${sample}*_R1*.fastq.gz -o $project_path/fastqc/$sample
fi

done < $sample_list
