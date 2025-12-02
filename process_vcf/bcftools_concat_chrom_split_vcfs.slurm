#!/bin/bash

usage() {
    echo "Concatenate chromosome-split vcf files for each sample"
    echo "Usage: $0 --sample_list SAMPLE_LIST --input_dir INPUT_DIR --prefix PREFIX"
    exit 1
}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample_list) sample_list="$2"; shift ;;
        --input_dir) input_dir="$2"; shift ;;
	--prefix) prefix="$2"; shift ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done

if [[ -z "$sample_list" || -z "$input_dir" ]]; then
    usage
fi

while IFS= read -r sample; do
    vcf_files=$(find "$input_dir" -name "*${sample}*.vcf.gz" ! -name "*merged*")
    if [[ -z "$vcf_files" ]]; then
        echo "No vcf found for sample: $sample"
        continue
    else
	echo "vcf files found for sample ${sample}"
	echo "$vcf_files"
    fi

    for vcf_file in $vcf_files; do
        bcftools index -t "$vcf_file"
    done

    output_file="${sample}.${prefix}.merged.vcf.gz"
    temp_vcf="${sample}.${prefix}.temp.vcf.gz"
    header_file="${sample}.${prefix}.header.txt"
    
    bcftools concat -a -Ou $vcf_files | bcftools sort -Oz -o "$temp_vcf"
    bcftools view -h "$temp_vcf" | sed "s/unknown/${sample}/" > "$header_file"
    bcftools reheader -h "$header_file" "$temp_vcf" -o "$output_file"
    rm "$temp_vcf" "$header_file" &

done < "$sample_list"

wait
