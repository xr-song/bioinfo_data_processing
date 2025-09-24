#!/bin/bash
# Concatenate (merge) fastq files from different lanes of the same sample

rawdatadir="$PWD/rawdata_all"
catdatadir="$PWD/rawdata"
mkdir -p ${catdatadir}
cd ${catdatadir}


ls ${rawdatadir}/*R1* | awk -F'_S' '{print $1}' | sort | uniq | \
parallel -j 6 ' \
    prefix={}; \
    id=$(basename "${prefix}"); \
    if [ -z "$id" ]; then exit 1; fi; \
    echo $id; \
    cat ${prefix}*R1*.fastq.gz > "${id}_R1_001.fastq.gz"; \
    cat ${prefix}*R2*.fastq.gz > "${id}_R2_001.fastq.gz"; \
'
