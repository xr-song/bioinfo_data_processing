#!/bin/bash
# Concate fastq from different folders (runs) with same file names

rawdatadir1="$PWD/rawdata1"
rawdatadir2="$PWD/rawdata2"
catdatadir="$PWD/rawdata"
mkdir -p "${catdatadir}"

export rawdatadir1
export rawdatadir2
export catdatadir

find ${rawdatadir1} -name '*_R1*.fastq.gz' | awk -F'_S' '{print $1}' | sort | uniq | \
parallel -j 6 ' \
	base={}; \
	id=$(basename "${base}"); \
	if [ -z "$id" ]; then exit 1; fi; \
	echo "Concatenating $id"; \
	cat ${rawdatadir1}/${id}_S*R1*.fastq.gz ${rawdatadir2}/${id}_S*R1*.fastq.gz > "${catdatadir}/${id}_R1_001.fastq.gz"; \
	cat ${rawdatadir1}/${id}_S*R2*.fastq.gz ${rawdatadir2}/${id}_S*R2*.fastq.gz > "${catdatadir}/${id}_R2_001.fastq.gz"; \
'
#rename "_001." "." *
