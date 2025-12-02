#!/bin/bash

rawdatadir="$1"
catdatadir="$2"

mkdir -p "$catdatadir"

export rawdatadir
export catdatadir

ls "${rawdatadir}"/*R1* \
  | sed -E 's/^(.*)_S[0-9]+.*/\1/' \
  | sort | uniq \
  | parallel -j 6 --env rawdatadir --env catdatadir "
      prefix={};
      id=\$(basename \"\$prefix\");

      echo \"\$id\"

      cat \"\${prefix}\"*R1*.fastq.gz > \"\$catdatadir/\${id}_R1_001.fastq.gz\"
      cat \"\${prefix}\"*R2*.fastq.gz > \"\$catdatadir/\${id}_R2_001.fastq.gz\"
    "

