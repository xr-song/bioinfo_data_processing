#!/bin/bash
# Extract specific tsv column based on column name
# Usage: ./select_cols.sh input.tsv col1 [col2]

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 input.tsv col1 [col2]"
    exit 1
fi

input_file="$1"
col1="$2"
col2="$3"  # may be empty if only 1 column

awk -F'\t' -v c1="$col1" -v c2="$col2" '
NR==1 {
    col1_idx = col2_idx = -1
    for (i=1; i<=NF; i++) {
        if ($i==c1) col1_idx=i
        if (c2!="" && $i==c2) col2_idx=i
    }
    if (col1_idx == -1) {
        print "Error: Column \"" c1 "\" not found in header" > "/dev/stderr"
        exit 1
    }
    if (c2!="" && col2_idx == -1) {
        print "Error: Column \"" c2 "\" not found in header" > "/dev/stderr"
        exit 1
    }
    if (c2=="") {
        print $col1_idx
    } else {
        print $col1_idx, $col2_idx
    }
    next
}
{
    if (c2=="") {
        print $col1_idx
    } else {
        print $col1_idx, $col2_idx
    }
}
' OFS='\t' "$input_file"

