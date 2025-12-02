#!/bin/bash

# rename FASTQ files by replacing _S###_L###_ with _
# Usage: ./rename_fastq.sh /path/to/fastq/folder [--dry-run]

DRY_RUN=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -*)
            echo "Unknown option $1"
            exit 1
            ;;
        *)
            FOLDER="$1"
            shift
            ;;
    esac
done

if [ -z "$FOLDER" ]; then
    echo "Error: Please provide the folder path containing FASTQ files"
    echo "Usage: $0 /path/to/fastq/folder [--dry-run]"
    echo "  --dry-run: Show what would be renamed without actually renaming"
    exit 1
fi

if [ ! -d "$FOLDER" ]; then
    echo "Error: Folder '$FOLDER' does not exist"
    exit 1
fi

cd "$FOLDER" || exit 1

echo "Checking FASTQ files in: $(pwd)"
echo "----------------------------------------"

echo "Checking for different lane numbers..."
lane_numbers=$(ls *_S[0-9]*_L[0-9]*_*.fastq.gz 2>/dev/null | sed 's/.*_L\([0-9]*\)_.*/\1/' | sort -u)

if [ -z "$lane_numbers" ]; then
    echo "No FASTQ files with _S###_L###_ pattern found"
    exit 0
fi

lane_count=$(echo "$lane_numbers" | wc -l)

echo "Found lane numbers: $lane_numbers"

if [ "$lane_count" -gt 1 ]; then
    echo "WARNING: Multiple different lane numbers detected!"
    echo "Lane numbers found: $(echo $lane_numbers | tr '\n' ' ')"
    exit 1
fi

echo "Single lane number detected: $lane_numbers"
echo "Proceeding with renaming..."

if [ "$DRY_RUN" = true ]; then
    echo "DRY RUN MODE - No files will be renamed"
fi

echo "----------------------------------------"

count=0

for file in *_S[0-9]*_L[0-9]*_*.fastq.gz; do
    if [ ! -e "$file" ]; then
        echo "No matching files found"
        break
    fi
    
    # replace _S###_L###_ with _
    new_name=$(echo "$file" | sed 's/_S[0-9]*_L[0-9]*_/_/g')
    
    # check if target file already exists
    if [ -e "$new_name" ] && [ "$file" != "$new_name" ]; then
        echo "ERROR: Target file '$new_name' already exists! Skipping '$file'"
        continue
    fi
    
    if [ "$file" != "$new_name" ]; then
        if [ "$DRY_RUN" = true ]; then
            echo "Would rename: $file -> $new_name"
        else
            echo "Renaming: $file -> $new_name"
            mv "$file" "$new_name"
        fi
        ((count++))
    else
        echo "Skipping: $file (no change needed)"
    fi
done

echo "----------------------------------------"
if [ "$DRY_RUN" = true ]; then
    echo "Would rename $count files"
else
    echo "Renamed $count files"
fi
echo "Done!"
