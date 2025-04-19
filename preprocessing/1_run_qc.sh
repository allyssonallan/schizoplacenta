#!/bin/bash

fastqc *.gz -t 4

if [[ ! "$filename" =~ ^[A-Za-z0-9_-]+_R[12]_001\.fastq\.gz$ ]]; then
    echo "Invalid filename: $filename" >&2
    exit 1
fi

# Loop through all PE files
for file in *_R1_001.fastq.gz; do

    sample_name=$(basename "$file" _R1_001.fastq.gz)

    input_r1="${sample_name}_R1_001.fastq.gz"
    input_r2="${sample_name}_R2_001.fastq.gz"
    output_r1="${sample_name}_R1_trimmed.fastq.gz"
    output_r2="${sample_name}_R2_trimmed.fastq.gz"

    log_file="${sample_name}.log"

    fastp -i "$input_r1" -I "$input_r2" -o "$output_r1" -O "$output_r2" &> "$log_file"
done

fastqc *_trimmed.fastq.gz -t 4
