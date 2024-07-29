#!/bin/bash

index="Homo_sapiens.GRCh38.cdna.all.index"

for file in *_R1_trimmed.fastq.gz; do

    sample_name=$(basename "$file" _R1_trimmed.fastq.gz)

    input_r1="${sample_name}_R1_trimmed.fastq.gz"
    input_r2="${sample_name}_R2_trimmed.fastq.gz"
    output_dir="${sample_name}"
    
    log_file="${sample_name}_kallisto.log"
    
    kallisto quant -i "$index" -o "$output_dir" --verbose -t 4 "$input_r1" "$input_r2" &> "$log_file"
done

multiqc -d .
