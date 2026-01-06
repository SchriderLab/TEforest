#!/bin/bash
# extract_and_fastp.sh
# This script sorts an input BAM file by read name, extracts FASTQ files,
# and then processes the FASTQ files with fastp for deduplication and adapter detection.
#
# Usage: ./extract_and_fastp.sh <input_bam> <threads>

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_bam> <threads>"
    exit 1
fi

input_bam="$1"
threads="$2"
#
## Extract sample name from input BAM file (removes path and .bam extension)
sample_name=$(basename "$input_bam" .bam)
echo "Sample name: $sample_name"
#
# Define intermediate and output file names
sorted_bam="${sample_name}.sorted.bam"
fq1="${sample_name}_1.fq"
fq2="${sample_name}_2.fq"
singletons="${sample_name}_singletons.fq"

# Define fastp output files
fastp_fq1="${sample_name}.R1.fastp.fq"
fastp_fq2="${sample_name}.R2.fastp.fq"
fastp_html="${sample_name}.fastp.html"
fastp_json="${sample_name}.fastp.json"

#echo "Sorting BAM file by read name..."
#samtools sort -n -@ "$threads" -o "$sorted_bam" "$input_bam"
#if [ $? -ne 0 ]; then
#    echo "Error: Sorting BAM file failed."
#    exit 1
#fi
#
#echo "Extracting FASTQ files from sorted BAM..."
#samtools fastq -@ "$threads" -1 "$fq1" -2 "$fq2" -s "$singletons" "$sorted_bam"
#if [ $? -ne 0 ]; then
#    echo "Error: Converting BAM to FASTQ failed."
#    exit 1
#fi
#
#echo "FASTQ extraction complete:"
#echo "Paired read file 1: $fq1"
#echo "Paired read file 2: $fq2"
#echo "Singleton reads:   $singletons"

echo "Running fastp for deduplication and adapter detection..."
fastp -w "$threads" --dedup --detect_adapter_for_pe \
    -i "$fq1" -I "$fq2" \
    -o "$fastp_fq1" -O "$fastp_fq2" \
    -h "$fastp_html" -j "$fastp_json" -z 4

if [ $? -ne 0 ]; then
    echo "Error: fastp processing failed."
    exit 1
fi

echo "fastp processing complete:"
echo "Deduplicated FASTQ file 1: $fastp_fq1"
echo "Deduplicated FASTQ file 2: $fastp_fq2"
echo "fastp HTML report: $fastp_html"
echo "fastp JSON report: $fastp_json"
