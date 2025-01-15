#!/bin/bash

# for the trimmed reads experiment, downsample a bam file to 30x covg

# Check for the correct number of arguments
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <input_bam_file> <threads> <output_dir> <sample> <target_coverage>"
    exit 1
fi

input_bam=$1
threads=$2
output_dir=$3
sample=$4
target_coverage=$5

# Create output directory if it doesn't exist
mkdir -p $output_dir

# Step 1: Calculate average coverage of the input BAM file
echo "Calculating average coverage for $input_bam"
coverage=$(samtools depth $input_bam | awk '{sum+=$3; count++} END {print sum/count}')

# Step 2: Calculate scale factor for downsampling
scale_factor=$(echo $coverage | awk -v target=$target_coverage '{print target/$1}')
echo "Scale factor: $scale_factor"

# Step 3: Downsample the BAM file
downsampled_bam="${output_dir}/${sample}_downsampled.bam"
echo "Downsampling $input_bam to $target_coverage x coverage"
samtools view -@ $threads -s $scale_factor -b -o $downsampled_bam $input_bam

# Step 4: Sort the downsampled BAM file
sorted_downsampled_bam="${output_dir}/${sample}_downsampled_sorted.bam"
echo "Sorting the downsampled BAM file"
samtools sort -n -@ $threads -o $sorted_downsampled_bam $downsampled_bam

# Step 5: Extract FASTQ reads from the sorted downsampled BAM file
downsampled_fastq_1="${output_dir}/${sample}_downsampled_R1.fastq"
downsampled_fastq_2="${output_dir}/${sample}_downsampled_R2.fastq"

echo "Extracting FASTQ reads from $sorted_downsampled_bam"
samtools fastq -@ $threads -1 $downsampled_fastq_1 -2 $downsampled_fastq_2 -0 /dev/null -s /dev/null -n $sorted_downsampled_bam

echo "FASTQ files generated:"
echo " - $downsampled_fastq_1"
echo " - $downsampled_fastq_2"