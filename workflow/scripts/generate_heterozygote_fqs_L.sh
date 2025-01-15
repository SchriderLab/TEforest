#!/bin/bash
#goal. Take two bam files as input, 
#downsample both to 30x coverage,
# make chromosome 2L and 3L 15x covg, 
# combine both genomes (and first half of X)
# others just keep genome 1.


# Check for the correct number of arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <input_bam_file1> <input_bam_file2> <threads> <output_dir> <sample1> <sample2> <target_coverage>"
    exit 1
fi

#module load samtools
# Input BAM files
input_bam1=$1
input_bam2=$2
threads=$3
output_dir=$4
sample1=$5
sample2=$6
target_coverage=$7

if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

cd $output_dir
# Calculate average coverage by chromosome for the first BAM file
echo "Calculating average coverage by chromosome for $input_bam1"
coverage1=$(samtools depth $input_bam1 | awk '{sum[$1]+=$3; count[$1]++} END {for(chr in sum) print chr, sum[chr]/count[chr]}')
echo "Calculating average coverage for $input_bam1"
coverage1=$(samtools depth $input_bam1 | awk '{sum+=$3; count++} END {print sum/count}')
#coverage1=49.4922
echo "Average coverage $coverage1"

echo "Calculating average coverage for $input_bam2"
coverage2=$(samtools depth $input_bam2 | awk '{sum+=$3; count++} END {print sum/count}')
#coverage1=49.4922
echo "Average coverage $coverage2"
find_min() {
    local coverage1="$1"
    local coverage2="$2"
    # Compare the coverages and assign the smaller value
    if (( $(echo "$coverage1 < $coverage2" | bc -l) )); then
        echo "$coverage1"
    else
        echo "$coverage2"
    fi
}
# Update target_coverage if it's greater than the minimum true coverage
min_true_coverage=$(find_min "$coverage1" "$coverage2")
target_coverage=$(find_min "$min_true_coverage" "$target_coverage")


# Calculate target_coverage2 as half of the updated target_coverage
# Usage: divide_by_two <value>
divide_by_two() {
    local input="$1"
    echo "scale=10; $input / 2" | bc
}

target_coverage2=$(divide_by_two "$target_coverage")

# Print the results
echo "Updated target_coverage: $target_coverage"
echo "target_coverage2: $target_coverage2"

# Extract the scale factor for downsampling
scale_factor1=$(echo $coverage1 | awk -v target=$target_coverage '{print target/$1}')
echo "Scale factor $scale_factor1"
scale_factor1_50=$(echo $coverage1 | awk -v target=$target_coverage2 '{print target/$1}')
echo "Scale factor $scale_factor1_50"


 #Extract the scale factor for downsampling
scale_factor2=$(echo $coverage2 | awk -v target=$target_coverage '{print target/$1}')
echo "Scale factor $scale_factor2"
scale_factor2_50=$(echo $coverage2 | awk -v target=$target_coverage2 '{print target/$1}')
echo "Scale factor $scale_factor2_50"

# Downsample the first BAM file using the calculated scale factor
output_bam1="30x_downsampled_genome1.bam"
output_bam1_50="15x_downsampled_genome1.bam"

echo "Downsampling $input_bam1 to $target_coverage x coverage"
samtools view -@ $threads -s $scale_factor1 -b -o $output_bam1 $input_bam1
samtools view -@ $threads -s $scale_factor1_50 -b -o $output_bam1_50 $input_bam1

# Downsample the second BAM file using the calculated scale factor
output_bam2="30x_downsampled_genome2.bam"
output_bam2_50="15x_downsampled_genome2.bam"

echo "Downsampling $input_bam2 to $target_coverage x coverage"
samtools view -@ $threads -s $scale_factor2 -b -o $output_bam2 $input_bam2
samtools view -@ $threads -s $scale_factor2_50 -b -o $output_bam2_50 $input_bam2

## Index the output BAM files
samtools index -@ $threads $output_bam1
samtools index -@ $threads $output_bam2
samtools index -@ $threads $output_bam1_50
samtools index -@ $threads $output_bam2_50
echo "Downsampling complete."

# Extract 1/2 reads for chromosome 2R
samtools view -@ $threads -b $output_bam1_50 2R -o 2R_reads_genome1.bam
samtools view -@ $threads -b $output_bam2_50 2R -o 2R_reads_genome2.bam
samtools view -@ $threads -b $output_bam1_50 3R -o 3R_reads_genome1.bam
samtools view -@ $threads -b $output_bam2_50 3R -o 3R_reads_genome2.bam

samtools view -@ $threads -b $output_bam1_50 2L -o 2L_reads_genome1.bam
samtools view -@ $threads -b $output_bam2_50 2L -o 2L_reads_genome2.bam
samtools view -@ $threads -b $output_bam1_50 3L -o 3L_reads_genome1.bam
samtools view -@ $threads -b $output_bam2_50 3L -o 3L_reads_genome2.bam
samtools view -@ $threads -b $output_bam1_50 X -o X_reads_genome1.bam
samtools view -@ $threads -b $output_bam2_50 X -o X_reads_genome2.bam
#
## Combine BAM files from genome 1 and genome 2 for each matching chromosome
#old version
#samtools cat -@ $threads -o combined_2R.bam 2R_reads_genome1.bam 2R_reads_genome2.bam
#samtools cat -@ $threads -o combined_3R.bam 3R_reads_genome1.bam 3R_reads_genome2.bam
#
samtools cat -o combined_2R3R_genome1.bam 2R_reads_genome1.bam 3R_reads_genome1.bam 2L_reads_genome1.bam 3L_reads_genome1.bam X_reads_genome1.bam
#find reads specific to 2R and 3R of genome2
#this is different than genome1 because we already have these reads for genome1
samtools cat -o combined_2R3R_genome2.bam 2R_reads_genome2.bam 3R_reads_genome2.bam 2L_reads_genome2.bam 3L_reads_genome2.bam X_reads_genome2.bam
samtools view -@ ${threads} -h combined_2R3R_genome2.bam | awk '{if ($1 ~ /^@/) next; print $1}' > genome2_readnames_2R3R.txt
#grab those reads from the bam file (gets TE reads that map to other chroms)
samtools view -H $input_bam2 > header.sam
samtools view $input_bam2 | grep -wFf genome2_readnames_2R3R.txt > genome2_2R3R.sam
cat header.sam genome2_2R3R.sam | samtools view -bS - > genome2_2R3R.bam
output_bam2R3R="30x_downsampled_genome1_without2R3R.bam"
#samtools view -@ $threads -o $output_bam2R3R -b $output_bam1 2R 3R
#samtools view -@ $threads -o $output_bam2R3R -b $output_bam1 | grep -v -E "2R|3R"
#samtools view -h $output_bam1 | awk '{if($3 != "2R" && $3 != "3R"){print $0}}' | samtools view -Sb - > $output_bam2R3R
samtools idxstats $output_bam1 | cut -f 1 | grep -v -E "2R|3R|2L|3L|X" | xargs samtools view -b $output_bam1 > $output_bam2R3R
#
final_output="final_product.bam"
echo $output_bam2R3R
##samtools cat -b combined_2R.bam combined_3R.bam > $final_output
#samtools cat -@ $threads -o $final_output combined_2R.bam combined_3R.bam 30x_downsampled_genome1_without2R3R.bam
samtools cat -o $final_output combined_2R3R_genome1.bam genome2_2R3R.bam 30x_downsampled_genome1_without2R3R.bam

final_output_sort="final_product.sort.bam"
final_output_name_sort="final_product.name.sort.bam"

#sorting by read name so I can extract fq's in proper format. 
samtools sort -@ $threads -o $final_output_sort $final_output
samtools index $final_output_sort
samtools sort -n -@ $threads -o $final_output_name_sort $final_output

#final_dir_name=$(basename "$output_dir")
sample_name="${sample1}_${sample2}"
echo $final_output_sort

samtools fastq -@ $threads -1 "${sample_name}_1.fq" -2 "${sample_name}_2.fq" -s singletons.fq $final_output_name_sort

