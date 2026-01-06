#!/usr/bin/env bash

# Usage: downsample_to_fastq.sh <input_bam> <threads> <output_dir> <sample> <reference_te_bed>
#
# 1) Calculates average and median coverage for <input_bam>, excluding regions in <reference_te_bed>.
# 2) From target coverage options [5, 10, 20, 30, 40, 50], it picks the one closest to the observed median.
#    - If the nearest target is higher than the observed coverage, no downsampling is performed.
# 3) Writes chosen coverage to <output_dir>/<sample>_coverage.txt.
# 4) Produces:
#    - coordinate-sorted downsampled BAM: <output_dir>/<sample>.bam
#    - its index: <output_dir>/<sample>.bam.bai
#    - FASTQs: <output_dir>/<sample>_1.fq, <output_dir>/<sample>_2.fq, <output_dir>/<sample>_singletons.fq.

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <input_bam> <threads> <output_dir> <sample> <reference_te_bed>"
    exit 1
fi

input_bam=$1
threads=$2
output_dir=$3
sample=$4
ref_tes=$5

mkdir -p "${output_dir}"

## Create a temporary file for the depth data in BED format:
#tmp_depth=$(mktemp -p "${output_dir}" "${sample}_depth.XXXXXX.bed")
## Convert samtools depth output into proper BED format (0-indexed start, end equals position, coverage in col4)
#samtools depth "${input_bam}" | awk '{print $1"\t"($2-1)"\t"$2"\t"$3}' > "${tmp_depth}"
#echo "Head of depth file:"
#head -n 5 "${tmp_depth}"
#
## Create a temporary file for filtered depth output
#tmp_depth_filtered=$(mktemp -p "${output_dir}" "${sample}_depth_filtered.XXXXXX.bed")
#
## (Optionally, reformat the reference BED file if needed)
#awk '{print $1"\t"$2"\t"$3}' "${ref_tes}" > "${output_dir}/${sample}_ref_temp.bed"
#
## Use bedtools to exclude positions in the filtered reference BED file.
#bedtools intersect -v -a "${tmp_depth}" -b "${output_dir}/${sample}_ref_temp.bed" > "${tmp_depth_filtered}"
#if [ $? -ne 0 ]; then
#    echo "Error: bedtools intersect failed. Check that ${ref_tes} is properly formatted (TAB-delimited, integer coordinates in cols 2 and 3)."
#    exit 1
#fi
#
## Calculate mean and median coverage using column 4 (the coverage)
#echo "Calculating average coverage for ${input_bam} excluding regions in ${ref_tes} ..."
#mean_coverage=$(awk '{sum+=$4; count++} END {if (count>0) print sum/count; else print 0}' "${tmp_depth_filtered}")
#echo "Mean coverage: ${mean_coverage}"
#
#echo "Calculating median coverage for ${input_bam} excluding regions in ${ref_tes} ..."
#coverage=$(awk '{print $4}' "${tmp_depth_filtered}" | sort -n | awk '{
#    a[NR] = $1
#}
#END {
#    if (NR % 2 == 1)
#        print a[(NR+1)/2];
#    else
#        print (a[NR/2] + a[NR/2+1]) / 2
#}')
#echo "Median coverage: ${coverage}"
#
## Clean up temporary depth file and temporary ref file
#rm -f "${tmp_depth}" "${output_dir}/${sample}_ref_temp.bed"
echo "Calculating average coverage for ${input_bam} ..."
coverage=$(samtools depth "${input_bam}" | awk '{sum+=$3; n++} END {if(n>0) printf "%.6f", sum/n; else print 0}')
echo "Average coverage: ${coverage}"

# Target coverage options (sorted ascending for ease)
coverage_levels=(5 10 20 30 40 50)
observed=${coverage}

# Find the target in the list that is numerically closest to the observed coverage.
nearest_target=0
min_diff=999999
for t in "${coverage_levels[@]}"; do
    # Calculate absolute difference between t and observed
    diff=$(awk -v a="$t" -v b="$observed" 'BEGIN {printf "%.2f", (a>=b)?a-b:b-a}')
    if (( $(echo "$diff < $min_diff" | bc -l) )); then
        min_diff=$diff
        nearest_target=$t
    fi
done
echo "Nearest target coverage is: ${nearest_target}"

# Decide: if the nearest target is greater than the observed coverage, do not downsample.
if (( $(echo "$nearest_target > $observed" | bc -l) )); then
    chosen_cov=${observed}
    echo "Observed coverage (${observed}x) is closer to a higher target (${nearest_target}x). No downsampling will be performed."
else
    chosen_cov=${nearest_target}
    echo "Downsampling to chosen coverage: ${chosen_cov}x"
fi

# Write out chosen coverage
coverage_file="${output_dir}/${sample}_coverage.txt"
echo "Using coverage: ${chosen_cov}" > "${coverage_file}"

# Temporary BAM paths in output_dir
tmp_bam="${output_dir}/tmp_${sample}.bam"
tmp_name_bam="${output_dir}/tmp_${sample}.name.bam"

# Downsample only if actual (observed) coverage is greater than the chosen target
scale_factor=1
cmp2=$(bc -l <<< "${observed} > ${chosen_cov}")
if [ "${cmp2}" -eq 1 ]; then
    scale_factor=$(bc -l <<< "${chosen_cov}/${observed}")
    echo "Coverage (${observed}x) is above the chosen target (${chosen_cov}x). Downsampling with scale factor: ${scale_factor}"
    samtools view -@ "${threads}" -b -s "${scale_factor}" "${input_bam}" > "${tmp_bam}"
else
    echo "Coverage (${observed}x) is less than or equal to chosen target (${chosen_cov}x). No downsampling performed."
    cp "${input_bam}" "${tmp_bam}"
fi

# Name-sort for FASTQ conversion
samtools sort -n -@ "${threads}" -o "${tmp_name_bam}" "${tmp_bam}"

# Generate FASTQs
fq1="${output_dir}/${sample}_1.fq"
fq2="${output_dir}/${sample}_2.fq"
fqs="${output_dir}/${sample}_singletons.fq"

echo "Converting to FASTQ ..."
samtools fastq -@ "${threads}" \
    -1 "${fq1}" \
    -2 "${fq2}" \
    -0 /dev/null \
    -s "${fqs}" \
    "${tmp_name_bam}"

# Create coordinate-sorted final BAM as output
final_bam="${output_dir}/${sample}.bam"
samtools sort -@ "${threads}" -o "${final_bam}" "${tmp_bam}"
samtools index -@ "${threads}" "${final_bam}"

# Clean up temporary BAM files
rm -f "${tmp_bam}" "${tmp_name_bam}"

echo "Done."
