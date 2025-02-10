#!/usr/bin/env bash

# Usage: downsample_to_fastq.sh <input_bam> <threads> <output_dir> <sample>
#
# 1) Calculates average coverage for <input_bam>.
# 2) Picks next-lower coverage from [5, 10, 20, 30, 40, 50].
#    - If coverage < 5 => no downsampling.
# 3) Writes chosen coverage to <output_dir>/<sample>_coverage.txt.
# 4) Produces:
#    - coordinate-sorted downsampled BAM: <output_dir>/<sample>.bam
#    - its index: <output_dir>/<sample>.bam.bai
#    - FASTQs: <output_dir>/<sample>_1.fq, <output_dir>/<sample>_2.fq, <output_dir>/<sample>_singletons.fq.

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_bam> <threads> <output_dir> <sample>"
    exit 1
fi

input_bam=$1
threads=$2
output_dir=$3
sample=$4

# Ensure output directory exists
mkdir -p "${output_dir}"

echo "Calculating average coverage for ${input_bam} ..."
coverage=$(samtools depth "${input_bam}" | awk '{sum+=$3; count++} END {if (count>0) print sum/count; else print 0}')
echo "Average coverage detected: ${coverage}"

# Coverage thresholds in descending order
coverage_levels=(50 40 30 20 10 5)

# Determine chosen coverage
chosen_cov=0
for lvl in "${coverage_levels[@]}"; do
    cmp=$(bc -l <<< "${coverage} >= ${lvl}")
    if [ "${cmp}" -eq 1 ]; then
        chosen_cov=${lvl}
        break
    fi
done

# If coverage <5 => no downsampling
if [ "${chosen_cov}" -eq 0 ]; then
    chosen_cov=${coverage}
fi

# Write out chosen coverage
coverage_file="${output_dir}/${sample}_coverage.txt"
echo "Using coverage: ${chosen_cov}" > "${coverage_file}"

# Temporary BAM paths in output_dir
tmp_bam="${output_dir}/tmp_${sample}.bam"
tmp_name_bam="${output_dir}/tmp_${sample}.name.bam"

# Downsample only if actual coverage > chosen_cov
scale_factor=1
cmp2=$(bc -l <<< "${coverage} > ${chosen_cov}")
if [ "${cmp2}" -eq 1 ]; then
    scale_factor=$(bc -l <<< "${chosen_cov}/${coverage}")
    echo "Coverage above ${chosen_cov}x. Downsampling with scale factor: ${scale_factor}"
    samtools view -@ "${threads}" -b -s "${scale_factor}" "${input_bam}" > "${tmp_bam}"
else
    echo "Coverage <= ${chosen_cov}x. No downsampling performed."
    cp "${input_bam}" "${tmp_bam}"
fi

# Name-sort for FASTQ
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

# Create coordinate-sorted BAM as final output
final_bam="${output_dir}/${sample}.bam"
samtools sort -@ "${threads}" -o "${final_bam}" "${tmp_bam}"
samtools index -@ "${threads}" "${final_bam}"

# Clean up
rm -f "${tmp_bam}" "${tmp_name_bam}"

echo "Done."