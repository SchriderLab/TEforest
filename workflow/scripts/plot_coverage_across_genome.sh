#!/bin/bash
# coverage_analysis_parallel.sh
# This script processes chromosomes 2R, 2L, 3R, 3L, and X in parallel,
# computes per-base coverage using bedtools genomecov,
# concatenates the coverage data, and then creates a genome-wide coverage plot.

# Set input BAM file and output directory
bam_file="/nas/longleaf/home/adaigle/work/test_TEforest/basenorm_feats_50X_smaller_regions/aligned_het/A2_A3.bam"
output_dir="/nas/longleaf/home/adaigle/work/test_TEforest/coverage_debugging"
mkdir -p "$output_dir"

# Define the chromosomes to process in the desired order
chroms=("2R" "2L" "3R" "3L" "X")
combined_file="$output_dir/combined_coverage.txt"
rm -f "$combined_file"  # Remove any existing combined file

# Loop over each chromosome in parallel
for chrom in "${chroms[@]}"; do
    (
      echo "Processing chromosome $chrom ..."
      # Use multiple threads for samtools view (if your system supports it)
      samtools view -@ 4 -b "$bam_file" "$chrom" > "$output_dir/${chrom}.bam"
      
      # Generate per-base coverage with bedtools genomecov.
      bedtools genomecov -d -ibam "$output_dir/${chrom}.bam" > "$output_dir/${chrom}_coverage.txt"
      
      # Append the chromosome's coverage data to the combined file.
      cat "$output_dir/${chrom}_coverage.txt" >> "$combined_file"
      echo "Finished processing $chrom"
    ) &
done

# Wait for all background processes to complete
wait
echo "Combined coverage data saved to $combined_file"

# Use Python to plot genome-wide coverage.
python3 << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt

# Read the combined coverage file (columns: chrom, pos, coverage)
df = pd.read_csv(r'./coverage_analysis/combined_coverage.txt', sep="\t", header=None, names=["chrom", "pos", "coverage"])
df['pos'] = df['pos'].astype(int)
df['coverage'] = df['coverage'].astype(float)

# Define the desired chromosome order for concatenation.
chrom_order = ["2R", "2L", "3R", "3L", "X"]

# Determine maximum positions as a proxy for chromosome lengths.
chrom_lengths = {chrom: df.loc[df['chrom'] == chrom, 'pos'].max() for chrom in chrom_order}

# Compute cumulative offset for each chromosome.
offset = 0
offsets = {}
for chrom in chrom_order:
    offsets[chrom] = offset
    offset += chrom_lengths.get(chrom, 0)

# Create a new column for the genome-wide concatenated position.
df['genome_pos'] = df.apply(lambda row: row['pos'] + offsets.get(row['chrom'], 0), axis=1)

# Plot the genome-wide coverage profile.
plt.figure(figsize=(12, 6))
plt.plot(df['genome_pos'], df['coverage'], lw=0.5)
plt.xlabel("Concatenated Genomic Position")
plt.ylabel("Coverage")
plt.title("Genome-wide Coverage Profile (Chromosomes 2R, 2L, 3R, 3L, X)")
plt.tight_layout()
plt.savefig("./coverage_analysis/coverage_plot.png", dpi=300)
plt.show()
EOF

echo "Plot saved to $output_dir/coverage_plot.png"
