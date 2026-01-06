#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Configuration
input_file = "/nas/longleaf/home/adaigle/work/test_TEforest/coverage_debugging/combined_coverage.txt"
chunksize = 10**6          # Process file in chunks
bin_size = 10000           # Bin size in base pairs (10 kb)
chrom_order = ["2R", "2L", "3R", "3L", "X"]

# First pass: Compute maximum position for each chromosome
max_positions = {chrom: 0 for chrom in chrom_order}
print("Determining chromosome lengths...")
for chunk in pd.read_csv(input_file, sep="\t", header=None,
                         names=["chrom", "pos", "coverage"],
                         chunksize=chunksize):
    for chrom in chrom_order:
        subset = chunk.loc[chunk['chrom'] == chrom, 'pos']
        if not subset.empty:
            current_max = subset.max()
            if current_max > max_positions[chrom]:
                max_positions[chrom] = current_max

# Compute cumulative offsets for each chromosome
offsets = {}
offset = 0
for chrom in chrom_order:
    offsets[chrom] = offset
    offset += max_positions[chrom]
print("Chromosome offsets:", offsets)

# Second pass: Process file in chunks and bin the data
# We'll accumulate for each bin: sum of genome positions, sum of coverage, and count
bin_data = {}  # key: bin number, value: [sum_genome_pos, sum_coverage, count]

print("Binning coverage data...")
for chunk in pd.read_csv(input_file, sep="\t", header=None,
                         names=["chrom", "pos", "coverage"],
                         chunksize=chunksize):
    # Convert chromosome-specific position to a concatenated genome position using offsets.
    # Use .apply with axis=1 (small overhead per chunk, but chunksize limits memory use).
    chunk['genome_pos'] = chunk.apply(lambda row: row['pos'] + offsets.get(row['chrom'], 0), axis=1)
    # Determine which bin each genome position belongs to.
    chunk['bin'] = (chunk['genome_pos'] // bin_size).astype(int)
    # Group the data by bin and accumulate sums.
    for b, sub in chunk.groupby('bin'):
        s_pos = sub['genome_pos'].sum()
        s_cov = sub['coverage'].sum()
        n = len(sub)
        if b in bin_data:
            bin_data[b][0] += s_pos
            bin_data[b][1] += s_cov
            bin_data[b][2] += n
        else:
            bin_data[b] = [s_pos, s_cov, n]

# Prepare binned data for plotting
bins_sorted = sorted(bin_data.keys())
binned_positions = [bin_data[b][0] / bin_data[b][2] for b in bins_sorted]  # Average genome position per bin
binned_coverage = [bin_data[b][1] / bin_data[b][2] for b in bins_sorted]   # Average coverage per bin

# Plot the binned coverage profile
plt.figure(figsize=(12, 6))
plt.plot(binned_positions, binned_coverage, lw=0.5)
plt.xlabel("Concatenated Genomic Position (binned)")
plt.ylabel("Average Coverage")
plt.title("Genome-wide Binned Coverage Profile")
plt.tight_layout()
plt.savefig("/nas/longleaf/home/adaigle/work/test_TEforest/coverage_debugging/coverage_plot_binned.png", dpi=300)
plt.show()
