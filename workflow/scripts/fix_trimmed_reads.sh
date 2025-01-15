#!/bin/bash

# Extract read IDs from R1
#seqkit -j 32 seq -n -i AKA-017_GIM-024_1.fq > R1_ids.txt

# Extract read IDs from R2
#seqkit -j 32 seq -n -i AKA-017_GIM-024_2.fq > R2_ids.txt
# Sort the read IDs
#sort R1_ids.txt > R1_ids_sorted.txt
#sort R2_ids.txt > R2_ids_sorted.txt

# Find common read IDs between R1 and R2
#comm -12 R1_ids_sorted.txt R2_ids_sorted.txt > common_ids.txt
# Add "@" symbol to the read IDs to match FASTQ headers
#awk '{print "@"$0}' common_ids.txt > common_ids_lookup.txt
# Filter R1 FASTQ based on common IDs
seqkit -j 32 grep -f common_ids.txt AKA-017_GIM-024_1.fq > R1_paired.fq

# Filter R2 FASTQ based on common IDs
seqkit -j 32 grep -f common_ids.txt AKA-017_GIM-024_2.fq > R2_paired.fq