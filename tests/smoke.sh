#!/usr/bin/env bash
set -euo pipefail

echo "Python:"
python --version

echo "Snakemake:"
snakemake --version

echo "R:"
R --version | head -n 2

echo "Key tools:"
bwa-mem2 version || true
samtools --version | head -n 2
bedtools --version
fastp --version || true
seqkit version || true
