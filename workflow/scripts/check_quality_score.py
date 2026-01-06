#!/usr/bin/env python3
import sys
import argparse

def detect_phred_encoding(fastq_file, sample_size=10000):
    quality_chars = []
    try:
        with open(fastq_file) as f:
            for i, line in enumerate(f):
                # Every 4th line in a FASTQ file is the quality string
                if (i + 1) % 4 == 0:
                    quality_chars.extend(line.strip())
                    if len(quality_chars) >= sample_size:
                        break
    except Exception as e:
        sys.exit(f"Error reading file {fastq_file}: {e}")
    
    if not quality_chars:
        sys.exit("No quality characters found in the file.")
    
    ascii_values = [ord(q) for q in quality_chars]
    # If any quality score is below 64, assume Phred33
    if min(ascii_values) < 64:
        return "Phred33"
    else:
        return "Phred64"

def main():
    parser = argparse.ArgumentParser(
        description="Quickly detect FASTQ quality score encoding (Phred33 or Phred64) by sampling a subset of reads."
    )
    parser.add_argument("fastq", help="Path to the FASTQ file")
    parser.add_argument(
        "--sample_size", type=int, default=10000,
        help="Number of quality characters to sample (default: 10000)"
    )
    args = parser.parse_args()

    encoding = detect_phred_encoding(args.fastq, args.sample_size)
    print(f"Detected encoding: {encoding}")

if __name__ == "__main__":
    main()
