#!/usr/bin/env python3

import argparse
import gzip
import os
import re
import sys
import subprocess
import textwrap
import yaml  # Ensure PyYAML is installed: pip install pyyaml
import logging

FASTQ_EXTS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def normalize_te_id(raw_id):
    return raw_id.strip().split()[0].replace("-", "_")


def _strip_fastq_ext(name):
    lower = name.lower()
    for ext in FASTQ_EXTS:
        if lower.endswith(ext):
            return name[: -len(ext)]
    return name


def _sample_remainder(name, sample):
    if name.startswith(sample):
        return name[len(sample):], 0
    pattern = re.compile(rf"(?:^|[_.-]){re.escape(sample)}(?:[_.-]|$)")
    match = pattern.search(name)
    if not match:
        return None, None
    start = match.start()
    if name[start] in "._-":
        start += 1
    return name[start + len(sample):], 1


def _infer_read_token(remainder):
    if remainder is None:
        return None
    rem = remainder.lstrip("._-")
    patterns = [
        (1, r"(?:^|[_.-])R1(?:[_.-]|$)"),
        (2, r"(?:^|[_.-])R2(?:[_.-]|$)"),
        (1, r"(?:^|[_.-])1(?:[_.-]|$)"),
        (2, r"(?:^|[_.-])2(?:[_.-]|$)"),
    ]
    for read, pattern in patterns:
        if re.search(pattern, rem):
            return read
    return None


def _list_fastq_files(base_path):
    files = []
    for name in os.listdir(base_path):
        lower = name.lower()
        if any(lower.endswith(ext) for ext in FASTQ_EXTS):
            files.append(name)
    return files


def find_fastq_path(base_path, sample, read):
    matches = []
    for name in _list_fastq_files(base_path):
        stem = _strip_fastq_ext(name)
        remainder, priority = _sample_remainder(stem, sample)
        if remainder is None:
            continue
        read_token = _infer_read_token(remainder)
        if read_token is None:
            continue
        if str(read_token) == str(read):
            matches.append((priority, name))
    if not matches:
        raise ValueError(
            f"No valid FASTQ found for sample '{sample}' read {read} in '{base_path}'."
        )
    matches.sort(key=lambda item: (item[0], item[1]))
    if len(matches) > 1:
        names = ", ".join(m[1] for m in matches)
        raise ValueError(
            f"Multiple FASTQs matched sample '{sample}' read {read}: {names}. "
            "Please merge lane files or provide one FASTQ per read."
        )
    return os.path.join(base_path, matches[0][1])


def validate_fastq_file(path):
    if not os.path.isfile(path):
        raise ValueError(f"FASTQ file not found: {path}")
    if os.path.getsize(path) == 0:
        raise ValueError(f"FASTQ file is empty: {path}")

    opener = gzip.open if path.lower().endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            lines = [handle.readline() for _ in range(4)]
    except OSError as exc:
        raise ValueError(f"FASTQ file is unreadable or invalid gzip: {path}") from exc

    if any(line == "" for line in lines):
        raise ValueError(f"FASTQ appears truncated (missing first record): {path}")

    header, seq, plus, qual = [line.rstrip("\r\n") for line in lines]
    if not header.startswith("@"):
        raise ValueError(f"FASTQ first header does not start with '@': {path}")
    if not plus.startswith("+"):
        raise ValueError(f"FASTQ third line does not start with '+': {path}")
    if len(seq) == 0 or len(qual) == 0:
        raise ValueError(f"FASTQ first record has empty sequence/quality: {path}")


def parse_consensus_ids(consensus_path):
    opener = gzip.open if consensus_path.lower().endswith(".gz") else open
    ids = set()
    with opener(consensus_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                te_id = normalize_te_id(line[1:].strip())
                if te_id:
                    ids.add(te_id)
    if not ids:
        raise ValueError(
            f"No FASTA headers were found in consensus TE file: {consensus_path}"
        )
    return ids


def parse_reference_contigs(ref_genome_path):
    opener = gzip.open if ref_genome_path.lower().endswith(".gz") else open
    contigs = set()
    with opener(ref_genome_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                header = line[1:].strip()
                if header:
                    contig = header.split()[0]
                    if contig:
                        contigs.add(contig)
    if not contigs:
        raise ValueError(
            f"No FASTA headers were found in reference genome file: {ref_genome_path}"
        )
    return contigs


def validate_reference_te_ids(ref_bed_path, consensus_path, ref_genome_path):
    consensus_ids = parse_consensus_ids(consensus_path)
    ref_contigs = parse_reference_contigs(ref_genome_path)
    bed_ids = set()
    missing_contigs = set()
    with open(ref_bed_path, "r", encoding="utf-8", errors="replace") as handle:
        for line_num, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            # Accept both tab-delimited BED and whitespace-delimited BED-like files.
            cols = stripped.split()
            if len(cols) < 7:
                raise ValueError(
                    f"Reference BED must have at least 7 columns (line {line_num}): {ref_bed_path}"
                )
            try:
                int(cols[1])
                int(cols[2])
            except ValueError as exc:
                raise ValueError(
                    f"Reference BED columns 2 and 3 must be integers (line {line_num}): {ref_bed_path}"
                ) from exc

            if cols[0] not in ref_contigs:
                missing_contigs.add(cols[0])

            te_id = normalize_te_id(cols[6])
            if not te_id:
                raise ValueError(
                    f"Reference BED column 7 TE ID is empty (line {line_num}): {ref_bed_path}"
                )
            bed_ids.add(te_id)

    if missing_contigs:
        missing_contigs = sorted(missing_contigs)
        preview = ", ".join(missing_contigs[:10])
        extra = (
            ""
            if len(missing_contigs) <= 10
            else f" ... (+{len(missing_contigs) - 10} more)"
        )
        raise ValueError(
            "Reference BED chromosomes not found in reference genome FASTA headers: "
            f"{preview}{extra}"
        )

    if not bed_ids:
        raise ValueError(f"No TE IDs were found in column 7 of BED: {ref_bed_path}")

    missing_ids = sorted(bed_ids - consensus_ids)
    if missing_ids:
        preview = ", ".join(missing_ids[:10])
        extra = "" if len(missing_ids) <= 10 else f" ... (+{len(missing_ids) - 10} more)"
        raise ValueError(
            "Reference BED TE IDs (column 7) not found in consensusTEs FASTA: "
            f"{preview}{extra}"
        )


def validate_inputs(args):
    required_files = {
        "consensusTEs": args.consensusTEs,
        "ref_genome": args.ref_genome,
        "ref_te_locations": args.ref_te_locations,
        "euchromatin": args.euchromatin,
    }
    for label, path in required_files.items():
        if not os.path.isfile(path):
            raise ValueError(f"{label} file not found: {path}")

    if not os.path.isdir(args.fq_base_path):
        raise ValueError(f"FASTQ base directory not found: {args.fq_base_path}")

    validate_reference_te_ids(
        args.ref_te_locations, args.consensusTEs, args.ref_genome
    )

    for sample in args.samples:
        fq1 = find_fastq_path(args.fq_base_path, sample, 1)
        fq2 = find_fastq_path(args.fq_base_path, sample, 2)
        validate_fastq_file(fq1)
        validate_fastq_file(fq2)
        logging.info(f"Validated FASTQs for {sample}: {fq1} | {fq2}")


def main():
    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    parser = argparse.ArgumentParser(
        description="Launcher for the TEforest Snakemake pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
            Example usage:
              python launcher.py \\
                --workflow_dir /path/to/TEforest/workflow \\
                --workdir /path/to/working_directory \\
                --threads 128 \\
                --consensusTEs /path/to/consensusTEs.fasta \\
                --ref_genome /path/to/ref_genome.fasta \\
                --ref_te_locations /path/to/ref_te_locations.bed \\
                --euchromatin /path/to/euchromatin.txt \\
                --model /path/to/svrf_classifier_all.pkl \\
                --ref_model /path/to/svrf_classifier_all_reference.pkl \\
                --fq_base_path /path/to/fastq_directory \\
                --samples A1 A2 A3 B1 B2 B3
        ''')
    )

    # Arguments for Snakemake environment
    parser.add_argument(
        "--workflow_dir",
        required=True,
        help="Directory containing the Snakefile and workflow scripts"
    )
    parser.add_argument(
        "--workdir",
        required=True,
        help="Working directory where Snakemake should execute and store outputs"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=128,
        help="Number of threads to use (defaults to 128)"
    )
    
    # Instead of --outputs, we use --samples
    parser.add_argument(
        "--samples",
        nargs='+',
        required=True,
        help="List of sample names to process"
    )
    
    # Config file name (will be created inside workdir)
    parser.add_argument(
        "--config_file",
        default="config.yaml",
        help="Name of the config file to be created inside the working directory"
    )

    # Arguments for pipeline configuration
    parser.add_argument(
        "--consensusTEs",
        required=True,
        help="Path to consensusTEs.fasta"
    )
    parser.add_argument(
        "--ref_genome",
        required=True,
        help="Path to reference genome fasta"
    )
    parser.add_argument(
        "--ref_te_locations",
        required=True,
        help="Path to reference TE locations bed file"
    )
    parser.add_argument(
        "--euchromatin",
        required=True,
        help="Path to euchromatin region file"
    )
    parser.add_argument(
        "--model",
        required=False,
        help="Path to model pkl file (optional; if omitted, TEforest auto-selects based on coverage)"
    )
    parser.add_argument(
        "--ref_model",
        required=False,
        help="Path to reference model pkl file (optional; if omitted, TEforest auto-selects based on coverage)"
    )
    #parser.add_argument(
    #    "--bam_path",
    #    required=True,
    #    help="Path to directory containing BAM files"
    #)
    parser.add_argument(
        "--fq_base_path",
        required=True,
        help="Path to directory containing fastq files"
    )
    parser.add_argument(
        "--cleanup_intermediates",
        action="store_true",
        help="Delete large intermediate files (fastp/, aligned/, downsampled/, candidate_regions_data/) after use."
    )
    parser.add_argument(
        "--disable_reference_detection",
        action="store_true",
        help="Disable reference TE detection (skip reference feature vectors and reference model prediction)."
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Run Snakemake in dry-run mode to validate DAG and inputs without executing jobs."
    )
    #parser.add_argument(
    #    "--target_coverage",
    #    type=int,
    #    required=True,
    #    help="Target coverage for the run"
    #)

    args = parser.parse_args()

    # Validate workflow directory
    if not os.path.isdir(args.workflow_dir):
        logging.error(f"Workflow directory '{args.workflow_dir}' does not exist.")
        sys.exit(1)
    
    # Path to the Snakefile
    snakefile_path = os.path.join(args.workflow_dir, "Snakefile")
    if not os.path.isfile(snakefile_path):
        logging.error(f"Snakefile not found in workflow directory '{args.workflow_dir}'. Expected at '{snakefile_path}'.")
        sys.exit(1)

    try:
        validate_inputs(args)
    except ValueError as exc:
        logging.error(f"Input validation failed: {exc}")
        sys.exit(1)

    # Ensure working directory exists
    if not os.path.isdir(args.workdir):
        logging.info(f"Working directory '{args.workdir}' does not exist. Creating it.")
        os.makedirs(args.workdir, exist_ok=True)

    # Path to the config file inside the working directory
    config_path = os.path.join(args.workdir, args.config_file)

    # 1) Write the config file (YAML) dynamically, based on user arguments
    config_data = {
        "threads": args.threads,
        "consensusTEs": args.consensusTEs,
        "ref_genome": args.ref_genome,
        "ref_te_locations": args.ref_te_locations,
        "euchromatin": args.euchromatin,
        "model": args.model if args.model else "auto",
        "ref_model": args.ref_model if args.ref_model else "auto",
        #"bam_path:": args.bam_path,
        "bam_path": "/na/",
        "fq_base_path": args.fq_base_path,
        #"target_coverage": args.target_coverage,
        "target_coverage": 50,
        "cleanup_intermediates": args.cleanup_intermediates,
        "disable_reference_detection": args.disable_reference_detection,

    }

    with open(config_path, 'w') as f:
        yaml.dump(config_data, f, sort_keys=False)

    logging.info(f"Config file created at: {config_path}")

    # 2) Generate output targets based on samples
    output_targets = []
    for sample in args.samples:
        output_targets.append(f"output/{sample}_TEforest_nonredundant.bed")
        output_targets.append(f"output/{sample}_TEforest_bps_nonredundant.vcf")

    # Ensure the output directory exists
    output_dir = os.path.join(args.workdir, "output")
    os.makedirs(output_dir, exist_ok=True)


    # 3) Prepare the Snakemake command
    snakemake_cmd = [
        "snakemake",
    ] + output_targets + [
        "--cores", str(args.threads),
        "-s", snakefile_path,
        "--configfile", os.path.abspath(config_path)
    ]

    if args.dry_run:
        snakemake_cmd.append("--dry-run")
        logging.info("Running Snakemake dry-run command:\n" + " ".join(snakemake_cmd))
        try:
            subprocess.run(
                snakemake_cmd,
                cwd=args.workdir,
                check=True
            )
        except subprocess.CalledProcessError as e:
            logging.error(f"Snakemake dry-run failed with exit code {e.returncode}")
            sys.exit(e.returncode)
        return

    # 4) Unlock working directory
    snakemake_cmd_unlock = [
        "snakemake",
    ] + output_targets + [
        "--cores", str(args.threads),
        "-s", snakefile_path,
        "--configfile", os.path.abspath(config_path),
        "--unlock"
    ]

    try:
        subprocess.run(
            snakemake_cmd_unlock,
            cwd=args.workdir,      # Set workflow_dir as current directory
            check=True
        )
    except subprocess.CalledProcessError as e:
        logging.error(f"Snakemake failed with exit code {e.returncode}")
        sys.exit(e.returncode)

    # 5) Print the command for user visibility
    logging.info("Running Snakemake command:\n" + " ".join(snakemake_cmd))

    # 6) Execute the Snakemake command from the workflow directory
    try:
        subprocess.run(
            snakemake_cmd,
            cwd=args.workdir,      # Set workdir as current directory
            check=True
        )
    except subprocess.CalledProcessError as e:
        logging.error(f"Snakemake failed with exit code {e.returncode}")
        sys.exit(e.returncode)

if __name__ == "__main__":
    main()
