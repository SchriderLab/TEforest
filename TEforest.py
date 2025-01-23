#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import textwrap
import yaml  # Ensure PyYAML is installed: pip install pyyaml
import logging

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
        required=True,
        help="Path to model pkl file"
    )
    parser.add_argument(
        "--ref_model",
        required=True,
        help="Path to reference model pkl file"
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
        "model": args.model,
        "ref_model": args.ref_model,
        #"bam_path:": args.bam_path,
        "bam_path": "/na/",
        "fq_base_path": args.fq_base_path,
        #"target_coverage": args.target_coverage,
        "target_coverage": 50

    }

    with open(config_path, 'w') as f:
        yaml.dump(config_data, f, sort_keys=False)

    logging.info(f"Config file created at: {config_path}")

    # 2) Generate output targets based on samples
    output_targets = [
        os.path.join(f"output/{sample}_TEforest_nonredundant.bed") 
        for sample in args.samples
    ]

    # Ensure the output directory exists
    output_dir = os.path.join(args.workdir, "output")
    os.makedirs(output_dir, exist_ok=True)


    # 3) Unlock working directory
    snakemake_cmd_unlock = [
        "snakemake",
    ] + output_targets + [
        "-j", str(args.threads),
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

    # 4) Prepare the Snakemake command
    snakemake_cmd = [
        "snakemake",
    ] + output_targets + [
        "-j", str(args.threads),
        "-s", snakefile_path,
        "--configfile", os.path.abspath(config_path)
    ]

    # 4) Print the command for user visibility
    logging.info("Running Snakemake command:\n" + " ".join(snakemake_cmd))

    # 5) Execute the Snakemake command from the workflow directory
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
