#!/usr/bin/env python3

import argparse
import logging
import os
import subprocess
import sys
import textwrap
import yaml


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    parser = argparse.ArgumentParser(
        description="Launcher for the TEforest training pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Example usage:
              python TEforest_train.py \\
                --workflow_dir /path/to/TEforest/workflow/training \\
                --workdir /path/to/working_directory \\
                --threads 32 \\
                --consensusTEs /path/to/consensusTEs.fasta \\
                --ref_genome /path/to/ref_genome.fasta \\
                --ref_te_locations /path/to/ref_te_locations.bed \\
                --euchromatin /path/to/euchromatin.txt \\
                --truth_bed /path/to/truth.bed \\
                --fq_base_path /path/to/fastq_directory \\
                --samples A1 A2 A3
            """
        ),
    )

    parser.add_argument("--workflow_dir", required=True, help="Directory containing training Snakefile")
    parser.add_argument("--workdir", required=True, help="Working directory for training outputs")
    parser.add_argument("--threads", type=int, default=32, help="Threads to use")
    parser.add_argument("--samples", nargs="+", required=True, help="Sample names to process")
    parser.add_argument("--config_file", default="config.yaml", help="Config filename in workdir")

    parser.add_argument("--consensusTEs", required=True, help="Path to consensusTEs.fasta")
    parser.add_argument("--ref_genome", required=True, help="Path to reference genome fasta")
    parser.add_argument("--ref_te_locations", required=True, help="Path to reference TE BED")
    parser.add_argument("--euchromatin", required=True, help="Path to euchromatin regions")
    parser.add_argument("--truth_bed", required=True, help="Truth BED with zygosity labels")
    parser.add_argument("--fq_base_path", required=True, help="Path to FASTQ directory")
    parser.add_argument(
        "--target_coverage", type=int, default=50, help="Target coverage for downsampling"
    )

    parser.add_argument(
        "--label_mode",
        choices=["classifier", "regressor"],
        default="classifier",
        help="Label mapping mode (classifier maps 1/0.5/0 to 2/1/0).",
    )
    parser.add_argument(
        "--train_mode",
        choices=["classifier", "regressor"],
        default="classifier",
        help="Training mode for the model.",
    )
    parser.add_argument(
        "--test_size",
        type=float,
        default=0.0,
        help="Holdout fraction for train/test split (0 disables split).",
    )
    parser.add_argument("--random_state", type=int, default=42, help="Random seed")
    parser.add_argument("--n_estimators", type=int, default=500, help="LightGBM estimators")
    parser.add_argument("--n_jobs", type=int, default=8, help="LightGBM threads")
    parser.add_argument(
        "--bam_to_fvec_processes",
        type=int,
        default=None,
        help="Override processes used by bam_to_fvec (defaults to threads).",
    )
    parser.add_argument(
        "--model_out",
        default="training/teforest_classifier.pkl",
        help="Output path for trained model (relative to workdir unless absolute).",
    )

    args = parser.parse_args()

    if not os.path.isdir(args.workflow_dir):
        logging.error(f"Workflow directory '{args.workflow_dir}' does not exist.")
        sys.exit(1)

    snakefile_path = os.path.join(args.workflow_dir, "Snakefile")
    if not os.path.isfile(snakefile_path):
        logging.error(f"Snakefile not found in '{args.workflow_dir}'.")
        sys.exit(1)

    os.makedirs(args.workdir, exist_ok=True)

    model_out = args.model_out
    if not os.path.isabs(model_out):
        model_out_dir = os.path.join(args.workdir, os.path.dirname(model_out))
    else:
        model_out_dir = os.path.dirname(model_out)
    if model_out_dir:
        os.makedirs(model_out_dir, exist_ok=True)

    config_path = os.path.join(args.workdir, args.config_file)

    config_data = {
        "threads": args.threads,
        "consensusTEs": args.consensusTEs,
        "ref_genome": args.ref_genome,
        "ref_te_locations": args.ref_te_locations,
        "euchromatin": args.euchromatin,
        "truth_bed": args.truth_bed,
        "fq_base_path": args.fq_base_path,
        "target_coverage": args.target_coverage,
        "samples": args.samples,
        "label_mode": args.label_mode,
        "train_mode": args.train_mode,
        "test_size": args.test_size,
        "random_state": args.random_state,
        "n_estimators": args.n_estimators,
        "n_jobs": args.n_jobs,
        "model_out": model_out,
    }
    if args.bam_to_fvec_processes is not None:
        config_data["bam_to_fvec_processes"] = args.bam_to_fvec_processes

    with open(config_path, "w") as fh:
        yaml.dump(config_data, fh, sort_keys=False)

    logging.info(f"Config file created at: {config_path}")

    output_targets = [model_out]

    snakemake_cmd_unlock = [
        "snakemake",
    ] + output_targets + [
        "--cores",
        str(args.threads),
        "-s",
        snakefile_path,
        "--configfile",
        os.path.abspath(config_path),
        "--unlock",
    ]

    try:
        subprocess.run(snakemake_cmd_unlock, cwd=args.workdir, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Snakemake unlock failed with exit code {e.returncode}")
        sys.exit(e.returncode)

    snakemake_cmd = [
        "snakemake",
    ] + output_targets + [
        "--cores",
        str(args.threads),
        "-s",
        snakefile_path,
        "--configfile",
        os.path.abspath(config_path),
    ]

    logging.info("Running Snakemake command:\n" + " ".join(snakemake_cmd))

    try:
        subprocess.run(snakemake_cmd, cwd=args.workdir, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Snakemake failed with exit code {e.returncode}")
        sys.exit(e.returncode)


if __name__ == "__main__":
    main()
