# TEforest

This repository is designed for detecting transposable element (TE) insertions using short-read sequencing data. It enables both the assessment of annotated TE presence/absence in the reference genome and the detection and genotyping of novel insertions not represented in the reference. The repository includes a Snakemake pipeline used for training, deploying, and benchmarking TE insertion detection models based on LightGBM classifiers.

## Overview

- **Snakemake Pipeline**  
  - The main pipeline is controlled by `workflow/Snakefile`.  
  - It orchestrates rules for:
    - Data preprocessing
    - Model training
    - Model deployment (prediction)
    - Benchmarking
  - The new user-facing training pipeline lives in `workflow/training` and is launched via `TEforest_train.py`. Some legacy training/benchmarking scripts still contain hard-coded paths.

- **Scripts**  
  - Supporting scripts for training and model deployment are located in `workflow/scripts`.

- **Pre-Trained Models**  
  - Trained Random Forest models are stored in `workflow/models`.

- **Outputs**  
  - Final output files (predictions) will be generated in your chosen work directory under the path:
    ```bash
    outputs/{sample}_TEforest_bps_nonredundant.bed
    outputs/{sample}_TEforest_bps_nonredundant.vcf
    ```
  - Candidate regions prior to precise breakpoint detection are also avaliable: 
    
    ```bash
    outputs/{sample}_TEforest_nonredundant.bed
    ```

## Environment Setup

Two supported options are available:

### Conda

Create and activate the environment with the provided YAML:

```bash
conda env create -f TEforest.yml
conda activate TEforest
```

### Docker (GitHub container)

We publish a container image via GitHub Container Registry (GHCR). Use the image from the repository to run TEforest without managing local dependencies. See the repository’s GHCR workflow for the current image name and tags.

## How to Launch

A Python script named **`TEforest.py`** is provided for launching the pipeline. This script will run the trained models on one or more genomes of your choice.

> **Note**: Ongoing usability developments may change the command-line arguments or environment requirements in the future.

### Basic Usage

```bash
python TEforest.py \
    --workflow_dir <path/to/workflow> \
    --workdir <path/to/desired_workdir> \
    --threads 16 \
    --consensusTEs <path/to/consensusTEs.fasta> \
    --ref_genome <path/to/reference_genome.fasta> \
    --ref_te_locations <path/to/te_locations.bed> \
    --euchromatin <path/to/euchromatin.bed> \
    --model <path/to/pretrained_model.pkl> \
    --ref_model <path/to/reference_model.pkl> \
    --fq_base_path <path/to/fastq/files> \
    --cleanup_intermediates \
    --disable_reference_detection \
    --dry_run \
    --samples A1 A2 A3
```

- **`--workflow_dir`**: Directory containing the `Snakefile` (`workflow/Snakefile`).
- **`--workdir`**: Directory to store outputs and logs.
- **`--threads`**: Number of CPU threads to use. 16 per sample is recommended.
- **`--consensusTEs`, `--ref_genome`, `--ref_te_locations`, `--euchromatin`**: Input reference files for TE detection. All calls outside of the regions denoted in euchromatin will be filtered. Example files used for Drosophila melanogaster are located in example_files/. Be aware the BWA-mem2 will treat IUPAC bases as missing, so TEforest may have reduced performance on consensus sequences with high IUPAC content.  
  - Current reference BED usage in inference: columns 1/2/3 are used as genomic coordinates, and column 7 is used as the TE family ID.
  - Columns 4/5/6 and any trailing columns are accepted but are not used by the pipeline.
  - BED can be tab-delimited or whitespace-delimited.
- **`--model`**: Path to the non-reference model (optional). If omitted, TEforest auto-selects a model based on the pipeline’s chosen downsample target coverage (5X/10X/20X/30X/40X/50X).
- **`--ref_model`**: Path to the reference model (optional). Auto-selection follows the same coverage logic as above.
- **`--fq_base_path`**: Directory containing FASTQ files. TEforest will match common read naming conventions (e.g., `_R1/_R2`, `_1/_2`, `.1/.2`, `R1/R2`, lane tokens like `_L001_R1_001`) as long as the sample name appears in the filename.
- **`--cleanup_intermediates`**: Optional flag to delete large intermediate files after they are used (e.g., `fastp/`, `aligned/`, `downsampled/`, `candidate_regions_data/`). Omit this if you want to keep read alignments or candidate-region BAMs for debugging.
- **`--disable_reference_detection`**: Optional flag to skip reference TE feature-vector creation and reference model prediction. This can be useful for genomes with very large numbers of old reference TEs, where reference detection can dominate runtime.
- **`--dry_run`**: Optional flag to run `snakemake --dry-run` through the wrapper, so you can validate file naming, inputs, and DAG construction before launching compute-heavy jobs.
- **`--samples`**: List of sample identifiers to process (space-separated). Note more than one sample can be run in parallel. 

Input validation (runs for both normal execution and `--dry_run`):
- FASTQs are resolved per sample/read, must exist, be non-empty, and have a valid first FASTQ record.
- Reference BED must have at least 7 whitespace-delimited columns.
- Reference BED chromosome names (column 1) must match sequence headers in `ref_genome`.
- Every TE ID in reference BED column 7 must be present in `consensusTEs` FASTA headers.
- Extra TE families in the consensus FASTA are allowed.

The script will generate:
- A `config.yaml` in your specified `workdir` with all parameters.
- Intermediate files used to run the pipeline
- Output `.bed` files ({sample}_TEforest_bps_nonredundant.bed) for each sample in `output/` within the specified `workdir`.
- Output `.vcf` files ({sample}_TEforest_bps_nonredundant.vcf) for each sample in `output/` within the specified `workdir`.

Automatic model selection outputs:
- `selected_models/{sample}_nonreference.pkl` and `selected_models/{sample}_reference.pkl` (symlinks to the chosen models).
- `selected_models/{sample}_nonreference.txt` and `selected_models/{sample}_reference.txt` with the observed coverage, chosen coverage, and model path.

Coverage and downsampling behavior:
- TEforest first measures average coverage from the aligned BAM.
- If observed coverage is below 5X, the pipeline exits with an error.
- TEforest then always chooses a supported downsample target from `{5, 10, 20, 30, 40, 50}` that does not exceed observed coverage.
- Model auto-selection uses this chosen target coverage (for example, observed 18.3X -> chosen target 10X -> 10X model).

VCF notes:
- Non-reference calls are represented as INS records with ALT `<INS:ME:<TEFAM>>`.
- Present reference TEs are represented as DEL records with GT=0/0 (absence will be GT=1/1 once deletion calling is added).

## Training pipeline (user-provided truth BED)

Use `TEforest_train.py` with the training Snakefile at `workflow/training/Snakefile` to train a model on your own data. This pipeline labels candidate regions by overlap with a truth BED and trains a LightGBM classifier. Candidate regions that do not overlap a truth locus are labeled absent (0). If you are interested in training on one region of the genome,be sure to use the --euchromatin flag to denote the regions you want to include (otherwise all regions of the genome will be considered for candidate regions). 

Truth BED requirements:
- Columns: `seqnames`, `start`, `end`, `zygosity`, `sample` (and optional `te` family).
- `zygosity` values are interpreted as 1, 0.5, or 0; in classifier mode these are mapped to classes 2/1/0.

Example usage:

```bash
python TEforest_train.py \
    --workflow_dir /path/to/TEforest/workflow/training \
    --workdir /path/to/training_workdir \
    --threads 32 \
    --consensusTEs /path/to/consensusTEs.fasta \
    --ref_genome /path/to/reference_genome.fasta \
    --ref_te_locations /path/to/ref_te_locations.bed \
    --euchromatin /path/to/euchromatin.bed \
    --truth_bed /path/to/truth_labels.bed \
    --fq_base_path /path/to/fastqs \
    --samples SAMPLE1 SAMPLE2 \
    --label_mode classifier \
    --train_mode classifier \
    --test_size 0.1 \
    --model_out training/teforest_classifier.pkl
```

Notes:
- `--label_mode classifier` maps 1/0.5/0 to 2/1/0; `regressor` keeps the numeric zygosity values.
- `--test_size` controls the train/test split (0 disables the split).
- Output model path is relative to the workdir unless you give an absolute path.

## Quick install test (example dataset)

A small dataset is bundled under `tests/data/simulate_te_insertions` so you can sanity-check a new install after activating the `TEforest` conda environment.

```bash
# from the repo root
./tests/run_example.sh
```

This runs TEforest on `RL150IS400_rep0_fwd` and compares the breakpoint output to the expected call. The expected breakpoint is:

```text
3L    5001    5005    roo|non-reference|1|RL150IS400_rep0_fwd|TEforest|rp|1    0    .
```

You can also override the work directory:

```bash
./tests/run_example.sh /tmp/teforest_example
```

## Contributing

- We are actively working on improving performance across different use cases and on large datasets.
- Issues and pull requests are welcome.
