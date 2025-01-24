# TEforest

This repository contains a Snakemake pipeline for training, deploying, and benchmarking TE (transposable element) insertion models based on Random Forest classifiers. 

## Overview

- **Snakemake Pipeline**  
  - The main pipeline is controlled by `workflow/Snakefile`.  
  - It orchestrates rules for:
    - Data preprocessing
    - Model training
    - Model deployment (prediction)
    - Benchmarking
  - Only scripts used for model deployment are designed to be run on another machine. The training and benchmarking scripts contain hard-coded paths for now. 

- **Scripts**  
  - Supporting scripts for training and model deployment are located in `workflow/scripts`.

- **Pre-Trained Models**  
  - Trained Random Forest models are stored in `workflow/models`.

- **Outputs**  
  - Final output files (predictions) will be generated in your chosen work directory under the path:
    ```bash
    outputs/{sample}_TEforest_nonredundant.bed
    ```

## Environment Setup

A Conda environment YAML file (`TEforest.yml`) is included for convenience. It defines the base dependencies for running this pipeline. However, note that **additional R packages** such as **GenomicAlignments** and **GenomicRanges** should be installed to ensure full functionality. 

```bash
# Example environment creation
conda env create -f TEforest.yml

# Activate the environment
conda activate TEforest

# Then install additional R packages within this environment.
```

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
    --samples A1 A2 A3
```

- **`--workflow_dir`**: Directory containing the `Snakefile` (`workflow/Snakefile`).
- **`--workdir`**: Directory to store outputs and logs.
- **`--threads`**: Number of CPU threads to use. 16 per sample is recommended.
- **`--consensusTEs`, `--ref_genome`, `--ref_te_locations`, `--euchromatin`**: Input reference files for TE detection. All calls outside of the regions denoted in euchromatin will be filtered. 
- **`--model`**: Path to the main trained model (Random Forest).
- **`--ref_model`**: Path to the reference model (e.g., for control or comparison).
- **`--fq_base_path`**: Directory containing FASTQ files. Should contain sample in name like {sample}_1.fastq.gz and {sample}_2.fastq.gz or {sample}_1.fq.gz.
- **`--samples`**: List of sample identifiers to process (space-separated).

The script will generate:
- A `config.yaml` in your specified `workdir` with all parameters.
- Intermediate files used to run the pipeline
- Output `.bed` files ({sample}_TEforest_bps_nonredundant.bed) for each sample in `outputs/` within the specified `workdir`.

## Contributing

- We are actively working on **usability improvements** (e.g., streamlined CLI arguments, a fully functional conda installation).
- This pipeline has only been tested in Drosophila melanogaster. Future tests and updates will ensure useability in other species. 
- Issues and pull requests are welcome.
