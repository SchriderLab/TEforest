#!/bin/bash
# grid_simulation_slurm.sh
# This script creates a grid of TE insertion simulations for a fake genome (small_chrom.fasta)
# using a single roo TE insertion at the same candidate site.
# For each combination of read length and insert size, a simulation (and analysis)
# is submitted as a separate Slurm job.
#
# Each job will request 16 cores and run for up to 12 hours.
# This version runs all callers except ngs_te_mapper2 and relocate2.
#
# Before running, edit the following paths and file names as needed.

# Set paths to required files and directories
MC_DIR="/nas/longleaf/home/adaigle/work/mcclintock"   # absolute path to your local mcclintock repo
OUT_DIR="/nas/longleaf/home/adaigle/work/mcclintock_stuff/simulate_te_insertions_sd200"  # base directory for simulation output
REF="/nas/longleaf/home/adaigle/work/mcclintock_stuff/simulate_te_insertions/small_chrom.fasta"   # your small fake genome
CONSENSUS="/nas/longleaf/home/adaigle/work/mcclintock_stuff/simulate_te_insertions/roo.fasta"       # your roo TE consensus sequence
BED="/nas/longleaf/home/adaigle/work/mcclintock_stuff/simulate_te_insertions/targets.bed"           # candidate insertion site (tab-delimited)
GFF="/nas/longleaf/home/adaigle/work/mcclintock_stuff/simulate_te_insertions/refTEs.gff"            # TE locations GFF
TAXONOMY="/nas/longleaf/home/adaigle/work/mcclintock_stuff/simulate_te_insertions/taxonomy.tsv"       # TE taxonomy file

# Define the parameter grid:
#read_lengths=(50 75 100 125 150)      # for testing; adjust as needed
insert_sizes=(200 250 300 400 500 600)     # adjust as needed
read_lengths=(50 75 100 125 150)      # for testing; adjust as needed
#insert_sizes=(200 350 450)     # adjust as needed
#read_lengths=(150)      # for testing; adjust as needed
#insert_sizes=(400)     # adjust as needed
num_reps=10             # number of replicates per combination

# Loop over each read length and insert size combination
for rl in "${read_lengths[@]}"; do
  for ins in "${insert_sizes[@]}"; do

    # Create a separate output directory for this simulation grid point.
    SIM_DIR="${OUT_DIR}/RL${rl}_IS${ins}"
    mkdir -p "${SIM_DIR}"

    # Define the path to the config file for this grid point.
    CONFIG="${SIM_DIR}/config.json"

    # Build a one-line job command string.
    # After generating the config file with make_snakemake_config.py,
    # we use sed to remove ngs_te_mapper2 and relocate2 from the methods field.
    JOB_CMD="source \$(conda info --base)/etc/profile.d/conda.sh; conda activate mcc_sim; \
python ${MC_DIR}/simulation/make_snakemake_config.py --covrange 5 10 20 30 40 50 --error 0 --family roo --bed ${BED} --mcc ${MC_DIR} --out ${CONFIG} --outdir ${SIM_DIR} --reference ${REF} --consensus ${CONSENSUS} --locations ${GFF} --taxonomy ${TAXONOMY} --length ${rl} --insert ${ins} --threads 16 --numrep ${num_reps}; \
sed -i 's/ngs_te_mapper2,//g; s/relocate2,//g' ${CONFIG}; \
snakemake --snakefile ${MC_DIR}/simulation/Snakefile_sim --configfile ${CONFIG} --cores 16 --use-conda; \
snakemake --snakefile ${MC_DIR}/simulation/Snakefile_analysis --configfile ${CONFIG} --cores 16 --use-conda"

    # Submit the job to Slurm.
    sbatch -p general -n 16 -t 2-12:00:00 \
      --mail-type=BEGIN,END,FAIL \
      --mail-user=adaigle@email.unc.edu \
      --job-name="RL${rl}_IS${ins}" \
      --wrap "$JOB_CMD"

    echo "Submitted Slurm job for RL=${rl}, IS=${ins}."
  done
done
