#!/usr/bin/env bash
#SBATCH -J bench_3x
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --cpus-per-task=8
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adaigle@email.unc.edu
#SBATCH -o slurm-%x-%j.out
#SBATCH -e slurm-%x-%j.err

set -euo pipefail

TOOL="${1:-mcclintock}"         # mcclintock | teforest
MC_MODE="${2:-temp2}"           # only used for mcclintock
REPS=3
THREADS=8
TIME_BIN="/usr/bin/time"

# ----- Make conda available for `conda run` -----
if ! command -v conda &>/dev/null; then
  for p in "$HOME/miniconda3/etc/profile.d/conda.sh" "$HOME/anaconda3/etc/profile.d/conda.sh"; do
    [[ -f "$p" ]] && source "$p" && break
  done
fi
command -v conda >/dev/null 2>&1 || { echo "ERROR: 'conda' not found in PATH."; exit 2; }

# ----- Sanity checks -----
"$TIME_BIN" -v true >/dev/null 2>&1 || { echo "ERROR: $TIME_BIN -v not usable"; exit 3; }

# ----- Threaded libs respect our 8 CPUs per replicate -----
export OMP_NUM_THREADS=$THREADS
export OPENBLAS_NUM_THREADS=$THREADS
export MKL_NUM_THREADS=$THREADS
export NUMEXPR_NUM_THREADS=$THREADS

# ====== Paths ======
BASE_SPEED="/work/users/a/d/adaigle/mcclintock_stuff/speed_test"

# McClintock
MC_CONDA_ENV="mcclintock"
MC_SCRIPT="/nas/longleaf/home/adaigle/work/mcclintock/mcclintock.py"
MC_REF="/nas/longleaf/home/adaigle/work/mcclintock_stuff/african_sfs/data/ISO1_GCF_000001215.4_Release_6_prepped.fasta"
MC_CONS="/nas/longleaf/home/adaigle/work/mcclintock_stuff/african_sfs/data/consensusTEs.fasta"
MC_GFF="/nas/longleaf/home/adaigle/work/mcclintock_stuff/african_sfs/data/ISO1.gff3"
MC_TAX="/nas/longleaf/home/adaigle/work/mcclintock_stuff/african_sfs/data/taxonomy.tsv"
MC_R1="/nas/longleaf/home/adaigle/users/test_teforest/hetrefte_30X_otherscript/het_reads/JUT-008_MUN-009/JUT-008_MUN-009_1.fq"
MC_R2="/nas/longleaf/home/adaigle/users/test_teforest/hetrefte_30X_otherscript/het_reads/JUT-008_MUN-009/JUT-008_MUN-009_2.fq"
MC_BASE_OUT="${BASE_SPEED}/mcclintock_out"
mkdir -p "$MC_BASE_OUT"

# TEforest
TE_CONDA_ENV="TEforest"
TE_SCRIPT="/nas/longleaf/home/adaigle/TEforest/TEforest.py"  # set full path if not on PATH
TE_WORKFLOW_DIR="/nas/longleaf/home/adaigle/TEforest/workflow"
TE_CONS="/nas/longleaf/home/adaigle/work/mcclintock_stuff/from_the_ashes/data/consensusTEs.fasta"
TE_REF_GENOME="/nas/longleaf/home/adaigle/work/mcclintock_stuff/ISO1_GCF_000001215.4_Release_6_prepped/genome_fasta/ISO1_GCF_000001215.4_Release_6_prepped_unaugmented.fasta"
TE_REF_TE="/nas/longleaf/home/adaigle/Rech_updated_supplemental/DeNovoCoordinates/ISO1.bed"
TE_EUCH="/nas/longleaf/home/adaigle/work/mcclintock_stuff/fullchrom.txt"
TE_MODEL="/nas/longleaf/home/adaigle/work/test_TEforest/test_lightgbm_30X_newfilter_dynamiclength/3L3RX/svrf_classifier_all.pkl"
TE_REF_MODEL="/nas/longleaf/home/adaigle/TEforest/workflow/models/teforest_reference_30X.pkl"
TE_FQ_BASE="/nas/longleaf/home/adaigle/users/test_teforest/hetrefte_30X_otherscript/het_reads/JUT-008_MUN-009/"
TE_SAMPLES="JUT-008_MUN-009"
TE_BASE_OUT="${BASE_SPEED}/teforest_out"
mkdir -p "$TE_BASE_OUT"

launch_mcclintock_rep() {
  local rep="$1"
  local outdir="${MC_BASE_OUT}/3x_${MC_MODE}_rep${rep}_${SLURM_JOB_ID}"
  local logdir="${outdir}/_logs"
  local timefile="${logdir}/time.txt"
  mkdir -p "$outdir" "$logdir"
  touch "${logdir}/stdout.log" "${logdir}/stderr.log" "${timefile}"

  srun --exclusive -N1 -n1 -c "$THREADS" bash -lc "
    set -euo pipefail
    export OMP_NUM_THREADS=$THREADS OPENBLAS_NUM_THREADS=$THREADS MKL_NUM_THREADS=$THREADS NUMEXPR_NUM_THREADS=$THREADS

    { echo -e 'key\tvalue';
      echo -e 'tool\tmcclintock';
      echo -e 'mode\t${MC_MODE}';
      echo -e 'replicate\t${rep}';
      echo -e 'threads\t${THREADS}';
      echo -e 'outdir\t${outdir}';
      echo -e 'start_iso\t'\"$(date -Is)\"; } > '${logdir}/meta.tsv'

    { echo -e 'key\tvalue';
      echo -e 'hostname\t'\"\$(hostname)\"; 
      echo -e 'date_iso\t'\"$(date -Is)\"; 
      echo -e 'kernel\t'\"\$(uname -r)\"; 
      os='unknown'; [[ -f /etc/os-release ]] && os=\$(grep PRETTY_NAME /etc/os-release | cut -d= -f2 | tr -d '\"'); echo -e 'os\t'\$os;
      echo -e 'slurm_jobid\t${SLURM_JOB_ID}';
      echo -e 'slurm_stepid\t'\${SLURM_STEP_ID:-};
      echo -e 'slurm_node\t'\${SLURMD_NODENAME:-\$(hostname)};
      echo -e 'cpus_on_node\t'\${SLURM_CPUS_ON_NODE:-};
      lscpu | awk -F: 'BEGIN{OFS=\"\t\"} \$1~/Model name|Socket|Core|Thread|CPU\\(s\\)/{gsub(/^[ \\t]+|[ \\t]+$/,\"\",\$1); gsub(/^[ \\t]+|[ \\t]+$/,\"\",\$2); print tolower(gensub(/[ ()]/,\"_\",\"g\",\$1)),\$2}';
      awk 'BEGIN{OFS=\"\t\"} /MemTotal:/{print \"memtotal_kb\",\$2}' /proc/meminfo;
      # Env-specific bits
      conda run -n '${MC_CONDA_ENV}' python -V 2>&1 | awk '{print \"python_version\\t\"\$0}';
      conda run -n '${MC_CONDA_ENV}' which python 2>/dev/null | awk '{print \"python_path\\t\"\$0}';
    } > '${logdir}/specs.tsv' || true

    set +e
    $TIME_BIN -v -o '${timefile}' conda run -n '${MC_CONDA_ENV}' --no-capture-output \
      python '${MC_SCRIPT}' \
        -p '${THREADS}' \
        -r '${MC_REF}' \
        -c '${MC_CONS}' \
        -g '${MC_GFF}' \
        -t '${MC_TAX}' \
        -1 '${MC_R1}' \
        -2 '${MC_R2}' \
        -o '${outdir}' \
        -a '${MC_CONS}' \
        -m '${MC_MODE}' \
        > '${logdir}/stdout.log' 2> '${logdir}/stderr.log'
    exit_code=\$?
    set -e

    echo -e 'end_iso\t'\"$(date -Is)\" >> '${logdir}/meta.tsv'
    echo -e 'exit_code\t'\$exit_code >> '${logdir}/meta.tsv'
    exit \$exit_code
  " &
}

launch_teforest_rep() {
  local rep="$1"
  local outdir="${TE_BASE_OUT}/3x_teforest_rep${rep}_${SLURM_JOB_ID}"
  local logdir="${outdir}/_logs"
  local timefile="${logdir}/time.txt"
  mkdir -p "$outdir" "$logdir"
  touch "${logdir}/stdout.log" "${logdir}/stderr.log" "${timefile}"

  srun --exclusive -N1 -n1 -c "$THREADS" bash -lc "
    set -euo pipefail
    export OMP_NUM_THREADS=$THREADS OPENBLAS_NUM_THREADS=$THREADS MKL_NUM_THREADS=$THREADS NUMEXPR_NUM_THREADS=$THREADS

    { echo -e 'key\tvalue';
      echo -e 'tool\tteforest';
      echo -e 'mode\tteforest';
      echo -e 'replicate\t${rep}';
      echo -e 'threads\t${THREADS}';
      echo -e 'outdir\t${outdir}';
      echo -e 'start_iso\t'\"$(date -Is)\"; } > '${logdir}/meta.tsv'

    { echo -e 'key\tvalue';
      echo -e 'hostname\t'\"\$(hostname)\"; 
      echo -e 'date_iso\t'\"$(date -Is)\"; 
      echo -e 'kernel\t'\"\$(uname -r)\"; 
      os='unknown'; [[ -f /etc/os-release ]] && os=\$(grep PRETTY_NAME /etc/os-release | cut -d= -f2 | tr -d '\"'); echo -e 'os\t'\$os;
      echo -e 'slurm_jobid\t${SLURM_JOB_ID}';
      echo -e 'slurm_stepid\t'\${SLURM_STEP_ID:-};
      echo -e 'slurm_node\t'\${SLURMD_NODENAME:-\$(hostname)};
      echo -e 'cpus_on_node\t'\${SLURM_CPUS_ON_NODE:-};
      lscpu | awk -F: 'BEGIN{OFS=\"\t\"} \$1~/Model name|Socket|Core|Thread|CPU\\(s\\)/{gsub(/^[ \\t]+|[ \\t]+$/,\"\",\$1); gsub(/^[ \\t]+|[ \\t]+$/,\"\",\$2); print tolower(gensub(/[ ()]/,\"_\",\"g\",\$1)),\$2}';
      awk 'BEGIN{OFS=\"\t\"} /MemTotal:/{print \"memtotal_kb\",\$2}' /proc/meminfo;
      # Env-specific bits
      conda run -n '${TE_CONDA_ENV}' python -V 2>&1 | awk '{print \"python_version\\t\"\$0}';
      conda run -n '${TE_CONDA_ENV}' which python 2>/dev/null | awk '{print \"python_path\\t\"\$0}';
    } > '${logdir}/specs.tsv' || true

    set +e
    $TIME_BIN -v -o '${timefile}' conda run -n '${TE_CONDA_ENV}' --no-capture-output \
      python '${TE_SCRIPT}' \
        --workflow_dir '${TE_WORKFLOW_DIR}' \
        --workdir '${outdir}' \
        --threads '${THREADS}' \
        --consensusTEs '${TE_CONS}' \
        --ref_genome '${TE_REF_GENOME}' \
        --ref_te_locations '${TE_REF_TE}' \
        --euchromatin '${TE_EUCH}' \
        --model '${TE_MODEL}' \
        --ref_model '${TE_REF_MODEL}' \
        --fq_base_path '${TE_FQ_BASE}' \
        --samples '${TE_SAMPLES}' \
        > '${logdir}/stdout.log' 2> '${logdir}/stderr.log'
    exit_code=\$?
    set -e

    echo -e 'end_iso\t'\"$(date -Is)\" >> '${logdir}/meta.tsv'
    echo -e 'exit_code\t'\$exit_code >> '${logdir}/meta.tsv'
    exit \$exit_code
  " &
}


echo "Launching 3x '${TOOL}' (mode: ${MC_MODE:-n/a}) with ${THREADS} threads each..."

case "$TOOL" in
  mcclintock)
    for rep in $(seq 1 "$REPS"); do launch_mcclintock_rep "$rep"; done
    ;;
  teforest)
    for rep in $(seq 1 "$REPS"); do launch_teforest_rep "$rep"; done
    ;;
  *)
    echo "Unknown tool: $TOOL (expected 'mcclintock' or 'teforest')" >&2
    exit 10
    ;;
esac

wait
echo "All replicates finished."
