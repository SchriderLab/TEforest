#!/usr/bin/env bash
# SBATCH HEADER
#SBATCH -J te_benchmark
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=168
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adaigle@email.unc.edu
#SBATCH -o slurm-%x-%j.out
#SBATCH -e slurm-%x-%j.err

set -uo pipefail

########################
# Config (paths/threads)
########################
PYTHON_BIN="python"                       # or "python3" if needed
THREADS_PER_RUN=8
TOTAL_CORES="${SLURM_CPUS_ON_NODE:-168}"
MAX_JOBS=$(( TOTAL_CORES / THREADS_PER_RUN ))
if (( MAX_JOBS < 1 )); then MAX_JOBS=1; fi

# Force libraries to respect our 8 threads/run
export OMP_NUM_THREADS=$THREADS_PER_RUN
export OPENBLAS_NUM_THREADS=$THREADS_PER_RUN
export MKL_NUM_THREADS=$THREADS_PER_RUN
export NUMEXPR_NUM_THREADS=$THREADS_PER_RUN

# Conda env names
MC_CONDA_ENV="mcclintock"
TE_CONDA_ENV="TEforest"

# Base output & logs
BASE_SPEED="/work/users/a/d/adaigle/mcclintock_stuff/speed_test"
MC_OUT_BASE="${BASE_SPEED}/mcclintock_out"
TE_OUT_BASE="${BASE_SPEED}/teforest_out"
LOG_DIR="${BASE_SPEED}/bench_logs"
RESULTS_FILE="${BASE_SPEED}/benchmark_results.tsv"

mkdir -p "$MC_OUT_BASE" "$TE_OUT_BASE" "$LOG_DIR"

# Ensure GNU time is available
TIME_BIN="/usr/bin/time"
if ! "$TIME_BIN" -v true &>/dev/null; then
  echo "ERROR: GNU time not found at /usr/bin/time (or -v not supported)."
  echo "Please load the correct module (e.g., 'module load anaconda; module load gnu time') or adjust TIME_BIN."
  exit 1
fi

# Make sure 'conda' command is usable (for conda run)
if ! command -v conda &>/dev/null; then
  for p in "$HOME/miniconda3/etc/profile.d/conda.sh" "$HOME/anaconda3/etc/profile.d/conda.sh"; do
    [[ -f "$p" ]] && source "$p" && break
  done
fi
if ! command -v conda &>/dev/null; then
  echo "ERROR: 'conda' command not found. Load your conda module or source conda.sh before running."
  exit 2
fi

#################################
# mcclintock fixed input paths
#################################
MC_SCRIPT="/nas/longleaf/home/adaigle/work/mcclintock/mcclintock.py"
MC_REF="/nas/longleaf/home/adaigle/work/mcclintock_stuff/african_sfs/data/ISO1_GCF_000001215.4_Release_6_prepped.fasta"
MC_CONS="/nas/longleaf/home/adaigle/work/mcclintock_stuff/african_sfs/data/consensusTEs.fasta"
MC_GFF="/nas/longleaf/home/adaigle/work/mcclintock_stuff/african_sfs/data/ISO1.gff3"
MC_TAX="/nas/longleaf/home/adaigle/work/mcclintock_stuff/african_sfs/data/taxonomy.tsv"
MC_R1="/nas/longleaf/home/adaigle/work/test_TEforest/test_lightgbm_50X_newfilter_dynamiclength/fastp_het/A2_A3.R1.fq"
MC_R2="/nas/longleaf/home/adaigle/work/test_TEforest/test_lightgbm_50X_newfilter_dynamiclength/fastp_het/A2_A3.R2.fq"

# mcclintock modes (include temp2 plus the rest)
MC_MODES=( temp2 temp retroseq popoolationte popoolationte2 teflon tebreak )

#################################
# TEforest fixed input paths
#################################
TE_SCRIPT="TEforest.py"  # assumes it's on PATH; set full path if not
TE_WORKFLOW_DIR="/nas/longleaf/home/adaigle/TEforest/workflow"
TE_CONS="/nas/longleaf/home/adaigle/work/mcclintock_stuff/from_the_ashes/data/consensusTEs.fasta"
TE_REF_GENOME="/nas/longleaf/home/adaigle/work/mcclintock_stuff/ISO1_GCF_000001215.4_Release_6_prepped/genome_fasta/ISO1_GCF_000001215.4_Release_6_prepped_unaugmented.fasta"
TE_REF_TE="/nas/longleaf/home/adaigle/Rech_updated_supplemental/DeNovoCoordinates/ISO1.bed"
TE_EUCH="/nas/longleaf/home/adaigle/work/mcclintock_stuff/fullchrom.txt"
TE_MODEL="/nas/longleaf/home/adaigle/work/test_TEforest/test_lightgbm_50X_newfilter_dynamiclength/3L3RX/svrf_classifier_all.pkl"
TE_REF_MODEL="/nas/longleaf/home/adaigle/TEforest/workflow/models/teforest_reference_50X.pkl"
TE_FQ_BASE="/nas/longleaf/home/adaigle/work/test_TEforest/test_lightgbm_50X_newfilter_dynamiclength/fastp_het/"
TE_SAMPLES="A2_A3"

#################################
# Results header
#################################
if [[ ! -s "$RESULTS_FILE" ]]; then
  echo -e "tool\tmode\treplicate\tthreads\tstart_iso\tend_iso\telapsed_seconds\telapsed_hms\tmax_rss_kb\tmax_rss_gb\texit_code\tworkdir" > "$RESULTS_FILE"
fi

#################################
# Helpers
#################################
hmstoseconds() {
  local t="$1"
  if [[ "$t" =~ ^([0-9]+):([0-9]{2}):([0-9]{2}(\.[0-9]+)?)$ ]]; then
    awk -v h="${BASH_REMATCH[1]}" -v m="${BASH_REMATCH[2]}" -v s="${BASH_REMATCH[3]}" 'BEGIN{printf "%.0f", h*3600 + m*60 + s}'
  elif [[ "$t" =~ ^([0-9]+):([0-9]{2}(\.[0-9]+)?)$ ]]; then
    awk -v m="${BASH_REMATCH[1]}" -v s="${BASH_REMATCH[2]}" 'BEGIN{printf "%.0f", m*60 + s}'
  else
    awk -v s="$t" 'BEGIN{printf "%.0f", s}'
  fi
}

# Now accepts the conda env name as the first arg
run_and_record() {
  local env_name="$1"
  local tool="$2"
  local mode="$3"
  local rep="$4"
  local workdir="$5"
  shift 5
  local cmd=( "$@" )

  mkdir -p "$workdir"

  local start_iso end_iso start_epoch end_epoch exit_code elapsed_hms elapsed_sec max_kb max_gb
  start_iso="$(date -Is)"
  start_epoch="$(date +%s)"

  local time_log="${LOG_DIR}/${tool}_${mode}_rep${rep}.time.txt"
  local out_log="${LOG_DIR}/${tool}_${mode}_rep${rep}.stdout.log"
  local err_log="${LOG_DIR}/${tool}_${mode}_rep${rep}.stderr.log"

  # Run inside the requested conda env
  # --no-capture-output avoids extra buffering that can interfere with logs
  $TIME_BIN -v -o "$time_log" conda run -n "$env_name" --no-capture-output "${cmd[@]}" >"$out_log" 2>"$err_log"
  exit_code=$?

  end_iso="$(date -Is)"
  end_epoch="$(date +%s)"

  if [[ -f "$time_log" ]]; then
    elapsed_hms="$(grep -F 'Elapsed (wall clock) time' "$time_log" | awk -F': ' '{print $2}')"
    max_kb="$(grep -F 'Maximum resident set size (kbytes)' "$time_log" | awk -F': ' '{gsub(/[[:space:]]*/,"",$2); print $2}')"
  fi

  [[ -z "${elapsed_hms:-}" ]] && elapsed_hms=""
  [[ -z "${max_kb:-}" ]] && max_kb=""

  if [[ -n "$elapsed_hms" ]]; then
    elapsed_sec="$(hmstoseconds "$elapsed_hms")"
  else
    elapsed_sec=$(( end_epoch - start_epoch ))
  fi

  if [[ -n "$max_kb" ]]; then
    max_gb="$(awk -v kb="$max_kb" 'BEGIN{printf "%.3f", kb/1048576.0}')"
  else
    max_gb=""
  fi

  echo -e "${tool}\t${mode}\t${rep}\t${THREADS_PER_RUN}\t${start_iso}\t${end_iso}\t${elapsed_sec}\t${elapsed_hms}\t${max_kb}\t${max_gb}\t${exit_code}\t${workdir}" >> "$RESULTS_FILE"

  # Clean up output dir on success
  if [[ "$exit_code" -eq 0 ]]; then
    rm -rf "$workdir"
  else
    echo "Run ${tool}/${mode} rep ${rep} failed (exit ${exit_code}). Kept workdir: ${workdir}" >&2
  fi
}

export -f hmstoseconds run_and_record
export TIME_BIN LOG_DIR RESULTS_FILE THREADS_PER_RUN

########################
# Queue the jobs
########################
running=0

# 1) mcclintock modes (3 reps each)
for mode in "${MC_MODES[@]}"; do
  for rep in 1 2 3; do
    workdir="${MC_OUT_BASE}/${mode}_${rep}"
    cmd=( "$PYTHON_BIN" "$MC_SCRIPT"
          -p "$THREADS_PER_RUN"
          -r "$MC_REF"
          -c "$MC_CONS"
          -g "$MC_GFF"
          -t "$MC_TAX"
          -1 "$MC_R1"
          -2 "$MC_R2"
          -o "$workdir"
          -a "$MC_CONS"
          -m "$mode"
        )
    bash -c 'run_and_record "$@"' _ "$MC_CONDA_ENV" mcclintock "$mode" "$rep" "$workdir" "${cmd[@]}" &
    ((++running >= MAX_JOBS)) && wait -n && ((running--))
  done
done

# 2) TEforest (3 reps)
for rep in 1 2 3; do
  workdir="${TE_OUT_BASE}/teforest_${rep}"
  cmd=( "$PYTHON_BIN" "$TE_SCRIPT"
        --workflow_dir "$TE_WORKFLOW_DIR"
        --workdir "$workdir"
        --threads "$THREADS_PER_RUN"
        --consensusTEs "$TE_CONS"
        --ref_genome "$TE_REF_GENOME"
        --ref_te_locations "$TE_REF_TE"
        --euchromatin "$TE_EUCH"
        --model "$TE_MODEL"
        --ref_model "$TE_REF_MODEL"
        --fq_base_path "$TE_FQ_BASE"
        --samples "$TE_SAMPLES"
      )
  bash -c 'run_and_record "$@"' _ "$TE_CONDA_ENV" teforest teforest "$rep" "$workdir" "${cmd[@]}" &
  ((++running >= MAX_JOBS)) && wait -n && ((running--))
done

# Wait for all jobs to finish
wait

echo "All done. Results written to: ${RESULTS_FILE}"
echo "Per-run logs in: ${LOG_DIR}"
echo "Max concurrent jobs: ${MAX_JOBS} (derived from ${TOTAL_CORES} cores / ${THREADS_PER_RUN} threads each)"
