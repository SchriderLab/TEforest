#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT}/tests/data/simulate_te_insertions"
WORKDIR="${1:-${ROOT}/tests/workdir/simulate_te_insertions}"

mkdir -p "${WORKDIR}"

python "${ROOT}/TEforest.py" \
  --workflow_dir "${ROOT}/workflow" \
  --workdir "${WORKDIR}" \
  --threads 4 \
  --consensusTEs "${DATA_DIR}/roo.fasta" \
  --ref_genome "${DATA_DIR}/small_chrom.fasta" \
  --ref_te_locations "${DATA_DIR}/ref_tes.bed" \
  --euchromatin "${DATA_DIR}/small_chrom.txt" \
  --model "${ROOT}/workflow/models/teforest_nonreference_50X.pkl" \
  --ref_model "${ROOT}/workflow/models/teforest_reference_50X.pkl" \
  --fq_base_path "${DATA_DIR}/fastq/" \
  --samples RL150IS400_rep0_fwd

OUTPUT_BPS="${WORKDIR}/output/RL150IS400_rep0_fwd_TEforest_bps_nonredundant.bed"
EXPECTED_BPS="${DATA_DIR}/expected/RL150IS400_rep0_fwd_TEforest_bps_nonredundant.bed"

if [[ ! -f "${OUTPUT_BPS}" ]]; then
  echo "Missing output: ${OUTPUT_BPS}" >&2
  exit 1
fi

diff -u "${EXPECTED_BPS}" "${OUTPUT_BPS}"
echo "OK: output matches expected breakpoint call."
