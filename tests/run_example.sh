#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT}/tests/data/simulate_te_insertions"
WORKDIR="${1:-${ROOT}/tests/workdir/simulate_te_insertions}"
CONFIG="${WORKDIR}/config.yaml"
THREADS="${TEFOREST_TEST_THREADS:-1}"

mkdir -p "${WORKDIR}/output"

cat > "${CONFIG}" <<EOF
threads: ${THREADS}
consensusTEs: ${DATA_DIR}/roo.fasta
ref_genome: ${DATA_DIR}/small_chrom.fasta
ref_te_locations: ${DATA_DIR}/ref_tes.bed
euchromatin: ${DATA_DIR}/small_chrom.txt
model: ${ROOT}/workflow/models/teforest_nonreference_50X.pkl
ref_model: ${ROOT}/workflow/models/teforest_reference_50X.pkl
bam_path: /na/
fq_base_path: ${DATA_DIR}/fastq/
target_coverage: 50
EOF

TARGETS=(
  "output/RL150IS400_rep0_fwd_TEforest_nonredundant.bed"
  "output/RL150IS400_rep0_fwd_TEforest_bps_nonredundant.vcf"
)

snakemake "${TARGETS[@]}" \
  --cores "${THREADS}" \
  -s "${ROOT}/workflow/Snakefile" \
  --configfile "${CONFIG}" \
  --unlock \
  --profile none \
  --workflow-profile none \
  --software-deployment-method env-modules \
  --directory "${WORKDIR}"

snakemake "${TARGETS[@]}" \
  --cores "${THREADS}" \
  -s "${ROOT}/workflow/Snakefile" \
  --configfile "${CONFIG}" \
  --profile none \
  --workflow-profile none \
  --software-deployment-method env-modules \
  --directory "${WORKDIR}"

OUTPUT_BPS="${WORKDIR}/output/RL150IS400_rep0_fwd_TEforest_bps_nonredundant.bed"
EXPECTED_BPS="${DATA_DIR}/expected/RL150IS400_rep0_fwd_TEforest_bps_nonredundant.bed"
OUTPUT_VCF="${WORKDIR}/output/RL150IS400_rep0_fwd_TEforest_bps_nonredundant.vcf"
EXPECTED_VCF="${DATA_DIR}/expected/RL150IS400_rep0_fwd_TEforest_bps_nonredundant.vcf"

if [[ ! -f "${OUTPUT_BPS}" ]]; then
  echo "Missing output: ${OUTPUT_BPS}" >&2
  exit 1
fi

diff -u "${EXPECTED_BPS}" "${OUTPUT_BPS}"
diff -u "${EXPECTED_VCF}" "${OUTPUT_VCF}"
echo "OK: BED and VCF outputs match expected calls."
