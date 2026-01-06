#!/usr/bin/env bash
#SBATCH -p general
#SBATCH --cpus-per-task=64
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adaigle@email.unc.edu
#SBATCH -o logs/%x-%j.out
#SBATCH --begin=now+5hours

set -euo pipefail

# Args: BAM [COVERAGE]
BAM="${1:?Provide a BAM path}"
COVERAGE="${2:-50}"
THREADS="${SLURM_CPUS_PER_TASK:-64}"

module purge
module load samtools seqkit

# Paths
MC_PATH="/nas/longleaf/home/adaigle/work/mcclintock/mcclintock.py"
REF_FA="/nas/longleaf/home/adaigle/work/mcclintock_stuff/data/ISO1_GCF_000001215.4_Release_6_prepped.fasta"
TE_FA="/nas/longleaf/home/adaigle/work/mcclintock_stuff/data/consensusTEs.fasta"
GFF="/nas/longleaf/home/adaigle/work/mcclintock_stuff/data/ISO1.gff3"
TAX="/nas/longleaf/home/adaigle/work/mcclintock_stuff/data/taxonomy.tsv"
METHODS="trimgalore,temp2,teflon,popoolationte2"

OUTROOT="/nas/longleaf/home/adaigle/work/mcclintock_stuff/het_ref_model"
OUTDIR="${OUTROOT}/${COVERAGE}X"
mkdir -p "${OUTDIR}" logs

SAMPLE="$(basename "${BAM}" .bam)"
FQ1="${OUTDIR}/${SAMPLE}_1.fastq"
FQ2="${OUTDIR}/${SAMPLE}_2.fastq"
FQ1D="${OUTDIR}/${SAMPLE}_1.fq.gz"
FQ2D="${OUTDIR}/${SAMPLE}_2.fq.gz"
MCC_INTERMEDIATE_DIR="${OUTDIR}/${SAMPLE}_1/intermediate"

echo "[$(date)] Starting ${SAMPLE} at ${COVERAGE}X -> ${OUTDIR}"

# 1) name-sort
samtools sort -n -@ "${THREADS}" -o "${OUTDIR}/${SAMPLE}.name_sorted.bam" "${BAM}"

# 2) extract paired FASTQs (uncompressed)
samtools fastq -@ "${THREADS}" \
  -1 "${FQ1}" \
  -2 "${FQ2}" \
  -0 /dev/null -s /dev/null -n \
  "${OUTDIR}/${SAMPLE}.name_sorted.bam"
rm -f "${OUTDIR}/${SAMPLE}.name_sorted.bam"

# 3) dedupe -> gz
seqkit rmdup -n -o "${FQ1D}" "${FQ1}"
seqkit rmdup -n -o "${FQ2D}" "${FQ2}"
rm -f "${FQ1}" "${FQ2}"

# 4) McClintock
module purge

# Ensure the conda shell functions are available in non-interactive jobs
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mcclintock

set +e
python "${MC_PATH}" \
  -p "${THREADS}" \
  -r "${REF_FA}" \
  -c "${TE_FA}" \
  -g "${GFF}" \
  -t "${TAX}" \
  -1 "${FQ1D}" \
  -2 "${FQ2D}" \
  -o "${OUTDIR}" \
  -m "${METHODS}" \
  -a "${TE_FA}" --resume
MCC_STATUS=$?
set -e

if [[ "${MCC_STATUS}" -ne 0 ]]; then
  echo "[$(date)] McClintock FAILED for ${SAMPLE} (exit ${MCC_STATUS}). Skipping cleanup."
  exit "${MCC_STATUS}"
fi

# 5) Mandatory cleanup (only after success)
echo "[$(date)] Cleaning ${MCC_INTERMEDIATE_DIR} and FASTQs"
rm -rf "${MCC_INTERMEDIATE_DIR}" || true
rm -f "${FQ1D}" "${FQ2D}" || true

echo "[$(date)] Done: ${SAMPLE} @ ${COVERAGE}X"
