#!/bin/bash
#SBATCH -J deepmei
#SBATCH -p general
#SBATCH -N 1
#SBATCH -c 64
#SBATCH -t 12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adaigle@email.unc.edu
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

set -Eeuo pipefail

# --- modules
module purge
module load apptainer
module load samtools
module load bwa || true

# --- inputs (all under /nas -> /work)
APP_SANDBOX="/nas/longleaf/home/adaigle/work/deepmei/deepmei_app"  # writable apptainer sandbox
BAM="/nas/longleaf/home/adaigle/work/human_TEs/megane_output/mapping/NA12878_highcov.bam"
MEI_FA="/nas/longleaf/home/adaigle/work/human_TEs/mcclintock_inputs/human_tes.fa"
REF_FA="/nas/longleaf/home/adaigle/work/human_TEs/mcclintock_inputs/hg38.fa"

THREADS=64
RUNROOT="/nas/longleaf/home/adaigle/work/deepmei/runs/${SLURM_JOB_ID}"
OUTDIR="${RUNROOT}/out"
SCRATCH="${RUNROOT}/scratch"

BAM_DIR=$(dirname "$BAM")
MEI_DIR=$(dirname "$MEI_FA")
REF_DIR=$(dirname "$REF_FA")
BAM_BN=$(basename "$BAM")
MEI_BN=$(basename "$MEI_FA")
REF_BN=$(basename "$REF_FA")

# --- host prep
mkdir -p "${OUTDIR}" "${SCRATCH}/tmp" "${SCRATCH}/cache" "${SCRATCH}/home"
mkdir -p "${RUNROOT}/tmp" "${RUNROOT}/home" "${RUNROOT}/ref"

[ -s "$BAM" ]     || { echo "[FATAL] missing BAM: $BAM" >&2; exit 1; }
[ -s "${BAM}.bai" ] || { echo "[INFO] indexing BAM..."; samtools index "$BAM"; }
[ -s "$MEI_FA" ]  || { echo "[FATAL] missing MEI fasta: $MEI_FA" >&2; exit 1; }
[ -s "$REF_FA" ]  || { echo "[FATAL] missing reference: $REF_FA" >&2; exit 1; }

# keep apptainer off random binds and off /tmp
export TMPDIR="${SCRATCH}/tmp"
export APPTAINER_TMPDIR="${SCRATCH}/tmp"
export APPTAINER_CACHEDIR="${SCRATCH}/cache"
unset APPTAINER_BINDPATH SINGULARITY_BINDPATH SINGULARITY_BIND APPTAINER_BIND
export APPTAINER_BINDPATH= SINGULARITY_BINDPATH= SINGULARITY_BIND= APPTAINER_BIND=

trap '
  echo "[FATAL] line $LINENO failed"
  df -h "${RUNROOT}" "${OUTDIR}" "${SCRATCH}" || true
  df -ih "${RUNROOT}" "${OUTDIR}" "${SCRATCH}" || true
' ERR

# --- run
apptainer exec --cleanenv --no-home --writable \
  --bind "${BAM_DIR}:/in_bam","${MEI_DIR}:/in_mei","${REF_DIR}:/ref","${RUNROOT}:/opt/data" \
  --pwd /opt/data \
  "$APP_SANDBOX" /bin/bash -lc "
    set -Eeuo pipefail
    say(){ echo \"[\$(date \"+%F %T\")] \$*\"; }
    fail(){ echo \"[\$(date \"+%F %T\")] [FATAL] \$*\" >&2; exit 99; }
    need(){ [ -s \"\$1\" ] || fail \"missing file: \$1\"; }

    # dirs + tmp
    export HOME=/opt/data/home
    export TMPDIR=/opt/data/tmp
    export TMP=/opt/data/tmp
    export TEMP=/opt/data/tmp
    export SORT_TMPDIR=/opt/data/tmp
    mkdir -p /opt/data/tmp /opt/data/home /opt/data/out /opt/data/ref

    # tools
    command -v DeepMEI >/dev/null || fail \"DeepMEI not on PATH\"
    command -v samtools >/dev/null || fail \"samtools not in container\"
    command -v bwa >/dev/null || fail \"bwa not in container\"

    # reference + sidecars (local, writable)
    cp -f /ref/${REF_BN} /opt/data/ref/hg38.fa
    need /opt/data/ref/hg38.fa

    if [ ! -s /opt/data/ref/hg38.fa.fai ]; then
      say \"faidx ref\"
      samtools faidx /opt/data/ref/hg38.fa
    fi
    if [ ! -s /opt/data/ref/hg38.fa.dict ]; then
      say \"dict ref\"
      samtools dict /opt/data/ref/hg38.fa -o /opt/data/ref/hg38.fa.dict
      ln -sf hg38.fa.dict /opt/data/ref/hg38.dict
    fi
    # recommended: BWA index ref
    if [ ! -s /opt/data/ref/hg38.fa.bwt ]; then
      say \"bwa index ref (once)\"
      bwa index /opt/data/ref/hg38.fa
    fi

    # MEI fasta -> fastq (+index)
    cp -f /in_mei/${MEI_BN} /opt/data/human_tes.fa
    need /opt/data/human_tes.fa
    if [ ! -s /opt/data/human_tes.fastq ] || [ /opt/data/human_tes.fa -nt /opt/data/human_tes.fastq ]; then
      awk '
        /^>/ { if (seq) { print \"@\"name; print seq; print \"+\"; gsub(/./,\"~\",seq); print seq } name=substr(\$0,2); seq=\"\"; next }
        { seq=seq \$0 }
        END { if (seq) { print \"@\"name; print seq; print \"+\"; gsub(/./,\"~\",seq); print seq } }
      ' /opt/data/human_tes.fa > /opt/data/human_tes.fastq
    fi
    need /opt/data/human_tes.fastq
    [ -s /opt/data/human_tes.fastq.bwt ] || bwa index /opt/data/human_tes.fastq

    # bam sanity
    need /in_bam/${BAM_BN}
    [ -s /in_bam/${BAM_BN}.bai ] || { say \"index BAM (in-container)\"; samtools index /in_bam/${BAM_BN}; }
    if ! samtools quickcheck -v /in_bam/${BAM_BN} > /opt/data/quickcheck.log 2>&1; then
      sed \"s/^/[quickcheck] /\" /opt/data/quickcheck.log >&2
      fail \"BAM failed quickcheck\"
    fi

    # quick contig peek (no awk quoting hell)
    samtools view -H /in_bam/${BAM_BN} | grep -m1 \"^@SQ\" | sed -E 's/.*SN:([^ \\t]+).*/[SN] BAM: \\1/' || true
    grep -m1 \"^>\" /opt/data/ref/hg38.fa | sed 's/^/[SN] REF: /' || true

    # run DeepMEI
    say \"DeepMEI start\"
    set +e
    DeepMEI \
      -i /in_bam/${BAM_BN} \
      -r /opt/data/ref/hg38.fa \
      -m /opt/data/human_tes.fastq \
      -w /opt/data/out \
      -p ${THREADS}
    rc=\$?
    set -e
    say \"DeepMEI exit=\$rc\"

    # must produce something
    OUTBASE=\"/opt/data/out/DeepMEI_output/${BAM_BN}\"
    if [ ! -d \"\$OUTBASE\" ] || ! compgen -G \"\$OUTBASE/*\" >/dev/null; then
      fail \"no outputs under \$OUTBASE (check logs above)\"
    fi

    # surface any final vcfs if DeepMEI uses its own dir
    [ -d /opt/DeepMEI/final_vcf ] && find /opt/DeepMEI/final_vcf -maxdepth 2 -type f -ls | sed \"s,^,[final_vcf] ,\" || true

    say \"done\"
  "

echo "[HOST] outputs -> ${OUTDIR}"
