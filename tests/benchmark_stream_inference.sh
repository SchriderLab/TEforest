#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOURCE_WORKDIR="${1:-${ROOT}/tests/workdir/simulate_te_insertions}"
SAMPLE="${2:-RL150IS400_rep0_fwd}"
REPEATS="${3:-10 100}"

CSV_SRC="${SOURCE_WORKDIR}/featvec_csvs/${SAMPLE}.csv"
BAM_DIR="${SOURCE_WORKDIR}/downsampled"
TE_BAM_DIR="${SOURCE_WORKDIR}/candidate_regions_data"
MODEL="${ROOT}/workflow/models/teforest_nonreference_50X.pkl"

OLD_BAM_TO_FVEC="${ROOT}/workflow/scripts/bam_to_fvec_cov_normalized_rescale_regions_majorchanges.py"
OLD_CONDENSE="${ROOT}/workflow/scripts/condense_training_data.py"
OLD_USE_MODEL="${ROOT}/workflow/scripts/use_model.py"
NEW_STREAM="${ROOT}/workflow/scripts/predict_from_regions_stream.py"

BENCH_DIR="${ROOT}/tests/workdir/stream_benchmark"
mkdir -p "${BENCH_DIR}"
SUMMARY="${BENCH_DIR}/summary.tsv"
printf "repeat\tmethod\trows\telapsed\tmax_rss_kb\n" > "${SUMMARY}"

if [[ ! -f "${CSV_SRC}" ]]; then
  echo "Missing source CSV: ${CSV_SRC}" >&2
  exit 1
fi

for n in ${REPEATS}; do
  REP_CSV="${BENCH_DIR}/${SAMPLE}.x${n}.csv"
  awk -F',' -v OFS=',' -v n="${n}" '
    NR==1 { print; next }
    NR==2 {
      for (i=0; i<n; i++) {
        print $1, $2, $3 + i, $4 + i, $5, $6
      }
      exit
    }
  ' "${CSV_SRC}" > "${REP_CSV}"

  OLD_DIR="${BENCH_DIR}/old_x${n}"
  NEW_DIR="${BENCH_DIR}/new_x${n}"
  rm -rf "${OLD_DIR}" "${NEW_DIR}"
  mkdir -p "${OLD_DIR}" "${NEW_DIR}"

  cat > "${OLD_DIR}/run_old.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
python "${OLD_BAM_TO_FVEC}" \\
  -i "${REP_CSV}" \\
  -od "${OLD_DIR}/feature_vectors" \\
  -bd "${BAM_DIR}" \\
  -tebd "${TE_BAM_DIR}" \\
  -p 1
python "${OLD_CONDENSE}" \\
  -i "${OLD_DIR}/feature_vectors/${SAMPLE}" \\
  -o "${OLD_DIR}/condensed.npz"
python "${OLD_USE_MODEL}" \\
  -n "${OLD_DIR}/condensed.npz" \\
  -m "${MODEL}" \\
  -o "${OLD_DIR}/predictions" \\
  -t 1
EOF
  chmod +x "${OLD_DIR}/run_old.sh"

  /usr/bin/time -v -o "${OLD_DIR}/time.txt" "${OLD_DIR}/run_old.sh"

  /usr/bin/time -v -o "${NEW_DIR}/time.txt" \
    python "${NEW_STREAM}" \
      -i "${REP_CSV}" \
      -bd "${BAM_DIR}" \
      -tebd "${TE_BAM_DIR}" \
      -m "${MODEL}" \
      -o "${NEW_DIR}/predictions.csv" \
      --mode nonreference \
      -t 1 \
      --batch_size 64

  python - <<PY
import pandas as pd
old_df = pd.read_csv("${OLD_DIR}/predictions/predictions.csv")
new_df = pd.read_csv("${NEW_DIR}/predictions.csv")
if len(old_df) != len(new_df):
    raise SystemExit(f"row count mismatch x${n}: old={len(old_df)} new={len(new_df)}")
if not old_df[["file", "true", "pred"]].equals(new_df[["file", "true", "pred"]]):
    raise SystemExit("prediction mismatch for x${n}")
print(f"x${n}: predictions match ({len(old_df)} rows)")
PY

  old_rows=$(($(wc -l < "${OLD_DIR}/predictions/predictions.csv") - 1))
  new_rows=$(($(wc -l < "${NEW_DIR}/predictions.csv") - 1))
  if [[ "${old_rows}" != "${new_rows}" ]]; then
    echo "row mismatch after comparison stage for x${n}" >&2
    exit 1
  fi

  old_elapsed=$(grep "Elapsed (wall clock) time" "${OLD_DIR}/time.txt" | sed 's/.*: //')
  old_rss=$(grep "Maximum resident set size" "${OLD_DIR}/time.txt" | awk -F': ' '{print $2}')
  new_elapsed=$(grep "Elapsed (wall clock) time" "${NEW_DIR}/time.txt" | sed 's/.*: //')
  new_rss=$(grep "Maximum resident set size" "${NEW_DIR}/time.txt" | awk -F': ' '{print $2}')

  printf "%s\told\t%s\t%s\t%s\n" "${n}" "${old_rows}" "${old_elapsed}" "${old_rss}" >> "${SUMMARY}"
  printf "%s\tnew\t%s\t%s\t%s\n" "${n}" "${new_rows}" "${new_elapsed}" "${new_rss}" >> "${SUMMARY}"
done

echo "Wrote benchmark summary: ${SUMMARY}"
cat "${SUMMARY}"
