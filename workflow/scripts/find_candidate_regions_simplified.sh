#!/usr/bin/env bash
# extract_te_reads_v11.sh  – robust against TE families with **no regions or reads**
# -----------------------------------------------------------------------------
# Purpose: Extract paired‑end reads (and their mates) for each TE family from a  
# BAM aligned to genome + TE consensus contigs, plus a BED of host TE loci.      
# Creates <family>_to_ISO1.bam for each family, skipping families that have zero 
# intervals or zero overlapping reads without aborting the whole run.           
# -----------------------------------------------------------------------------
set -euo pipefail

usage() {
  cat <<EOF
Usage: ${0##*/} \
  -b ALN_BAM            BAM aligned vs genome+TE contigs \
  -c TE_FASTA           FASTA of TE consensus \
  -r REF_TE_BED         BED of host TE annotations (col7 = family) \
  -o OUTDIR             output directory \
  -@ THREADS            samtools sort threads (default 4) \
  -j JOBS               parallel jobs (default 1) \
  -n SAMPLE_NAME        sample prefix (optional)
EOF
  exit 1
}

# defaults
threads=4
jobs=1

# --------------------
# argument processing
# --------------------
while getopts ":b:c:r:o:@:j:n:" opt; do
  case $opt in
    b) aln_bam=$(readlink -f "$OPTARG");;
    c) te_fa=$(readlink -f "$OPTARG");;
    r) ref_te_bed=$(readlink -f "$OPTARG");;
    o) outdir=$(readlink -m "$OPTARG");;
    @) threads=$OPTARG;;
    j) jobs=$OPTARG;;
    n) sample_name=$OPTARG;;
    *) usage;;
  esac
done

for v in aln_bam te_fa ref_te_bed outdir; do
  [[ -z ${!v-} ]] && { echo "Error: -$v is required"; usage; }
  [[ $v != outdir && ! -r ${!v} ]] && { echo "Error: ${!v} not found or unreadable"; exit 2; }
done

mkdir -p "$outdir"
: "${sample_name:=$(basename "${aln_bam%%.*}")}"  # default sample name
[[ -e ${aln_bam}.bai ]] || samtools index "$aln_bam"

# --------------------
# load TE families
# --------------------
mapfile -t families < <(grep '^>' "$te_fa" | sed 's/^>//')
echo "[INFO] Found ${#families[@]} TE families in $te_fa" >&2

tmpdir="$outdir/tmp"
export aln_bam te_fa ref_te_bed outdir threads sample_name tmpdir

work_family() {
  local fam="$1"
  local famtmp="$tmpdir/$fam"            # <<< isolate per‑family tmp to avoid race
  mkdir -p "$famtmp"
  local regions="$famtmp/${fam}.regions.bed"
  local region_bam="$famtmp/${fam}.region.bam"
  local names="$famtmp/${fam}.names.txt"
  local outbam="$outdir/${fam}_to_ISO1.bam"

  # ------------------------------------
  # 1) build interval list => $regions
  # ------------------------------------
  # (a) consensus contigs whose name *is* the family or ends with it after a dash)
  samtools view -H "$aln_bam" | \
    awk -v f="$fam" 'BEGIN{OFS="\t"} \
      $1=="@SQ" {split($2,a,":" ); ctg=a[2]; split($3,b,":" ); len=b[2]; norm=ctg; gsub(/-/,"_",norm); \
                   if(norm==f || norm ~ f"$") print ctg,0,len }' > "$regions"

  # (b) host TE intervals (col7 == fam)
  awk -v f="$fam" 'BEGIN{OFS="\t"} $7==f{print $1,$2,$3}' "$ref_te_bed" >> "$regions"

  # If no intervals at all, bail out early – nothing to do for this family.
  if [[ ! -s $regions ]]; then
    echo "[WARN] $fam: no intervals found – skipping." >&2
    rm -rf "$famtmp"
    return 0           # success so GNU parallel continues
  fi

  echo "[DEBUG] Regions for $fam:" >&2
  cat "$regions" >&2

  # ------------------------------------
  # 2) extract overlapping reads
  # ------------------------------------
  if ! samtools view -b -L "$regions" "$aln_bam" > "$region_bam" 2>/dev/null; then
    echo "[WARN] $fam: samtools view failed – skipping." >&2
    rm -rf "$famtmp"
    return 0
  fi

  # quick count – if zero, clean up and leave
  local read_count
  read_count=$(samtools view -c "$region_bam")
  if [[ $read_count -eq 0 ]]; then
    echo "[WARN] $fam: 0 overlapping reads – skipping." >&2
    rm -rf "$famtmp"
    return 0
  fi

  # ------------------------------------
  # 3) collect unique read names & mates
  # ------------------------------------
  samtools view "$region_bam" | cut -f1 | sort -u > "$names"
  echo "[INFO] $fam: $read_count overlapping reads found" >&2

  echo "[INFO] Fetching both mates for $fam" >&2
  samtools view -h "$aln_bam" | \
    awk 'NR==FNR{names[$1]; next} /^@/ {print; next} $1 in names{print}' "$names" - | \
    samtools view -b -@ "$threads" - | samtools sort -@ "$threads" -o "$outbam" -
  samtools index "$outbam"

  echo "[DONE] $fam -> $(basename "$outbam")" >&2
  rm -rf "$famtmp"
}
export -f work_family

echo "[INFO] Launching parallel jobs with -j $jobs" >&2
printf '%s\n' "${families[@]}" | parallel -j "$jobs" work_family {}

touch "$outdir/completed.txt"
echo "[ALL DONE] Results in $outdir" >&2
