#!/bin/bash
#SBATCH -J megane_dual
#SBATCH -p general
#SBATCH -N 1
#SBATCH -c 256
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adaigle@email.unc.edu

set -Eeuo pipefail
trap 'echo "[FATAL] Line $LINENO failed"; exit 1' ERR

module purge
module load apptainer
module load samtools
module load repeatmasker
module load blast/2.16.0
module load bwa-mem2 || module load bwa
module load blast/2.16.0
module load python
python - <<'PY'
import h5py, sys
print("h5py OK", h5py.version.info)
PY

BWA_DIR="$(dirname "$(which bwa-mem2 || true)")"
export BWA="${BWA_DIR}/bwa-mem2.avx2"
export BWA_MEM2_DISABLE_AVX512=1
$BWA version || true

############################################
# User-configurable paths
############################################
SIF="/nas/longleaf/home/adaigle/work/megane/MEGAnE.sif"
WORKDIR="/nas/longleaf/home/adaigle/work/human_TEs/megane_output"

REF="/nas/longleaf/home/adaigle/work/human_TEs/mcclintock_inputs/hg38.fa"
FQ1="/nas/longleaf/home/adaigle/users/test_TEforest/test_humans_repeatmasker_onlyconsensus/downsampled/NA12878_highcov_1.fq"
FQ2="/nas/longleaf/home/adaigle/users/test_TEforest/test_humans_repeatmasker_onlyconsensus/downsampled/NA12878_highcov_2.fq"
TE_FASTA="/nas/longleaf/home/adaigle/work/human_TEs/mcclintock_inputs/human_tes.fa"  # >ALU, >LINE1, >SVA

SAMPLE="NA12878_highcov"

# Threads
MAP_THREADS=96
RM_THREADS=96
MEGANE_THREADS=4

# Apptainer bind
BIND="-B /nas/longleaf:/nas/longleaf"

############################################
# Layout & logging
############################################
MAPDIR="${WORKDIR}/mapping"
KMERDIR="${WORKDIR}/megane_kmer_set"
CUSTOM_PREP="${WORKDIR}/custom_prep"
BUILTIN_OUT="${WORKDIR}/results_builtin"
CUSTOM_OUT="${WORKDIR}/results_custom"
PATCH_DIR="/nas/longleaf/home/adaigle/work/megane/patches"
LOG="${WORKDIR}/pipeline.log"

mkdir -p "${MAPDIR}" "${KMERDIR}" "${CUSTOM_PREP}" "${BUILTIN_OUT}" "${CUSTOM_OUT}" "${PATCH_DIR}"
echo "[$(date)] Starting pipeline" | tee -a "$LOG"

############################################
# Preflight checks & indexes
############################################
[ -s "$REF" ] || { echo "Missing REF: $REF" | tee -a "$LOG"; exit 1; }
[ -s "$FQ1" ] || { echo "Missing FQ1: $FQ1" | tee -a "$LOG"; exit 1; }
[ -s "$FQ2" ] || { echo "Missing FQ2: $FQ2" | tee -a "$LOG"; exit 1; }
[ -s "$TE_FASTA" ] || { echo "Missing TE_FASTA: $TE_FASTA" | tee -a "$LOG"; exit 1; }

if [ ! -f "${REF}.fai" ]; then
  echo "[$(date)] samtools faidx $REF" | tee -a "$LOG"
  samtools faidx "$REF"
fi

if command -v bwa-mem2 >/dev/null 2>&1; then
  if [ ! -f "${REF}.0123" ]; then
    echo "[$(date)] bwa-mem2 index $REF" | tee -a "$LOG"
    $BWA index "$REF"
  fi
else
  if [ ! -f "${REF}.bwt" ]; then
    echo "[$(date)] bwa index $REF" | tee -a "$LOG"
    bwa index "$REF"
  fi
fi

############################################
# Step 0: build k-mers (once) & custom prep
############################################
if [ ! -f "${KMERDIR}/hg38.mk" ]; then
  echo "[$(date)] build_kmerset" | tee -a "$LOG"
  apptainer exec ${BIND} "$SIF" build_kmerset -fa "$REF" -prefix hg38 -outdir "$KMERDIR" &
  job_kmer=$!
else
  job_kmer=
fi

pushd "$CUSTOM_PREP" >/dev/null
REPOUT_BASENAME="$(basename "$REF").out"            # e.g., hg38.fa.out
FAOUT="${CUSTOM_PREP}/${REPOUT_BASENAME}"

# RepeatMasker on REF with custom library (only if not already done)
if [ ! -f "$REPOUT_BASENAME" ]; then
  echo "[$(date)] RepeatMasker -lib (custom) -> .fa.out" | tee -a "$LOG"
  RepeatMasker -engine rmblast -lib "$TE_FASTA" -s -no_is -a -pa "$RM_THREADS" -dir . "$REF" &
  job_rm=$!
else
  job_rm=
fi

# PolyA lists
cat > ME_with_pA.txt <<'EOF'
ALU
LINE1
SVA
EOF

cat > ME_with_pA.families.txt <<'EOF'
SINE/Alu
LINE/L1
Retroposon/SVA
EOF

: > non_ME_rep.txt

# Main chromosomes (UCSC names)
cat > main_chrs.txt <<'EOF'
chr1
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
chr20
chr21
chr22
chrX
chrY
EOF

# BLAST DB for reference (host builds it; container reads it)
if [ ! -f ref_blastdb.nsq ]; then
  echo "[$(date)] makeblastdb (reference)" | tee -a "$LOG"
  makeblastdb -in "$REF" -out ref_blastdb -dbtype nucl -parse_seqids
fi
popd >/dev/null

############################################
# Step 1: mapping (skipped if bam exists)
############################################
if [ ! -f "${MAPDIR}/${SAMPLE}.bam" ]; then
  echo "[$(date)] Mapping $SAMPLE" | tee -a "$LOG"
  export BWA_MEM2_DISABLE_AVX512=1
  if command -v $BWA >/dev/null 2>&1; then
    $BWA mem -Y -K 100000000 -t "$MAP_THREADS" "$REF" "$FQ1" "$FQ2" \
      | samtools sort -@ "$MAP_THREADS" -o "${MAPDIR}/${SAMPLE}.bam" -
  else
    echo "[WARN] bwa-mem2 not found; falling back to bwa mem" | tee -a "$LOG"
    bwa mem -Y -K 100000000 -t "$MAP_THREADS" "$REF" "$FQ1" "$FQ2" \
      | samtools sort -@ "$MAP_THREADS" -o "${MAPDIR}/${SAMPLE}.bam" -
  fi
  samtools index -@ "$MAP_THREADS" "${MAPDIR}/${SAMPLE}.bam"
fi

# Wait for Step 0 jobs (if any)
[ -n "${job_kmer:-}" ] && wait "$job_kmer" || true
[ -n "${job_rm:-}" ] && wait "$job_rm" || true

############################################
# Step 1.5: Prep custom inputs (no touching originals)
############################################
# A) Patch .fa.out families
if [ -s "$FAOUT" ]; then
  echo "[$(date)] Patching $FAOUT families (Unspecified -> proper)" | tee -a "$LOG"
  cp "$FAOUT" "${FAOUT}.bak"
  awk 'BEGIN{OFS="\t"}
       NR<=3 {print; next}
       ($10=="ALU"   && $11=="Unspecified"){ $11="SINE/Alu" }
       ($10=="LINE1" && $11=="Unspecified"){ $11="LINE/L1" }
       ($10=="SVA"   && $11=="Unspecified"){ $11="Retroposon/SVA" }
       {print}
  ' "${FAOUT}.bak" > "$FAOUT"
else
  echo "[FATAL] RepeatMasker output $FAOUT missing or empty" | tee -a "$LOG"
  exit 1
fi

# B) Build a Dfam-style custom rep with **THREE TAB FIELDS**: NAME<TAB>FAMILY<TAB>.
TE_FASTA_DFAM3="${CUSTOM_PREP}/human_tes.dfam3.rep"
awk '
  BEGIN{OFS=""}
  /^>/{
    n=substr($0,2); gsub(/[ \t\r]/,"",n)
    if (n=="ALU")        fam="SINE/Alu";
    else if (n=="LINE1") fam="LINE/L1";
    else if (n=="SVA")   fam="Retroposon/SVA";
    else                 fam="";
    print ">", n, "\t", fam, "\t.";   # *** THREE FIELDS, like Dfam rep ***
    next
  }
  { gsub(/[ \t\r]/,""); print toupper($0) }
' "$TE_FASTA" > "$TE_FASTA_DFAM3"

# C) Run our patched reshaper as a PRE-STEP (NO module override)
PATCH="${PATCH_DIR}/reshape_rep.py"
cat > "$PATCH" <<'PY'
#!/usr/bin/env python3
import argparse, os, re, sys, subprocess
def parse_args():
    p = argparse.ArgumentParser(description="Lightweight reshaper for custom rep + RepeatMasker .fa.out")
    p.add_argument("-rep", required=True)
    p.add_argument("-repout", required=True)
    p.add_argument("-repremove", default=None)
    p.add_argument("-pA_ME", default=None)
    p.add_argument("-outdir", required=True)
    return p.parse_args()
_hdr = re.compile(r'^>(\S+?)[#\t ]+(\S+)')
def load_rep(rep_path):
    rep = {}; name=fam=None; seq=[]
    with open(rep_path) as f:
        for line in f:
            if line.startswith('>'):
                if name and fam and seq:
                    rep[(name,fam)]=''.join(seq).replace(' ','').replace('\t','').replace('\r','').replace('\n','').upper()
                m=_hdr.match(line.strip())
                if not m:
                    name=line.strip()[1:]; fam=''
                else:
                    name,fam=m.group(1),m.group(2)
                fam=fam.split()[0]
                seq=[]
            else:
                seq.append(line.strip())
        if name and fam and seq:
            rep[(name,fam)]=''.join(seq).replace(' ','').replace('\t','').replace('\r','').replace('\n','').upper()
    return rep
def load_pairs_from_faout(path):
    pairs=set()
    with open(path) as f:
        for i,line in enumerate(f,1):
            if i<=3 or not line.strip(): continue
            cols=line.strip().split()
            if len(cols)<11: continue
            pairs.add((cols[9], cols[10]))
    return pairs
def load_family_list(path):
    if not path: return set()
    s=set()
    with open(path) as f:
        for line in f:
            t=line.strip()
            if not t or t.startswith('#'): continue
            s.add(t)
    return s
def write_reshaped(outdir, chosen):
    os.makedirs(outdir, exist_ok=True)
    outfa=os.path.join(outdir,"reshaped_repbase.fa")
    with open(outfa,"w") as w:
        for (name,fam),seq in sorted(chosen.items()):
            w.write(f">{name}\t{fam}\t.\n")
            for i in range(0,len(seq),60):
                w.write(seq[i:i+60]+"\n")
    return outfa
def run_makeblastdb(fasta, outdir):
    db=os.path.join(outdir,"repdb")
    for cmd in [
        ["makeblastdb","-in",fasta,"-out",db,"-dbtype","nucl","-parse_seqids"],
        ["/usr/local/bin/ncbi-blast-2.12.0+/bin/makeblastdb","-in",fasta,"-out",db,"-dbtype","nucl","-parse_seqids"],
    ]:
        try: subprocess.run(cmd,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE); return
        except Exception: pass
    raise SystemExit(1)
def main():
    a=parse_args()
    repmap=load_rep(a.rep)
    pairs=load_pairs_from_faout(a.repout)
    exclude=load_family_list(a.repremove)
    chosen={}; miss=[]
    for (n,f) in pairs:
        if exclude and f in exclude: continue
        k=(n,f)
        if k in repmap: chosen[k]=repmap[k]
        else:
            k2=(n.upper(),f); k3=(n.capitalize(),f)
            if k2 in repmap: chosen[k2]=repmap[k2]
            elif k3 in repmap: chosen[k3]=repmap[k3]
            else: miss.append(k)
    if not chosen:
        print("[ERROR] zero sequences matched; examples:", *[f"{n}\t{f}" for n,f in miss[:10]], sep="\n", file=sys.stderr)
        sys.exit(1)
    outfa=write_reshaped(a.outdir, chosen)
    run_makeblastdb(outfa, a.outdir)
    print(f"[OK] Wrote {outfa}", file=sys.stderr)
if __name__=="__main__": main()
PY
chmod +x "$PATCH"

# Pre-generate reshaped FASTA + BLAST DB in CUSTOM_OUT
rm -f "${CUSTOM_OUT}/reshaped_repbase.fa" "${CUSTOM_OUT}"/repdb.*
python "$PATCH" \
  -rep "${TE_FASTA_DFAM3}" \
  -repout "${FAOUT}" \
  -repremove "${CUSTOM_PREP}/non_ME_rep.txt" \
  -pA_ME "${CUSTOM_PREP}/ME_with_pA.families.txt" \
  -outdir "${CUSTOM_OUT}"

############################################
# Step 2: run MEGAnE (built-in in bg, custom in fg)
############################################
echo "[$(date)] Launching MEGAnE runs" | tee -a "$LOG"

# Built-in (hg38) — background
apptainer exec ${BIND} "$SIF" call_genotype_38 \
  -i "${MAPDIR}/${SAMPLE}.bam" \
  -fa "$REF" \
  -mk "${KMERDIR}/hg38.mk" \
  -outdir "${BUILTIN_OUT}" \
  -sample_name "${SAMPLE}" \
  -p "${MEGANE_THREADS}" &

# Custom — no overrides; -rep now passes header check
apptainer exec ${BIND} "$SIF" call_genotype \
  -i "${MAPDIR}/${SAMPLE}.bam" \
  -fa "$REF" \
  -mk "${KMERDIR}/hg38.mk" \
  -rep "${TE_FASTA_DFAM3}" \
  -repout "${FAOUT}" \
  -repremove "${CUSTOM_PREP}/non_ME_rep.txt" \
  -pA_ME "${CUSTOM_PREP}/ME_with_pA.txt" \
  -mainchr "${CUSTOM_PREP}/main_chrs.txt" \
  -fadb "${CUSTOM_PREP}/ref_blastdb" \
  -outdir "${CUSTOM_OUT}" \
  -sample_name "${SAMPLE}" \
  -p "${MEGANE_THREADS}"

wait
echo "[$(date)] DONE. Built-in: ${BUILTIN_OUT} | Custom: ${CUSTOM_OUT}" | tee -a "$LOG"
