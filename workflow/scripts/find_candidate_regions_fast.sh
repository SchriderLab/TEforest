#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# usage / arg‑parse helpers
###############################################################################
usage() {
    echo "Usage: $0 -o <output_dir> -c <consensus.fasta> -r <ref_te_locations.bed> \\
-g <reference_genome.fasta> -@ <threads> -1 <R1.fastq.gz> -2 <R2.fastq.gz> \\
[-e] -n <genome_name>

  -e   map only to TE consensus sequences (skip genomic TE copies)
"
    exit 1
}
[[ $# -eq 0 ]] && usage

# ──────────────── NEW: add “e” to the optstring ────────────────
while getopts o:c:@:g:1:2:r:n:e flag; do
    case $flag in
        o) output_directory=${OPTARG} ;;
        c) fasta_file=${OPTARG} ;;
        @) threads=${OPTARG} ;;
        g) ref_genome=${OPTARG} ;;
        1) fq1=${OPTARG} ;;
        2) fq2=${OPTARG} ;;
        r) ref_te_locations=${OPTARG} ;;
        n) genome_name=${OPTARG} ;;
        e) consensus_only=1 ;;                     # NEW
        *) usage ;;
    esac
done
: "${consensus_only:=0}"                          # default = OFF
[[ ! $threads =~ ^[0-9]+$ ]] && { echo "ERROR: –@ must be an integer"; exit 1; }

###############################################################################
# sanity checks & I/O prep
###############################################################################
mkdir -p "$output_directory"
# BED file is optional when ‑e is set
for f in "$fasta_file" "$ref_genome" "$fq1" "$fq2"; do
    [[ -f $f ]] || { echo "ERROR: $f not found"; exit 1; }
done
if (( ! consensus_only )); then
    [[ -f $ref_te_locations ]] || { echo "ERROR: $ref_te_locations not found"; exit 1; }
fi

###############################################################################
# 0. Ensure the reference genome FASTA indexes cleanly (faidx)
###############################################################################
fix_and_index_fasta() {
    local fasta=$1
    local tmp=${fasta}.tmp

    echo "[INFO] Cleaning and re‑wrapping $fasta"
    awk '
        BEGIN { ORS=""; }
        /^>/ {
            if (seq) {
                while (length(seq) > 60) { print substr(seq,1,60) "\n"; seq = substr(seq,61) }
                print seq "\n"; seq = ""
            }
            print $0 "\n"; next
        }
        {
            gsub(/[^A-Za-z]/,""); seq = seq toupper($0)
        }
        END {
            if (seq) {
                while (length(seq) > 60) { print substr(seq,1,60) "\n"; seq = substr(seq,61) }
                print seq "\n"
            }
        }
    ' "$fasta" > "$tmp" && mv "$tmp" "$fasta"
    samtools faidx "$fasta" || { echo "ERROR: cannot index $fasta"; exit 1; }
}
if [ ! -f "${ref_genome}.fai" ]; then
    echo "[INFO] Indexing $ref_genome"
    samtools faidx "$ref_genome" 2>/dev/null || fix_and_index_fasta "$ref_genome"
fi

###############################################################################
# 1. Build TE reference (consensus + genomic TE copies)   <── MODIFIED
###############################################################################
if (( consensus_only )); then
    echo "[INFO] ‑e set: using consensus FASTA only"
    cp "$fasta_file" "$output_directory/te.tmp.fasta"
    sed 's/-/_/g' "$output_directory/te.tmp.fasta" > "$output_directory/te.fasta"
    rm -f "$output_directory/te.tmp.fasta"
else
    echo "[INFO] building consensus + genomic‑copy reference"
    awk '
        NF >= 3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ { OFS="\t"; print $1,$2,$3,$4,$5,$6 }
    ' "$ref_te_locations" |
    bedtools getfasta -fi "$ref_genome" -bed - -fo "$output_directory/extracted_sequences.fasta" -nameOnly

    cat "$output_directory/extracted_sequences.fasta" "$fasta_file" > "$output_directory/te.tmp.fasta"
    sed 's/-/_/g' "$output_directory/te.tmp.fasta" > "$output_directory/te.fasta"
    #rm "$output_directory"/{extracted_sequences.fasta,te.tmp.fasta}
fi

###############################################################################
# 2. Initial map: reads → TE reference
###############################################################################
bwa-mem2 index "$output_directory/te.fasta"
bwa-mem2 mem -t "$threads" "$output_directory/te.fasta" "$fq1" "$fq2" |
    samtools view -@ "$threads" -b - |
    samtools sort -@ "$threads" -o "$output_directory/${genome_name}_to_alltes.bam" -
samtools index "$output_directory/${genome_name}_to_alltes.bam"

###############################################################################
# 3. Collect list of TE families present in the consensus FASTA
###############################################################################
mapfile -t names < <(grep '^>' "$fasta_file" | cut -c2-)

###############################################################################
# 4. Automatic thread allocation
###############################################################################
num_fam=${#names[@]}
if (( num_fam <= threads )); then
    family_threads=$(( threads / num_fam ))
    jobs=$num_fam
else
    family_threads=1
    jobs=$threads
fi
(( family_threads == 0 )) && family_threads=1
echo "Running $jobs families at a time, $family_threads thread(s) each"

###############################################################################
# 5. Per‑family function (runs in GNU Parallel) — WITH DEBUG OUTPUT
###############################################################################
align_family() {
    local name=$1
    local tmp
    tmp=$(mktemp -d -p "$output_directory" "${name}.XXXX")

    echo "=== [$name] starting in $tmp ===" >&2

    ###########################################################################
    # a)  build contig list for this family
    ###########################################################################
    mapfile -t contigs < <(
        samtools view -H "$output_directory/${genome_name}_to_alltes.bam" |
        awk -v fam="$name" -F'\t' '
            /^@SQ/ {
                split($2,a,":"); gsub(/-/, "_", a[2]);
                if ( a[2] ~ ("^" fam "($|_)") || a[2] ~ ("(_|^)" fam "$") )
                    print a[2];
            }
        '
    )
    echo "[$name] contig count: ${#contigs[@]}  (${contigs[*]})" >&2

    if (( ${#contigs[@]} == 0 )); then
        echo "[$name] no contigs → skipping" >&2
        #rm -rf "$tmp"; return
    fi

    ###########################################################################
    # b)  Build a BED of whole-contig intervals for this family (samtools 1.6)
    ###########################################################################
    samtools view -H "$output_directory/${genome_name}_to_alltes.bam" |
    awk -v fam="$name" -F'\t' '
        $1=="@SQ" {
            sn=""; ln=""
            for (i=1; i<=NF; i++) {
                if ($i ~ /^SN:/) sn=substr($i,4)
                else if ($i ~ /^LN:/) ln=substr($i,4)
            }
            # Match contigs whose name starts/ends with the family token
            if (sn ~ ("^" fam "($|_)") || sn ~ ("(_|^)" fam "$")) {
                # BED is 0-based, half-open: [0, ln)
                print sn "\t0\t" ln
            }
        }
    ' > "$tmp/${name}.bed"

    # Sanity: if no regions, bail early
    if [ ! -s "$tmp/${name}.bed" ]; then
        echo "[$name] no regions in BED → skipping" >&2
        rm -rf "$tmp"
        return
    fi

    # Filter using regions (requires 3-column BED in samtools 1.6)
    samtools view -@ "$family_threads" -b -L "$tmp/${name}.bed" \
        "$output_directory/${genome_name}_to_alltes.bam" \
        > "$tmp/${name}.hits.bam"

    echo "[$name] hits.bam reads: $(samtools view -c "$tmp/${name}.hits.bam")" >&2

    ###########################################################################
    # c)  unique QNAMEs from hits
    ###########################################################################
    samtools view -@ "$family_threads" -h "$tmp/${name}.hits.bam" |
        awk '$1 !~ /^@/ {print $1}' |
        sort -u > "$tmp/${name}.qnames"

    echo "[$name] qname count:  $(wc -l < "$tmp/${name}.qnames")" >&2

    ###########################################################################
    # d)  fetch both mates by QNAME (works on any samtools version)
    ###########################################################################
    #   stream all‑TE BAM through awk; keep header lines and any record whose
    #   QNAME appears in the name set built above
    samtools view -@ "$family_threads" -h \
        "$output_directory/${genome_name}_to_alltes.bam" |
    awk '
        BEGIN { while ((getline < "'"$tmp/${name}.qnames"'") > 0) seen[$1]=1 }
        /^@/ { print; next }
        ($1 in seen)
    ' |
    samtools view -@ "$family_threads" -b -o "$tmp/${name}.pair.bam" -

    echo "[$name] pair.bam reads: $(samtools view -c "$tmp/${name}.pair.bam")" >&2

    ###########################################################################
    # e)  bam → paired FASTQ
    ###########################################################################
    samtools sort -n -@ "$family_threads" -o "$tmp/${name}.name.bam" \
        "$tmp/${name}.pair.bam"

    samtools fastq -@ "$family_threads" \
        -1 "$tmp/${name}_1.fq" -2 "$tmp/${name}_2.fq" \
        -0 /dev/null -s /dev/null -n "$tmp/${name}.name.bam"

    echo "[$name] FASTQ1 lines: $(wc -l < "$tmp/${name}_1.fq")" >&2
    echo "[$name] FASTQ2 lines: $(wc -l < "$tmp/${name}_2.fq")" >&2

    ###########################################################################
    # f)  remap to reference genome
    ###########################################################################
    bwa-mem2 mem -t "$family_threads" "$ref_genome" \
        "$tmp/${name}_1.fq" "$tmp/${name}_2.fq" |
        samtools view -@ "$family_threads" -b - |
        samtools sort -@ "$family_threads" -o "$output_directory/${name}_firstmap.bam" -

    samtools index "$output_directory/${name}_firstmap.bam"
    ###############################################################################
    # g) keep discordant / orphaned pairs  OR  split-read evidence
    ###############################################################################
    samtools view -h "$output_directory/${name}_firstmap.bam" | \
    awk '
        BEGIN { FS = OFS = "\t"; kept = 0; removed = 0 }

        # ---- helpers ------------------------------------------------------------
        # True if <bit> is set in <val>   (works in every awk – no "&" operator)
        function hasBit(val, bit) { return int(val / bit) % 2 }

        # Does the CIGAR contain soft or hard clips?
        function hasClip(cigar)  { return cigar ~ /[0-9][SH]/ }

        # -------------------------------------------------------------------------
        /^@/ { print; next }                           # pass header lines

        {
            flag = $2 + 0

            properPair   = hasBit(flag,   2)           # 0x002
            mateUnmapped = hasBit(flag,   8)           # 0x008
            supplementary= hasBit(flag,2048)           # 0x800

            isSplit   = hasClip($6) || supplementary || index($0,"SA:Z:")
            discord   = !properPair                   # not a proper pair

            if ( discord || mateUnmapped || isSplit ) {
                print
                ++kept
            } else ++removed
        }

        END {
            printf "## filter-stats: kept %d; removed %d\n", kept, removed > "/dev/stderr"
        }
    ' | samtools view -b -o "$output_directory/${name}_to_ISO1.bam" -

    samtools index "$output_directory/${name}_to_ISO1.bam"

    # optional: comment-out next line if you still want the unfiltered BAM around
    #rm "$output_directory/${name}_to_ISO1.bam"*
    echo "[$name] final BAM reads: $(samtools view -c "$output_directory/${name}_to_ISO1.bam")" >&2
    echo "=== [$name] done ===" >&2

    rm -rf "$tmp"
}

export -f align_family
export output_directory ref_genome genome_name family_threads

###############################################################################
# 6. Launch parallel run
###############################################################################
printf '%s\n' "${names[@]}" | parallel -j "$jobs" align_family {}

###############################################################################
# 7. Clean global intermediates
###############################################################################
#rm "$output_directory/${genome_name}_to_alltes.bam"*
#rm "$output_directory/te.fasta"*          # bwa‑mem2 adds .{amb,ann,bwt,pac,sa}
touch ${output_directory}/completed.txt

echo "All done ✔  (outputs: *_to_ISO1.bam + .bai for each family)"
