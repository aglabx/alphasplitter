#!/usr/bin/env bash
# =============================================================================
# satkmer_compare.sh  --  satellite k-mer comparison of two BAM files
#
# Given two BAM files (e.g. a matched normal/tumour pair, or any two samples),
# this:
#   1. counts canonical k-mers directly from each BAM with KMC          (no FASTQ unpack)
#   2. builds a coverage-robust RANK comparison of the high-copy k-mers
#   3. (optional, with a reference pack) computes a per-chromosome
#      "satellite karyotype": for each chromosome, the relative abundance of
#      its centromeric alpha-satellite array-specific marker k-mers, and the
#      sample2 / sample1 ratio + rank shift  -> picks up whole-chromosome
#      gains/losses (e.g. the KICH chromosome-loss signature) from satellite
#      k-mers alone.
#
# Why ranks / fractions, not raw counts?  Two libraries usually differ in
# depth and prep; raw k-mer counts are not comparable.  Rank order and
# normalised fractions (a chromosome's marker-sum / total) are coverage-robust.
#
# Requirements: kmc, kmc_tools (>=3.x), awk, sort, comm  (samtools NOT needed -
# KMC reads BAM natively).  blastn only needed for the optional alpha-extract.
#
# Reference pack (for step 3): a directory containing
#   marker_chr_map.tsv          # cols: <canonical_23mer> <field_idx> <array_id> <chr> [...]
#   canonical_alpha.fa          # one canonical alpha-satellite monomer  (optional)
#   clean_alpha_k23_DF*.txt     # clean alpha 23-mer sets, fwd+RC, one per line  (optional)
# Download from: <fill in your hosting URL>  or build with build_refpack.sh
# =============================================================================
set -euo pipefail

# ---------- defaults ----------
K=23
CI=3                 # min k-mer count to keep (drops most sequencing errors)
HICOPY=50000         # for the rank comparison: only consider k-mers with count >= this
TOPN=2000            # union of top-N high-copy k-mers from each sample to tabulate
THREADS=$(command -v nproc >/dev/null && nproc || echo 8)
OUTDIR=satkmer_out
NAMES=""
REFPACK=""
TMPDIR_BASE="${TMPDIR:-/tmp}"

usage() {
  cat <<EOF
Usage: $0 [options] <sample1.bam> <sample2.bam>

  -o DIR        output directory                       (default: $OUTDIR)
  -t N          threads                                (default: $THREADS)
  -k N          k-mer size                             (default: $K)
  -c N          min k-mer count (KMC -ci)              (default: $CI)
  --hicopy N    high-copy threshold for rank compare   (default: $HICOPY)
  --topn N      union of top-N per sample to tabulate  (default: $TOPN)
  --names A,B   sample names (default: BAM basenames)
  --refpack DIR reference pack dir (enables per-chromosome satellite karyotype)
  -h            this help

Outputs in DIR/:
  <name>_k<k>.kmc_pre / .kmc_suf      KMC databases (kept; reusable)
  <name>_hicopy_ranked.tsv            high-copy k-mers, sorted by count, ranked
  rank_compare.tsv                    union top-N: kmer, count1,rank1, count2,rank2, rank_delta
  per_chromosome_satellite_cn.tsv     (if --refpack) chr, n_markers, sum1,sum2, frac1,frac2, ratio(2/1), rank1,rank2,rank_delta
  SUMMARY.md                          human-readable summary
EOF
  exit "${1:-0}"
}

# ---------- parse args ----------
POS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o) OUTDIR="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    -k) K="$2"; shift 2;;
    -c) CI="$2"; shift 2;;
    --hicopy) HICOPY="$2"; shift 2;;
    --topn) TOPN="$2"; shift 2;;
    --names) NAMES="$2"; shift 2;;
    --refpack) REFPACK="$2"; shift 2;;
    -h|--help) usage 0;;
    -*) echo "unknown option: $1" >&2; usage 1;;
    *) POS+=("$1"); shift;;
  esac
done
[[ ${#POS[@]} -eq 2 ]] || { echo "ERROR: need exactly two BAM files" >&2; usage 1; }
BAM1="${POS[0]}"; BAM2="${POS[1]}"
for b in "$BAM1" "$BAM2"; do [[ -s "$b" ]] || { echo "ERROR: not found / empty: $b" >&2; exit 1; }; done
for t in kmc kmc_tools awk sort comm; do command -v "$t" >/dev/null || { echo "ERROR: '$t' not in PATH" >&2; exit 1; }; done

if [[ -n "$NAMES" ]]; then N1="${NAMES%%,*}"; N2="${NAMES##*,}"; else
  N1=$(basename "$BAM1"); N1="${N1%.bam}"; N2=$(basename "$BAM2"); N2="${N2%.bam}"
fi
[[ "$N1" != "$N2" ]] || { N1="${N1}_1"; N2="${N2}_2"; }

mkdir -p "$OUTDIR"
LOG="$OUTDIR/run.log"
log(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG" >&2; }
log "satkmer_compare  k=$K ci=$CI threads=$THREADS"
log "  sample1: $N1  <- $BAM1"
log "  sample2: $N2  <- $BAM2"
[[ -n "$REFPACK" ]] && log "  refpack: $REFPACK"

# ---------- step 1: KMC ----------
run_kmc(){  # name bam -> $OUTDIR/<name>_k<k>.kmc_{pre,suf}
  local name="$1" bam="$2" db="$OUTDIR/${1}_k${K}"
  if [[ -s "${db}.kmc_pre" && -s "${db}.kmc_suf" ]]; then log "  KMC $name: DB exists, skipping"; return; fi
  local tmp; tmp=$(mktemp -d "${TMPDIR_BASE}/kmc.${name}.XXXXXX")
  log "  KMC $name ..."
  kmc -k"$K" -ci"$CI" -cs1000000000 -t"$THREADS" -fbam "$bam" "$db" "$tmp" >>"$LOG" 2>&1
  rm -rf "$tmp"
  log "  KMC $name done: $(stat -c %s "${db}.kmc_suf" 2>/dev/null || echo ?) B suffix"
}
run_kmc "$N1" "$BAM1"
run_kmc "$N2" "$BAM2"
DB1="$OUTDIR/${N1}_k${K}"; DB2="$OUTDIR/${N2}_k${K}"

# ---------- step 2: rank comparison of high-copy k-mers ----------
log "step 2: rank comparison (high-copy k-mers, count >= $HICOPY)"
for nm in "$N1" "$N2"; do
  db="$OUTDIR/${nm}_k${K}"
  kmc_tools -t"$THREADS" transform "$db" -ci"$HICOPY" dump "$OUTDIR/${nm}_hicopy.txt" >>"$LOG" 2>&1
  sort -t$'\t' -k2,2 -rn "$OUTDIR/${nm}_hicopy.txt" | awk -F'\t' '{print $1"\t"NR"\t"$2}' > "$OUTDIR/${nm}_hicopy_ranked.tsv"  # kmer rank count
  log "  $nm: $(wc -l < "$OUTDIR/${nm}_hicopy_ranked.tsv") high-copy k-mers"
done
head -"$TOPN" "$OUTDIR/${N1}_hicopy_ranked.tsv" | cut -f1  > "$OUTDIR/_top.ids"
head -"$TOPN" "$OUTDIR/${N2}_hicopy_ranked.tsv" | cut -f1 >> "$OUTDIR/_top.ids"
sort -u "$OUTDIR/_top.ids" > "$OUTDIR/_top_union.ids"
awk -F'\t' -v N1="$N1" -v N2="$N2" '
  FILENAME ~ /'"${N1}"'_hicopy_ranked/ { r1[$1]=$2; c1[$1]=$3; next }
  FILENAME ~ /'"${N2}"'_hicopy_ranked/ { r2[$1]=$2; c2[$1]=$3; next }
  { k=$1; cc1=(k in c1)?c1[k]:0; rr1=(k in r1)?r1[k]:"NA"; cc2=(k in c2)?c2[k]:0; rr2=(k in r2)?r2[k]:"NA";
    d=(rr1!="NA"&&rr2!="NA")?(rr2-rr1):"NA";
    print k"\t"cc1"\t"rr1"\t"cc2"\t"rr2"\t"d }
' "$OUTDIR/${N1}_hicopy_ranked.tsv" "$OUTDIR/${N2}_hicopy_ranked.tsv" "$OUTDIR/_top_union.ids" \
  | sort -t$'\t' -k3,3n > "$OUTDIR/_rc.tsv"
{ echo -e "kmer\tcount_${N1}\trank_${N1}\tcount_${N2}\trank_${N2}\trank_delta(${N2}-${N1})"; cat "$OUTDIR/_rc.tsv"; } > "$OUTDIR/rank_compare.tsv"
rm -f "$OUTDIR/_top.ids" "$OUTDIR/_top_union.ids" "$OUTDIR/_rc.tsv" "$OUTDIR/${N1}_hicopy.txt" "$OUTDIR/${N2}_hicopy.txt"
log "  -> $OUTDIR/rank_compare.tsv ($(($(wc -l < "$OUTDIR/rank_compare.tsv")-1)) k-mers)"

# ---------- step 3 (optional): per-chromosome satellite karyotype ----------
if [[ -n "$REFPACK" ]]; then
  MAP="$REFPACK/marker_chr_map.tsv"
  [[ -s "$MAP" ]] || { echo "ERROR: $MAP not found in refpack" >&2; exit 1; }
  log "step 3: per-chromosome satellite karyotype (markers: $MAP)"
  # build / reuse a KMC DB of the marker k-mers (each as its own k-bp record, count 1)
  MDB="$REFPACK/markers_k${K}"
  if [[ ! -s "${MDB}.kmc_pre" ]]; then
    log "  building markers KMC DB (first run) ..."
    tmp=$(mktemp -d "${TMPDIR_BASE}/kmc.markers.XXXXXX")
    awk -F'\t' '{print ">m"NR"\n"$1}' "$MAP" > "$REFPACK/_markers.fa"
    kmc -k"$K" -ci1 -cs2 -fm -t"$THREADS" "$REFPACK/_markers.fa" "$MDB" "$tmp" >>"$LOG" 2>&1
    rm -rf "$tmp" "$REFPACK/_markers.fa"
  fi
  # intersect each sample DB with the markers DB, keep the sample (left) counts
  for nm in "$N1" "$N2"; do
    tmp=$(mktemp -d "${TMPDIR_BASE}/kmc.isect.XXXXXX")
    kmc_tools -t"$THREADS" simple "$OUTDIR/${nm}_k${K}" -ci"$CI" "$MDB" -ci1 intersect "$OUTDIR/${nm}_markers" -ocleft >>"$LOG" 2>&1
    kmc_tools -t"$THREADS" transform "$OUTDIR/${nm}_markers" dump "$OUTDIR/${nm}_markers.txt" >>"$LOG" 2>&1
    rm -rf "$tmp"; rm -f "$OUTDIR/${nm}_markers.kmc_pre" "$OUTDIR/${nm}_markers.kmc_suf"
    log "  $nm: $(wc -l < "$OUTDIR/${nm}_markers.txt") marker k-mers detected"
  done
  # join: kmer -> chr (col4 of map) + count1 + count2
  LC_ALL=C sort -t$'\t' -k1,1 -u <(cut -f1,4 "$MAP") > "$OUTDIR/_mc"
  LC_ALL=C sort -t$'\t' -k1,1 "$OUTDIR/${N1}_markers.txt" > "$OUTDIR/_m1"
  LC_ALL=C sort -t$'\t' -k1,1 "$OUTDIR/${N2}_markers.txt" > "$OUTDIR/_m2"
  LC_ALL=C join -t$'\t' -a1 -e0 -o 1.1,1.2,2.2 "$OUTDIR/_mc" "$OUTDIR/_m1" > "$OUTDIR/_j1"; LC_ALL=C sort -t$'\t' -k1,1 "$OUTDIR/_j1" > "$OUTDIR/_j1s"
  LC_ALL=C join -t$'\t' -a1 -e0 -o 1.1,1.2,1.3,2.2 "$OUTDIR/_j1s" "$OUTDIR/_m2" > "$OUTDIR/marker_counts.tsv"   # kmer chr count1 count2
  # aggregate per chr + global totals + ranks
  awk -F'\t' '
    { c=$2; n[c]++; s1[c]+=$3; s2[c]+=$4; if($3>0)d1[c]++; if($4>0)d2[c]++; gt1+=$3; gt2+=$4 }
    END{
      m=0; for(c in n){ m++; chr[m]=c; nn[m]=n[c]; dd1[m]=d1[c]; dd2[m]=d2[c]; ss1[m]=s1[c]; ss2[m]=s2[c] }
      for(i=1;i<=m;i++){o1[i]=i;o2[i]=i}
      for(i=1;i<m;i++)for(j=1;j<=m-i;j++){if(ss1[o1[j]]<ss1[o1[j+1]]){t=o1[j];o1[j]=o1[j+1];o1[j+1]=t}}
      for(i=1;i<m;i++)for(j=1;j<=m-i;j++){if(ss2[o2[j]]<ss2[o2[j+1]]){t=o2[j];o2[j]=o2[j+1];o2[j+1]=t}}
      for(i=1;i<=m;i++){rk1[o1[i]]=i; rk2[o2[i]]=i}
      print "chr\tn_markers\tdet_'"${N1}"'\tdet_'"${N2}"'\tsum_'"${N1}"'\tsum_'"${N2}"'\tfrac_'"${N1}"'(%)\tfrac_'"${N2}"'(%)\tratio_'"${N2}"'/'"${N1}"'\trank_'"${N1}"'\trank_'"${N2}"'\trank_delta"
      for(i=1;i<=m;i++){k=o1[i]; f1=ss1[k]/gt1; f2=ss2[k]/gt2;
        printf "%s\t%d\t%d\t%d\t%.0f\t%.0f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\n", chr[k],nn[k],dd1[k],dd2[k],ss1[k],ss2[k],100*f1,100*f2,(f1>0?f2/f1:0),rk1[k],rk2[k],rk2[k]-rk1[k] }
    }' "$OUTDIR/marker_counts.tsv" > "$OUTDIR/per_chromosome_satellite_cn.tsv"
  rm -f "$OUTDIR/_mc" "$OUTDIR/_m1" "$OUTDIR/_m2" "$OUTDIR/_j1" "$OUTDIR/_j1s" "$OUTDIR/${N1}_markers.txt" "$OUTDIR/${N2}_markers.txt"
  log "  -> $OUTDIR/per_chromosome_satellite_cn.tsv"
fi

# ---------- summary ----------
{
  echo "# satkmer_compare summary"
  echo
  echo "- sample1 (\`$N1\`): \`$BAM1\`"
  echo "- sample2 (\`$N2\`): \`$BAM2\`"
  echo "- k = $K, min count (KMC -ci) = $CI"
  echo
  echo "## KMC k-mer databases"
  for nm in "$N1" "$N2"; do
    pre="$OUTDIR/${nm}_k${K}.kmc_pre"
    [[ -s "$pre" ]] && echo "- \`${nm}_k${K}.kmc_{pre,suf}\` ($(du -Lch "$OUTDIR/${nm}_k${K}".kmc_* 2>/dev/null | tail -1 | cut -f1))"
  done
  echo
  echo "## Rank comparison (\`rank_compare.tsv\`)"
  echo "Top high-copy k-mers (count >= $HICOPY), ranked by count in each sample; union of top-$TOPN."
  echo "\`rank_delta = rank_${N2} - rank_${N1}\` — negative = relatively UP in $N2, positive = relatively DOWN in $N2."
  echo
  echo "Largest rank shifts:"
  echo '```'
  tail -n +2 "$OUTDIR/rank_compare.tsv" | awk -F'\t' '$6!="NA"{d=$6; if(d<0)d=-d; print d"\t"$0}' | sort -rn | head -15 | cut -f2- | column -t
  echo '```'
  if [[ -s "$OUTDIR/per_chromosome_satellite_cn.tsv" ]]; then
    echo
    echo "## Per-chromosome satellite karyotype (\`per_chromosome_satellite_cn.tsv\`)"
    echo "Relative abundance of each chromosome's centromeric alpha-satellite marker k-mers."
    echo "\`ratio = frac_${N2} / frac_${N1}\` — < 1: chr relatively lost in $N2; > 1: relatively gained."
    echo
    echo '```'
    column -t "$OUTDIR/per_chromosome_satellite_cn.tsv"
    echo '```'
    echo
    echo "Chromosomes with ratio < 0.85 (candidate losses in $N2):"
    echo '```'
    { head -1 "$OUTDIR/per_chromosome_satellite_cn.tsv"; tail -n +2 "$OUTDIR/per_chromosome_satellite_cn.tsv" | awk -F'\t' '$9<0.85'; } | column -t
    echo '```'
  fi
  echo
  echo "_Caveats: ranks/fractions are coverage-robust but cannot be turned into absolute copy numbers without a normalising assumption; tumour purity < 100% dilutes losses; CHM13-derived markers approximate the patient; satellite-region signal ≈ centromeric copy number, not necessarily whole-arm._"
} > "$OUTDIR/SUMMARY.md"
log "DONE. See $OUTDIR/SUMMARY.md"
cat "$OUTDIR/SUMMARY.md"
