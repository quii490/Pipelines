#!/usr/bin/env bash
set -euo pipefail

# Build the RepEnTools CHM13 reference profile from a complete local T2T FASTA
# and the matching RepeatMasker BED. Figshare's upstream downloads can be
# protected by a WAF in Pods, so the pipeline accepts audited local references
# and records their provenance instead of silently using a partial download.

FASTA=""
REPEAT_BED=""
OUT_DIR="/path/to/CUTnRUN/resources/RepEnTools/chm13v2"
THREADS=16

usage() {
  cat <<'USAGE'
Prepare a RepEnTools CHM13 reference profile.

  --fasta FILE       complete hs1/chm13v2 FASTA (plain or .gz)
  --repeat-bed FILE  RepeatMasker BED with columns chrom,start,end,repName,...
  --out-dir DIR      output profile (default: /path/to/CUTnRUN/resources/RepEnTools/chm13v2)
  --threads INT      hisat2-build threads (default: 16)
USAGE
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fasta) FASTA="$2"; shift 2 ;;
    --repeat-bed) REPEAT_BED="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage; exit 1 ;;
  esac
done
[[ -s "$FASTA" ]] || { echo "ERROR: FASTA not found: $FASTA" >&2; exit 1; }
[[ -s "$REPEAT_BED" ]] || { echo "ERROR: RepeatMasker BED not found: $REPEAT_BED" >&2; exit 1; }
if ! command -v hisat2-build >/dev/null 2>&1; then
  for candidate in /opt/conda/envs/hisat2.2.1/bin /path/to/.conda/envs/hisat2.2.1/bin; do
    if [[ -x "$candidate/hisat2-build" ]]; then PATH="$candidate:$PATH"; export PATH; break; fi
  done
fi
command -v hisat2-build >/dev/null 2>&1 || { echo "ERROR: hisat2-build is not in PATH" >&2; exit 1; }
command -v awk >/dev/null 2>&1 || { echo "ERROR: awk is required" >&2; exit 1; }
if [[ "$FASTA" == *.gz ]]; then gzip -t "$FASTA" || { echo "ERROR: corrupt gzip FASTA" >&2; exit 1; }; fi

mkdir -p "$OUT_DIR/annotation" "$OUT_DIR/indexes"
BUILD_FASTA="$FASTA"
if [[ "$FASTA" == *.gz ]]; then
  # HISAT2-build expects a plain FASTA (it does not transparently decompress
  # gzip input). Keep the decompressed copy under the reference profile so a
  # later run can reuse it without another 3-GB decompression pass.
  BUILD_FASTA="$OUT_DIR/reference/$(basename "${FASTA%.gz}")"
  mkdir -p "$(dirname "$BUILD_FASTA")"
  if [[ ! -s "$BUILD_FASTA" ]]; then gzip -cd "$FASTA" > "$BUILD_FASTA"; fi
fi
GTF="$OUT_DIR/annotation/rmsk.gtf"
awk 'BEGIN{OFS="\t"} !/^#/ && NF>=10 {
  id=$10; if (id=="" || id==".") id=$1"_"$2"_"$3"_"NR;
  gsub(/[";[:space:]]/,"_",id); gsub(/[";[:space:]]/,"_",$4);
  cls=$7; fam=$8; div=$9;
  if (cls=="") cls="Unknown"; if (fam=="") fam="Unknown";
  print $1,"RepeatMasker","exon",$2+1,$3,".",$6,".",
    "gene_id \""id"\"; transcript_id \""$4"\"; repName \""$4"\"; repClass \""cls"\"; repFamily \""fam"\"; milliDiv \""$9"\";"
}' "$REPEAT_BED" > "$GTF"
[[ -s "$GTF" ]] || { echo "ERROR: generated GTF is empty" >&2; exit 1; }

PREFIX="$OUT_DIR/indexes/chm13-2"
if [[ ! -s "${PREFIX}.1.ht2" || ! -s "${PREFIX}.8.ht2" ]]; then
  rm -f "${PREFIX}".*.ht2 "${PREFIX}".*.ht2l
  hisat2-build -p "$THREADS" "$BUILD_FASTA" "$PREFIX" > "$OUT_DIR/indexes/hisat2-build.log" 2>&1
fi
for suffix in 1 2 3 4 5 6 7 8; do
  [[ -s "${PREFIX}.${suffix}.ht2" ]] || { echo "ERROR: missing index shard ${PREFIX}.${suffix}.ht2" >&2; exit 1; }
done
cat > "$OUT_DIR/REFERENCE_METADATA.txt" <<EOF
assembly_source=$FASTA
build_fasta=$BUILD_FASTA
repeatmasker_bed_source=$REPEAT_BED
generated_gtf=$GTF
index_prefix=$PREFIX
gtf_feature=exon
gtf_group_attribute=gene_id
generated_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF
echo "[prepare_repentools_reference] ready: $OUT_DIR"
