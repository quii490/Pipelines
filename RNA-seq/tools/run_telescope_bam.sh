#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Run Telescope locus-level TE quantification from existing TE BAM files.

Required:
  --bam-dir DIR              Directory containing <sample>.te.bam
  --outdir DIR               Output directory, normally <results-dir>/12_telescope
  --te-gtf FILE              RepeatMasker TE GTF used to create the TE BAM

Optional:
  --samples CSV              Comma-separated sample names; default: every *.te.bam
  --threads INT              samtools collate threads per sample, default 4
  --jobs INT                 Concurrent samples, default 1
  --telescope-bin FILE       Telescope executable, default: telescope in PATH
  --keep-collated BOOL       Keep temporary name-collated BAMs, default false
  -h, --help                 Show this help
USAGE
}

BAM_DIR=""
OUTDIR=""
TE_GTF=""
SAMPLES=""
THREADS=4
JOBS=1
TELESCOPE_BIN="${TELESCOPE_BIN:-telescope}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
KEEP_COLLATED=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam-dir) BAM_DIR="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --te-gtf) TE_GTF="$2"; shift 2 ;;
    --samples) SAMPLES="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --telescope-bin) TELESCOPE_BIN="$2"; shift 2 ;;
    --keep-collated) KEEP_COLLATED="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -d "$BAM_DIR" ]] || { echo "--bam-dir is required and must exist" >&2; exit 1; }
[[ -n "$OUTDIR" ]] || { echo "--outdir is required" >&2; exit 1; }
[[ -f "$TE_GTF" ]] || { echo "--te-gtf is required and must exist" >&2; exit 1; }
[[ "$THREADS" =~ ^[0-9]+$ ]] && (( THREADS > 0 )) || { echo "--threads must be a positive integer" >&2; exit 1; }
[[ "$JOBS" =~ ^[0-9]+$ ]] && (( JOBS > 0 )) || { echo "--jobs must be a positive integer" >&2; exit 1; }
case "${KEEP_COLLATED,,}" in true|false) ;; *) echo "--keep-collated must be true or false" >&2; exit 1 ;; esac

if [[ "$TELESCOPE_BIN" == "telescope" ]] && ! command -v telescope >/dev/null 2>&1 && [[ -x /path/to/.conda/envs/rnaseq/bin/telescope ]]; then
  TELESCOPE_BIN=/path/to/.conda/envs/rnaseq/bin/telescope
fi
if [[ "$SAMTOOLS_BIN" == "samtools" ]] && ! command -v samtools >/dev/null 2>&1 && [[ -x /path/to/.conda/envs/rnaseq/bin/samtools ]]; then
  SAMTOOLS_BIN=/path/to/.conda/envs/rnaseq/bin/samtools
fi
if [[ "$TELESCOPE_BIN" == */* ]]; then
  [[ -x "$TELESCOPE_BIN" ]] || { echo "Telescope executable not found: $TELESCOPE_BIN" >&2; exit 1; }
else
  command -v "$TELESCOPE_BIN" >/dev/null 2>&1 || { echo "telescope not found in PATH" >&2; exit 1; }
fi
if [[ "$SAMTOOLS_BIN" == */* ]]; then
  [[ -x "$SAMTOOLS_BIN" ]] || { echo "samtools executable not found: $SAMTOOLS_BIN" >&2; exit 1; }
else
  command -v "$SAMTOOLS_BIN" >/dev/null 2>&1 || { echo "samtools not found in PATH" >&2; exit 1; }
fi

mkdir -p "$OUTDIR" "$OUTDIR/logs" "$OUTDIR/.work" "$OUTDIR/_reference"
ANNOTATION="$OUTDIR/_reference/telescope_locus_annotation.gtf"
if [[ ! -s "$ANNOTATION" ]]; then
  tmp_annotation="${ANNOTATION}.tmp.$$"
  awk '
    /locus "/ { print; next }
    {
      if (match($0, /transcript_id "[^"]+"/)) {
        locus = substr($0, RSTART + 15, RLENGTH - 16)
        print $0 " locus \"" locus "\";"
      } else {
        print
      }
    }
  ' "$TE_GTF" > "$tmp_annotation"
  mv -f "$tmp_annotation" "$ANNOTATION"
fi

declare -A requested=()
if [[ -n "$SAMPLES" ]]; then
  IFS=',' read -r -a requested_samples <<< "$SAMPLES"
  for sample in "${requested_samples[@]}"; do
    sample="${sample//[[:space:]]/}"
    [[ -n "$sample" ]] && requested["$sample"]=1
  done
fi

bams=()
while IFS= read -r -d '' bam; do
  sample="$(basename "$bam" .te.bam)"
  if (( ${#requested[@]} == 0 )) || [[ -n "${requested[$sample]:-}" ]]; then
    bams+=("$bam")
  fi
done < <(find "$BAM_DIR" -maxdepth 1 -type f -name '*.te.bam' -print0 | sort -z)

(( ${#bams[@]} > 0 )) || { echo "No matching *.te.bam files found" >&2; exit 1; }

run_one() {
  set -euo pipefail
  local bam="$1" sample work collated report report_legacy
  sample="$(basename "$bam" .te.bam)"
  work="$OUTDIR/.work/$sample"
  collated="$work/${sample}.namecollated.bam"
  report="$work/telescope_out/${sample}-telescope_report.tsv"
  report_legacy="$work/telescope_out/${sample}-telescope_report.tsv"
  if [[ -s "$OUTDIR/${sample}-telescope_report.tsv" ]]; then
    echo "[telescope] skip completed sample=$sample"
    return 0
  fi
  if [[ -s "$OUTDIR/${sample}-telescope_report.tsv" ]]; then
    shopt -s nullglob
    for f in "$OUTDIR/${sample}"-telescope_*; do
      base="$(basename "$f")"
      cp -f "$f" "$OUTDIR/${base/${sample}-telescope_/${sample}_telescope_}"
    done
    if [[ -s "$OUTDIR/${sample}-telescope_report.tsv" ]]; then
      echo "[telescope] normalized legacy outputs sample=$sample"
      return 0
    fi
  fi
  mkdir -p "$work/telescope_out"

  {
    echo "[telescope] sample=$sample"
    "$SAMTOOLS_BIN" collate -@ "$THREADS" -o "$collated" "$bam"
    "$TELESCOPE_BIN" assign \
      --attribute locus \
      --ncpu 1 \
      --outdir "$work/telescope_out" \
      --exp_tag "$sample" \
      "$collated" "$ANNOTATION"
    test -s "$report" || test -s "$report_legacy"
    shopt -s nullglob
    telescope_files=("$work/telescope_out/${sample}"_telescope_* "$work/telescope_out/${sample}"-telescope_*)
    if (( ${#telescope_files[@]} == 0 )); then
      echo "[telescope] ERROR: no Telescope output files found for sample=$sample" >&2
      exit 1
    fi
    for f in "${telescope_files[@]}"; do
      base="$(basename "$f")"
      cp -f "$f" "$OUTDIR/${base/${sample}-telescope_/${sample}_telescope_}"
    done
    [[ "$KEEP_COLLATED" == true ]] || rm -f "$collated"
    echo "[telescope] completed sample=$sample"
  } 2>&1 | tee -a "$OUTDIR/logs/${sample}.telescope.log"
}

export BAM_DIR OUTDIR ANNOTATION THREADS TELESCOPE_BIN SAMTOOLS_BIN KEEP_COLLATED
export -f run_one
printf '%s\0' "${bams[@]}" | xargs -0 -r -n 1 -P "$JOBS" bash -c 'run_one "$1"' _

printf '[telescope] reports=%s\n' "$(find "$OUTDIR" -maxdepth 1 -type f -name '*-telescope_report.tsv' | wc -l | tr -d ' ')"
