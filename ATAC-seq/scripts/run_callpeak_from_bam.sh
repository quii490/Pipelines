#!/usr/bin/env bash
set -euo pipefail

resolve_script_dir() {
  local src="${BASH_SOURCE[0]}"
  while [ -L "$src" ]; do
    local dir
    dir="$(cd -P "$(dirname "$src")" && pwd)"
    src="$(readlink "$src")"
    [[ "$src" != /* ]] && src="$dir/$src"
  done
  cd -P "$(dirname "$src")" && pwd
}

SCRIPT_DIR="$(resolve_script_dir)"
# shellcheck source=/dev/null
source "$SCRIPT_DIR/species_config_lib.sh"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="$ROOT_DIR/conf/species_refs.config"

usage() {
  cat <<'USAGE'
Usage:
  run_callpeak_from_bam.sh --bam-dir DIR --species hg38 --outdir DIR [options]
  run_callpeak_from_bam.sh --bam-glob "bam/*.clean.bam" --species hg38 --outdir DIR [options]
  run_callpeak_from_bam.sh --bam sample.clean.bam --species hg38 --outdir DIR [options]

Purpose:
  ATAC-seq peak calling from cleaned BAM files with MACS3.
  It auto-detects paired-end BAMs and uses BAMPE; single-end BAMs use the common
  ATAC shift/extsize model.

Options:
  --species STR                  hg38 | mm10 | mm39. Default: hg38.
  --genome-size STR              MACS3 -g value; default from conf/species_refs.config.
  --blacklist FILE               Optional blacklist BED; default from config when available.
  --format STR                   auto | BAMPE | BAM. Default: auto.
  --qvalue FLOAT                 MACS3 -q cutoff. Default: 0.05.
  --pvalue FLOAT                 MACS3 -p cutoff. If set, overrides --qvalue.
  --shift INT                    Single-end ATAC shift. Default: -100.
  --extsize INT                  Single-end ATAC extsize. Default: 200.
  --keep-dup STR                 MACS3 --keep-dup. Default: all.
  --broad                        Call broad peaks instead of narrow peaks.
  --no-summits                   Do not pass --call-summits.
  --cores INT                    Default: SLURM_CPUS_PER_TASK or 8.

Outputs:
  05_peaks/raw_macs3             Raw MACS3 outputs.
  05_peaks/filtered              Blacklist-filtered peaks when blacklist is available.
  05_peaks/summits               Standardized summit BED files.
  05_peaks/qc/frip.tsv           FRiP summary when bedtools/samtools are available.
USAGE
}

bam_dir=""
bam_glob=""
single_bam=""
species="hg38"
outdir=""
genome_size=""
blacklist=""
format="auto"
qvalue="0.05"
pvalue=""
shift_size="-100"
extsize="200"
keep_dup="all"
broad="false"
call_summits="true"
cores="${SLURM_CPUS_PER_TASK:-8}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam-dir) bam_dir="$2"; shift 2 ;;
    --bam-glob) bam_glob="$2"; shift 2 ;;
    --bam) single_bam="$2"; shift 2 ;;
    --species) species="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    --genome-size) genome_size="$2"; shift 2 ;;
    --blacklist) blacklist="$2"; shift 2 ;;
    --format) format="$2"; shift 2 ;;
    --qvalue|-q) qvalue="$2"; shift 2 ;;
    --pvalue|-p) pvalue="$2"; shift 2 ;;
    --shift) shift_size="$2"; shift 2 ;;
    --extsize) extsize="$2"; shift 2 ;;
    --keep-dup) keep_dup="$2"; shift 2 ;;
    --broad) broad="true"; shift ;;
    --no-summits) call_summits="false"; shift ;;
    --cores|--threads) cores="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "$outdir" ]] || { usage >&2; exit 1; }
[[ -f "$CONFIG_FILE" ]] || { echo "ERROR: config not found: $CONFIG_FILE" >&2; exit 1; }
[[ -n "$genome_size" ]] || genome_size="$(get_species_param "$species" genome_size "$CONFIG_FILE" || true)"
[[ -n "$blacklist" ]] || blacklist="$(get_species_param "$species" blacklist "$CONFIG_FILE" || true)"
[[ -n "$genome_size" ]] || { echo "ERROR: genome_size missing for species $species" >&2; exit 1; }

command -v macs3 >/dev/null || { echo "ERROR: macs3 not found in PATH" >&2; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not found in PATH" >&2; exit 1; }

shopt -s nullglob
if [[ -n "$single_bam" ]]; then
  bam_files=( "$single_bam" )
elif [[ -n "$bam_glob" ]]; then
  bam_files=( $bam_glob )
else
  [[ -d "$bam_dir" ]] || { echo "ERROR: bam dir not found: $bam_dir" >&2; exit 1; }
  bam_files=( "$bam_dir"/*.bam )
fi
[[ ${#bam_files[@]} -gt 0 ]] || { echo "ERROR: no BAM files found" >&2; exit 1; }

raw_dir="$outdir/05_peaks/raw_macs3"
filtered_dir="$outdir/05_peaks/filtered"
summit_dir="$outdir/05_peaks/summits"
qc_dir="$outdir/05_peaks/qc"
mkdir -p "$raw_dir" "$filtered_dir" "$summit_dir" "$qc_dir"
frip_tsv="$qc_dir/frip.tsv"
printf "sample\tbam\tpeak_file\ttotal_reads\tin_peak_reads\tfrip\tmacs_format\n" > "$frip_tsv"

for bam in "${bam_files[@]}"; do
  [[ -s "$bam" ]] || { echo "ERROR: BAM not found: $bam" >&2; exit 1; }
  sample="$(basename "$bam")"
  sample="${sample%.clean.bam}"
  sample="${sample%.bam}"

  if [[ ! -s "${bam}.bai" && ! -s "${bam%.bam}.bai" ]]; then
    samtools index -@ "$cores" "$bam"
  fi

  macs_format="$format"
  if [[ "$format" == "auto" ]]; then
    paired_n="$(samtools view -c -f 1 "$bam")"
    if [[ "$paired_n" -gt 0 ]]; then
      macs_format="BAMPE"
    else
      macs_format="BAM"
    fi
  fi

  macs_args=(callpeak -t "$bam" -f "$macs_format" -g "$genome_size" -n "$sample"
    --nomodel --keep-dup "$keep_dup" --outdir "$raw_dir")
  if [[ -n "$pvalue" ]]; then
    macs_args+=(-p "$pvalue")
  else
    macs_args+=(-q "$qvalue")
  fi
  [[ "$broad" == "true" ]] && macs_args+=(--broad)
  [[ "$call_summits" == "true" && "$broad" != "true" ]] && macs_args+=(--call-summits)
  if [[ "$macs_format" == "BAM" ]]; then
    macs_args+=(--shift "$shift_size" --extsize "$extsize")
  fi

  echo "[run_callpeak_from_bam] MACS3 sample=$sample format=$macs_format"
  macs3 "${macs_args[@]}"

  peak_file="$raw_dir/${sample}_peaks.narrowPeak"
  [[ "$broad" == "true" ]] && peak_file="$raw_dir/${sample}_peaks.broadPeak"
  [[ -s "$peak_file" ]] || { echo "ERROR: peak file missing: $peak_file" >&2; exit 1; }

  final_peak="$peak_file"
  if [[ -n "$blacklist" && -s "$blacklist" && "$(command -v bedtools || true)" ]]; then
    final_peak="$filtered_dir/${sample}.peaks.${broad/true/broadPeak}"
    [[ "$broad" != "true" ]] && final_peak="$filtered_dir/${sample}.peaks.narrowPeak"
    bedtools intersect -v -a "$peak_file" -b "$blacklist" > "$final_peak"
  else
    final_peak="$filtered_dir/$(basename "$peak_file")"
    cp "$peak_file" "$final_peak"
  fi

  if [[ "$broad" != "true" ]]; then
    if [[ -s "$raw_dir/${sample}_summits.bed" ]]; then
      cp "$raw_dir/${sample}_summits.bed" "$summit_dir/${sample}.summits.bed"
    else
      awk 'BEGIN{OFS="\t"} {c=int(($2+$3)/2); if(c<0)c=0; print $1,c,c+1,$4}' "$final_peak" > "$summit_dir/${sample}.summits.bed"
    fi
  fi

  total_reads="$(samtools view -c "$bam")"
  in_peak="NA"
  frip="NA"
  if command -v bedtools >/dev/null; then
    in_peak="$(bedtools intersect -u -abam "$bam" -b "$final_peak" | samtools view -c)"
    frip="$(awk -v total="$total_reads" -v hit="$in_peak" 'BEGIN{print total == 0 ? 0 : hit / total}')"
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$sample" "$bam" "$final_peak" "$total_reads" "$in_peak" "$frip" "$macs_format" >> "$frip_tsv"
done

printf "Finished ATAC peak calling: %s\n" "$outdir/05_peaks"
