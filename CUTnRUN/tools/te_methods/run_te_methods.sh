#!/usr/bin/env bash
set -euo pipefail

MANIFEST=""
BAM_DIR=""
OUT_DIR=""
SPECIES="hg38"
METHODS="t3e,allo,repentools"
REPEAT_BED=""
TE_SAF=""
TE_ANNO=""
T3E_DIR=""
ALLO_COMMAND=""
REPENTOOLS_COMMAND=""
EXECUTE="false"
THREADS=8
T3E_PYTHON="${T3E_PYTHON:-}"
T3E_ITERATIONS=100
T3E_ALPHA=0.05
T3E_ENRICHMENT=1.0
# Reuse completed BED conversions.  Converting multi-gigabyte TE BAMs is one
# of the most expensive steps and previously made a partial run impossible to
# resume.  Set T3E_REUSE_BED=false to force regeneration.
T3E_REUSE_BED="${T3E_REUSE_BED:-true}"
# T3E's reference implementation expands every read across its full read
# length and is not tractable on a BAM with tens of millions of multi-mappers.
# A deterministic stride sample keeps the estimator reproducible and bounded;
# set 0 to use every BED row.
T3E_MAX_BED_READS="${T3E_MAX_BED_READS:-1000000}"
ALLO_BIN=""
REPENTOOLS_BIN=""

usage() {
  cat <<'USAGE'
TE method adapters. These methods are optional and run in isolated output folders.

Required:
  --manifest FILE       resolved manifest with target/IgG mappings
  --bam-dir DIR         TE BAM directory, normally results/04_te_bam
  --out-dir DIR         output directory

Options:
  --methods LIST        t3e,allo,repentools (default: all)
  --species STR         hg38 or mm39 (default: hg38)
  --repeat-bed FILE     RepeatMasker BED for T3E/Allo metadata
  --te-saf FILE         TE SAF; used to build RepeatMasker BED when needed
  --te-anno FILE        TE annotation TSV; used to build RepeatMasker BED
  --t3e-dir DIR         cloned T3E repository
  --allo-command CMD    explicit Allo command template
  --repentools-command CMD  explicit RepEnTools command template
  --execute             execute external commands when a complete adapter exists
  --threads INT         threads to expose to command templates (default: 8)
  --t3e-python FILE     Python executable in the isolated T3E environment
  --t3e-iterations INT  T3E simulated input libraries (default: 100)
  --t3e-alpha FLOAT     T3E significance threshold (default: 0.05)
  --t3e-enrichment FLOAT  T3E log2FC threshold (default: 1.0)
  --t3e-reuse-bed BOOL    reuse existing non-empty T3E BED files (default: true)
  --t3e-max-bed-reads INT deterministic BED rows per sample/control (default: 1000000; 0=all)

Without --execute, the adapter only validates inputs and writes reproducible
commands/status. This is intentional because T3E and RepEnTools use their own
reference layouts, while Allo has a separate TensorFlow environment.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --bam-dir) BAM_DIR="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --methods) METHODS="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --repeat-bed) REPEAT_BED="$2"; shift 2 ;;
    --te-saf) TE_SAF="$2"; shift 2 ;;
    --te-anno) TE_ANNO="$2"; shift 2 ;;
    --t3e-dir) T3E_DIR="$2"; shift 2 ;;
    --allo-command) ALLO_COMMAND="$2"; shift 2 ;;
    --repentools-command) REPENTOOLS_COMMAND="$2"; shift 2 ;;
    --execute) EXECUTE="true"; shift ;;
    --threads) THREADS="$2"; shift 2 ;;
    --t3e-python) T3E_PYTHON="$2"; shift 2 ;;
    --t3e-iterations) T3E_ITERATIONS="$2"; shift 2 ;;
    --t3e-alpha) T3E_ALPHA="$2"; shift 2 ;;
    --t3e-enrichment) T3E_ENRICHMENT="$2"; shift 2 ;;
    --t3e-reuse-bed) T3E_REUSE_BED="$2"; shift 2 ;;
    --t3e-max-bed-reads) T3E_MAX_BED_READS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage; exit 1 ;;
  esac
done

[[ -s "$MANIFEST" ]] || { echo "ERROR: --manifest is required" >&2; exit 1; }
[[ -d "$BAM_DIR" ]] || { echo "ERROR: --bam-dir not found: $BAM_DIR" >&2; exit 1; }
[[ -n "$OUT_DIR" ]] || { echo "ERROR: --out-dir is required" >&2; exit 1; }
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$ROOT_DIR/../.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || command -v python || true)}"
[[ -n "$PYTHON_BIN" ]] || { echo "ERROR: python3/python not found" >&2; exit 1; }
SAMTOOLS_BIN="${SAMTOOLS_BIN:-$(command -v samtools || true)}"
BEDTOOLS_BIN="${BEDTOOLS_BIN:-$(command -v bedtools || true)}"
[[ -x "$SAMTOOLS_BIN" ]] || [[ -x /path/to/.conda/envs/chipseq/bin/samtools ]] && SAMTOOLS_BIN="${SAMTOOLS_BIN:-/path/to/.conda/envs/chipseq/bin/samtools}"
[[ -x "$BEDTOOLS_BIN" ]] || [[ -x /path/to/.conda/envs/chipseq/bin/bedtools ]] && BEDTOOLS_BIN="${BEDTOOLS_BIN:-/path/to/.conda/envs/chipseq/bin/bedtools}"
if [[ -z "$T3E_PYTHON" ]]; then
  for candidate in "/path/to/.conda/envs/cutrun-t3e/bin/python" "${HOME}/.conda/envs/cutrun-t3e/bin/python"; do
    if [[ -x "$candidate" ]]; then T3E_PYTHON="$candidate"; break; fi
  done
fi
for candidate in "/path/to/.conda/envs/cutrun-allo/bin/allo" "${HOME}/.conda/envs/cutrun-allo/bin/allo"; do
  if [[ -x "$candidate" ]]; then ALLO_BIN="$candidate"; break; fi
done
for candidate in "/path/to/.local/share/CUTnRUN-tools/RepEnTools/ret" "${HOME}/.local/share/CUTnRUN-tools/RepEnTools/ret"; do
  if [[ -x "$candidate" ]]; then REPENTOOLS_BIN="$candidate"; break; fi
done
# Allo accepts BAM input and writes SAM output. CUT&RUN TE BAMs can lack the
# optional PG/collate tag; --ignore makes Allo process those valid alignments
# instead of returning success without writing an output file. Keep this
# default explicit and overridable because paired-end is the normal layout.
if [[ -z "$ALLO_COMMAND" && -n "$ALLO_BIN" ]]; then
  ALLO_COMMAND="$ALLO_BIN {input_bam} -seq pe -o {output_bam} -p {threads} --ignore"
fi
mkdir -p "$OUT_DIR"
if [[ -x "$ROOT_DIR/verify_installation.sh" ]]; then
  "$ROOT_DIR/verify_installation.sh" "$OUT_DIR/tool_installation.tsv" >/dev/null || true
fi
"$PYTHON_BIN" "$ROOT_DIR/te_method_plan.py" --manifest "$MANIFEST" --output-dir "$OUT_DIR"
if [[ -z "$TE_SAF" ]]; then TE_SAF="$REPO_ROOT/resources/CUTRUN_analysis/anno/te_${SPECIES}.saf"; fi
if [[ -z "$TE_ANNO" ]]; then TE_ANNO="$REPO_ROOT/resources/CUTRUN_analysis/anno/te_anno_${SPECIES}.tsv"; fi
if [[ -z "$REPEAT_BED" && -s "$TE_SAF" && -s "$TE_ANNO" ]]; then
  REPEAT_BED="$OUT_DIR/repeat_annotations.bed"
  "$PYTHON_BIN" "$ROOT_DIR/build_repeat_bed.py" --te-saf "$TE_SAF" --te-anno "$TE_ANNO" --output "$REPEAT_BED"
fi
STATUS_FILE="$OUT_DIR/method_status.tsv"
if [[ -s "$STATUS_FILE" ]]; then
  # The adapter status is an attempt-level file.  Preserve a previous run as
  # legacy instead of appending to it; otherwise an old PASS would make a
  # current all-SKIP/FAIL attempt exit successfully.
  mv "$STATUS_FILE" "${STATUS_FILE%.tsv}.legacy_$(date -u +%Y%m%dT%H%M%SZ)_$$.tsv"
fi
printf 'method\tstatus\tdetail\n' > "$STATUS_FILE"

status() {
  printf '%s\t%s\t%s\n' "$1" "$2" "$3" >> "$STATUS_FILE"
  echo "[te_methods] $1: $2 — $3"
}

resolve_bam() {
  local sample="$1"
  for candidate in "$BAM_DIR/${sample}_te.bam" "$BAM_DIR/${sample}.bam" "$BAM_DIR/${sample}_clean.bam"; do
    if [[ -s "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

infer_read_length() {
  # awk exits after the first alignment; with `set -o pipefail` samtools then
  # receives SIGPIPE (141), which used to abort the whole adapter before any
  # method status was written.  Keep the early exit but isolate that expected
  # SIGPIPE in a subshell with pipefail disabled.
  (
    set +o pipefail
    "$SAMTOOLS_BIN" view "$1" | awk 'length($10) > 0 {print length($10); exit}'
  )
}

render_template() {
  local template="$1" sample="$2" control="$3" input_bam="$4" output_bam="$5" method_dir="$6"
  template="${template//\{sample\}/$sample}"
  template="${template//\{control\}/$control}"
  template="${template//\{input_bam\}/$input_bam}"
  template="${template//\{output_bam\}/$output_bam}"
  template="${template//\{repeat_bed\}/$REPEAT_BED}"
  template="${template//\{manifest\}/$MANIFEST}"
  template="${template//\{out_dir\}/$method_dir}"
  template="${template//\{species\}/$SPECIES}"
  template="${template//\{threads\}/$THREADS}"
  printf '%s\n' "$template"
}

IFS=',' read -r -a selected <<< "$METHODS"
for method in "${selected[@]}"; do
  method="$(echo "$method" | tr '[:upper:]' '[:lower:]' | xargs)"
  case "$method" in
    t3e)
      method_dir="$OUT_DIR/t3e"
      mkdir -p "$method_dir"
      {
        echo "# T3E command plan"
        echo "# species=${SPECIES}"
        echo "# bam_dir=${BAM_DIR}"
        echo "# repeat_bed=${REPEAT_BED:-MISSING}"
        echo "# T3E expects BAM files with secondary mappings and a control_sample.txt/parameters.txt layout."
        echo "# Install: conda env create -f ${ROOT_DIR}/environment-t3e.yml"
        echo "# Execute with --execute --t3e-python /path/to/cutrun-t3e/bin/python"
      } > "$method_dir/COMMAND_PLAN.txt"
      t3e_species="$SPECIES"
      [[ "$t3e_species" == mm39 ]] && t3e_species="mm10"
      if [[ ! -d "$T3E_DIR" || ! -f "$T3E_DIR/scripts/probabilities.py" || ! -f "$T3E_DIR/scripts/t3e.py" || ! -f "$T3E_DIR/scripts/enrichment.py" ]]; then
        status t3e SKIP "T3E repository not supplied; see t3e/COMMAND_PLAN.txt"
      elif [[ -z "$REPEAT_BED" || ! -s "$REPEAT_BED" ]]; then
        status t3e SKIP "--repeat-bed is required for T3E"
      elif [[ "$t3e_species" != hg38 && "$t3e_species" != mm10 ]]; then
        status t3e SKIP "T3E upstream reference profile supports hg38/mm10; got ${SPECIES}"
      elif [[ ! -x "$SAMTOOLS_BIN" || ! -x "$BEDTOOLS_BIN" ]]; then
        status t3e SKIP "samtools and bedtools are required in the execution environment"
      elif [[ "$EXECUTE" != true ]]; then
        status t3e READY "repository and RepeatMasker BED detected; rerun with --execute"
      else
        t3e_python="${T3E_PYTHON:-$PYTHON_BIN}"
        [[ -x "$t3e_python" ]] || { status t3e FAIL "T3E Python not executable: $t3e_python"; continue; }
        t3e_bed_dir="$method_dir/bed"
        t3e_results="$method_dir/results"
        t3e_prob="$method_dir/probabilities"
        t3e_repeat_bed="$method_dir/repeat_bed4.bed"
        awk 'BEGIN {OFS="\t"} !/^#/ && NF >= 4 {print $1, $2, $3, $4}' "$REPEAT_BED" > "$t3e_repeat_bed"
        [[ -s "$t3e_repeat_bed" ]] || { status t3e FAIL "RepeatMasker BED has no BED4 rows"; continue; }
        t3e_logs="$method_dir/logs"
        mkdir -p "$t3e_bed_dir" "$t3e_results" "$t3e_prob" "$t3e_logs"
        t3e_ok=true
        t3e_pass=false
        t3e_skip=false
        t3e_targets=0
        while IFS=$'\t' read -r sample is_igg control; do
          sample="${sample%$'\r'}"
          is_igg="${is_igg%$'\r'}"
          control="${control%$'\r'}"
          [[ "$is_igg" == "True" ]] && continue
          t3e_targets=$((t3e_targets + 1))
          [[ -n "$control" ]] || { status "t3e:${sample}" SKIP "no control mapping"; t3e_ok=false; continue; }
          sample_bam="$(resolve_bam "$sample")" || { status "t3e:${sample}" FAIL "sample BAM missing"; t3e_ok=false; continue; }
          control_bam="$(resolve_bam "$control")" || { status "t3e:${sample}" FAIL "control BAM missing: $control"; t3e_ok=false; continue; }
          sample_bed_raw="$t3e_bed_dir/${sample}.bed"
          control_bed_raw="$t3e_bed_dir/${control}.bed"
          # T3E's parsers require BED4 and use column 4 as the read id.
          if [[ "$T3E_REUSE_BED" != true || ! -s "$sample_bed_raw" ]]; then
            if ! "$BEDTOOLS_BIN" bamtobed -i "$sample_bam" | cut -f1-4 > "$sample_bed_raw" 2>"$t3e_logs/${sample}.bamtobed.log"; then
              status "t3e:${sample}" FAIL "sample BAM to BED failed; see ${sample}.bamtobed.log"
              t3e_ok=false
              continue
            fi
          fi
          if [[ "$T3E_REUSE_BED" != true || ! -s "$control_bed_raw" ]]; then
            if ! "$BEDTOOLS_BIN" bamtobed -i "$control_bam" | cut -f1-4 > "$control_bed_raw" 2>"$t3e_logs/${control}.bamtobed.log"; then
              status "t3e:${sample}" FAIL "control BAM to BED failed; see ${control}.bamtobed.log"
              t3e_ok=false
              continue
            fi
          fi
          sample_bed="$sample_bed_raw"
          control_bed="$control_bed_raw"
          if [[ "$T3E_MAX_BED_READS" =~ ^[0-9]+$ && "$T3E_MAX_BED_READS" -gt 0 ]]; then
            for bed_pair in sample control; do
              if [[ "$bed_pair" == sample ]]; then
                raw_bed="$sample_bed_raw"
                sampled_bed="$t3e_bed_dir/${sample}.sampled.bed"
              else
                raw_bed="$control_bed_raw"
                sampled_bed="$t3e_bed_dir/${control}.sampled.bed"
              fi
              # Reuse a previously generated sampled BED when it is newer than
              # the raw BED.  Counting lines in multi-GB BEDs over NFS is very
              # expensive and adds no value for a resumable run.
              if [[ -s "$sampled_bed" && ! "$raw_bed" -nt "$sampled_bed" ]]; then
                if [[ "$bed_pair" == sample ]]; then sample_bed="$sampled_bed"; else control_bed="$sampled_bed"; fi
                continue
              fi
              raw_rows="$(wc -l < "$raw_bed")"
              if (( raw_rows > T3E_MAX_BED_READS )); then
                stride=$(( (raw_rows + T3E_MAX_BED_READS - 1) / T3E_MAX_BED_READS ))
                if [[ ! -s "$sampled_bed" || "$raw_bed" -nt "$sampled_bed" ]]; then
                  awk -v stride="$stride" 'NR == 1 || (NR - 1) % stride == 0' "$raw_bed" > "$sampled_bed"
                fi
                if [[ "$bed_pair" == sample ]]; then sample_bed="$sampled_bed"; else control_bed="$sampled_bed"; fi
              fi
            done
          fi
          readlen="$(infer_read_length "$control_bam")"
          [[ "$readlen" =~ ^[0-9]+$ ]] || { status "t3e:${sample}" FAIL "could not infer read length"; t3e_ok=false; continue; }
          control_counts="$t3e_results/${control}_counts.txt"
          sample_counts="$t3e_results/${sample}_counts.txt"
          if [[ ! -s "$control_counts" ]]; then
            if ! "$BEDTOOLS_BIN" intersect -a "$t3e_repeat_bed" -b "$control_bed" -c 2>"$t3e_logs/${control}.intersect.log" | awk -F '\t' 'BEGIN {OFS="\t"} {x[$4]+=$NF} END {for (k in x) print k, x[k]}' | sort > "$control_counts"; then
              status "t3e:${sample}" FAIL "control repeat counts failed; see ${control}.intersect.log"
              t3e_ok=false
              continue
            fi
          fi
          if [[ ! -s "$sample_counts" ]]; then
            if ! "$BEDTOOLS_BIN" intersect -a "$t3e_repeat_bed" -b "$sample_bed" -c 2>"$t3e_logs/${sample}.intersect.log" | awk -F '\t' 'BEGIN {OFS="\t"} {x[$4]+=$NF} END {for (k in x) print k, x[k]}' | sort > "$sample_counts"; then
              status "t3e:${sample}" FAIL "sample repeat counts failed; see ${sample}.intersect.log"
              t3e_ok=false
              continue
            fi
          fi
          probability_dir="$t3e_prob/${control}"
          mkdir -p "$probability_dir"
          if ! "$t3e_python" "$T3E_DIR/scripts/probabilities.py" --control "$control_bed" --readlen "$readlen" --species "$t3e_species" --outputfolder "$probability_dir" >"$t3e_logs/${control}.probabilities.log" 2>&1; then
            status "t3e:${sample}" FAIL "probabilities.py failed; see adapter log"
            t3e_ok=false
            continue
          fi
          target_dir="$t3e_results/$sample"
          mkdir -p "$target_dir"
          if ! "$t3e_python" "$T3E_DIR/scripts/t3e.py" --repeat "$t3e_repeat_bed" --sample "$sample_bed" --readlen "$readlen" --control "$control_bed" --controlcounts "$control_counts" --probability "$probability_dir" --iter "$T3E_ITERATIONS" --species "$t3e_species" --outputfolder "$target_dir" --outputprefix "$sample" >"$t3e_logs/${sample}.t3e.log" 2>&1; then
            if grep -qiE 'error.*list of reads|unsupported|not compatible' "$t3e_logs/${sample}.t3e.log"; then
              status "t3e:${sample}" SKIP "upstream T3E input parser is incompatible with this dataset; no scientific result"
              t3e_skip=true
            else
              status "t3e:${sample}" FAIL "t3e.py failed; see adapter log"
            fi
            t3e_ok=false
            continue
          fi
          background="$target_dir/${sample}_background.txt"
          if ! "$t3e_python" "$T3E_DIR/scripts/enrichment.py" --background "$background" --signal "$sample_counts" --iter "$T3E_ITERATIONS" --alpha "$T3E_ALPHA" --enrichment "$T3E_ENRICHMENT" --outputfolder "$target_dir" --outputprefix "$sample" >"$t3e_logs/${sample}.enrichment.log" 2>&1; then
            status "t3e:${sample}" FAIL "enrichment.py failed; see adapter log"
            t3e_ok=false
            continue
          fi
          t3e_result_count="$(find "$target_dir" -maxdepth 1 -type f -size +0c ! -name '*.log' | wc -l | tr -d ' ')"
          if [[ "$t3e_result_count" -eq 0 ]]; then
            status "t3e:${sample}" FAIL "enrichment.py returned 0 but wrote no non-empty result file"
            t3e_ok=false
            continue
          fi
          status "t3e:${sample}" PASS "family/subfamily enrichment completed"
          t3e_pass=true
        done < <(tail -n +2 "$OUT_DIR/method_plan.tsv")
        if [[ "$t3e_targets" -eq 0 ]]; then
          status t3e SKIP "no target samples were available in method_plan.tsv"
        elif [[ "$t3e_ok" == true ]]; then
          status t3e PASS "all mapped targets completed"
        elif [[ "$t3e_pass" == false && "$t3e_skip" == true ]]; then
          status t3e SKIP "T3E upstream parser incompatible; no scientifically valid result"
        else
          status t3e FAIL "one or more targets failed"
        fi
      fi
      ;;
    allo)
      method_dir="$OUT_DIR/allo"
      mkdir -p "$method_dir"
      printf '%s\n' "# Allo command template" "# input_bam={input_bam}" "# output_bam={output_bam} (SAM output)" "# threads=${THREADS}" "# resolved=${ALLO_COMMAND:-MISSING}" > "$method_dir/COMMAND_PLAN.txt"
      if [[ -z "$ALLO_COMMAND" ]]; then
        status allo SKIP "Allo is not installed and --allo-command was not supplied"
      elif [[ "$EXECUTE" != true ]]; then
        status allo READY "Allo command is available; rerun with --execute after reviewing the template"
      else
        if [[ -z "$ALLO_COMMAND" ]]; then
          status allo SKIP "Allo flags are version-specific; supply an explicit --allo-command template"
        else
          allo_ok=true
          allo_targets=0
          allo_skip=false
          while IFS=$'\t' read -r sample is_igg control; do
            sample="${sample%$'\r'}"
            is_igg="${is_igg%$'\r'}"
            control="${control%$'\r'}"
            [[ "$is_igg" == "True" ]] && continue
            allo_targets=$((allo_targets + 1))
            input_bam="$(resolve_bam "$sample")" || { status "allo:${sample}" FAIL "sample BAM missing"; allo_ok=false; continue; }
            output_bam="$method_dir/${sample}_allo.sam"
            sorted_bam="$method_dir/${sample}_allo.bam"
            command_text="$(render_template "$ALLO_COMMAND" "$sample" "$control" "$input_bam" "$output_bam" "$method_dir")"
            printf '%s\n' "$command_text" >> "$method_dir/commands.sh"
            if [[ -s "$sorted_bam" && -s "$sorted_bam.bai" && "$output_bam" -ot "$sorted_bam" && "$output_bam" -ot "$sorted_bam.bai" ]]; then
              status "allo:${sample}" PASS "validated BAM/BAI already present (resume)"
              continue
            fi
            bash -lc "$command_text" >"$method_dir/${sample}.stdout.log" 2>"$method_dir/${sample}.stderr.log" || { status "allo:${sample}" FAIL "command failed; see stderr log"; allo_ok=false; continue; }
            [[ -s "$output_bam" ]] || { status "allo:${sample}" FAIL "command completed but output SAM is missing"; allo_ok=false; continue; }
            # Allo can return 0 after writing a diagnostic/empty file.  Before
            # accepting it, require at least one SAM alignment (or a header
            # plus an alignment) so a text log is never treated as a BAM.
            if ! awk 'BEGIN {ok=0} /^@/ {next} NF >= 11 {ok=1; exit} END {exit(ok ? 0 : 1)}' "$output_bam"; then
              status "allo:${sample}" FAIL "command completed but output is not a valid SAM alignment file"; allo_ok=false; continue
            fi
            if [[ -x "$SAMTOOLS_BIN" ]]; then
              if ! "$SAMTOOLS_BIN" view -bS "$output_bam" | "$SAMTOOLS_BIN" sort -@ "$THREADS" -o "$sorted_bam" -; then
                status "allo:${sample}" FAIL "SAM-to-BAM conversion failed; see stderr log"; allo_ok=false; continue
              fi
              if ! "$SAMTOOLS_BIN" index -@ "$THREADS" "$sorted_bam" || ! "$SAMTOOLS_BIN" flagstat "$sorted_bam" > "$method_dir/${sample}_allo.flagstat.txt"; then
                status "allo:${sample}" FAIL "BAM index/flagstat failed"; allo_ok=false; continue
              fi
              status "allo:${sample}" PASS "SAM, sorted BAM, BAI and flagstat validated"
            else
              status "allo:${sample}" SKIP "SAM validated; samtools unavailable for BAM conversion"
              allo_skip=true
            fi
          done < <(tail -n +2 "$OUT_DIR/method_plan.tsv")
          if [[ "$allo_targets" -eq 0 ]]; then
            status allo SKIP "no target samples were available in method_plan.tsv"
          elif [[ "$allo_ok" == true && "$allo_skip" == true ]]; then
            status allo SKIP "SAM outputs validated but sorted BAM/BAI require samtools"
          elif [[ "$allo_ok" == true ]]; then
            status allo PASS "all mapped targets completed"
          else
            status allo FAIL "one or more targets failed"
          fi
        fi
      fi
      ;;
    repentools|repen-tools)
      method_dir="$OUT_DIR/repentools"
      mkdir -p "$method_dir"
      {
        echo "# RepEnTools command plan"
        echo "# RepEnTools is a chm13v2 FASTQ workflow; current hg38/mm39 BAMs are not silently converted."
        echo "# command=${REPENTOOLS_COMMAND:-MISSING}"
        echo "# Use a dedicated chm13v2 reference profile and record its RMSK version."
      } > "$method_dir/COMMAND_PLAN.txt"
      if [[ -z "$REPENTOOLS_COMMAND" && -z "$REPENTOOLS_BIN" ]]; then
        status repentools SKIP "RepEnTools command is not installed; chm13v2 profile required"
      elif [[ "$SPECIES" != chm13v2 ]]; then
        status repentools SKIP "RepEnTools adapter requires an explicit chm13v2 reference"
      elif [[ "$EXECUTE" != true ]]; then
        status repentools READY "command detected; rerun with --execute after reviewing the chm13v2 plan"
      else
        if [[ -z "$REPENTOOLS_COMMAND" ]]; then
          status repentools SKIP "RepEnTools flags are workflow-specific; supply --repentools-command"
        else
          command_text="$(render_template "$REPENTOOLS_COMMAND" "all" "" "" "$method_dir" "$method_dir")"
          printf '%s\n' "$command_text" > "$method_dir/commands.sh"
          bash -lc "$command_text" >"$method_dir/repentools.stdout.log" 2>"$method_dir/repentools.stderr.log" || { status repentools FAIL "command failed; see stderr log"; continue; }
          status repentools PASS "chm13v2 command completed"
        fi
      fi
      ;;
    *) status "$method" ERROR "unknown method" ;;
  esac
done

# Make the standalone adapter produce the same figure contract as the main
# downstream entry point.  The renderer is non-fatal and records SKIP/FAIL in
# visualization_status.tsv when an upstream method has no scientific output.
if [[ -f "$ROOT_DIR/render_te_method_visuals.py" ]]; then
  "$PYTHON_BIN" "$ROOT_DIR/render_te_method_visuals.py" \
    --methods-dir "$OUT_DIR" --output-dir "$OUT_DIR/visuals" >/dev/null || true
fi

echo "[te_methods] status: $STATUS_FILE"
# Aggregate rows (t3e/allo) determine the attempt result.  Per-sample PASS
# rows must not hide an aggregate FAIL when only one replicate completed.
if awk -F '\t' 'NR > 1 && $1 !~ /:/ && $2 == "FAIL" {found=1} END {exit found ? 0 : 1}' "$STATUS_FILE"; then
  exit 1
fi
if awk -F '\t' 'NR > 1 && $1 !~ /:/ && $2 == "ERROR" {found=1} END {exit found ? 0 : 1}' "$STATUS_FILE"; then
  exit 1
fi
if awk -F '\t' 'NR > 1 && $1 !~ /:/ && $2 == "PASS" {found=1} END {exit found ? 0 : 1}' "$STATUS_FILE"; then
  exit 0
fi
# 42 is the shared optional-module SKIP code consumed by run_downstream.sh.
exit 42
