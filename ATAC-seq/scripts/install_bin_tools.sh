#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGET_BIN="${1:-${CONDA_PREFIX:-}/bin}"
[[ -n "$TARGET_BIN" ]] || { echo 'Please provide target bin dir or activate conda env first.' >&2; exit 1; }
mkdir -p "$TARGET_BIN"
chmod +x "$SCRIPT_DIR"/*.sh "$SCRIPT_DIR/atacseq" || true
ln -sf "$SCRIPT_DIR/atacseq" "$TARGET_BIN/atacseq"
echo "[install_bin_tools] linked atacseq -> $TARGET_BIN/atacseq"
for tool in run_tss_qc.sh run_gene_body_profile.sh run_te_heatmap.sh run_te_heatmap_batch.sh run_atac_te_tracks_from_bam.sh run_motif_homer.sh run_tobias.sh run_atac_downstream_only.sh run_cross_contrast_heatmap.sh run_atac_from_bam.sh run_fixedbin_from_bam.sh species_config_lib.sh run_nuc_phasing.sh; do
  ln -sf "$SCRIPT_DIR/$tool" "$TARGET_BIN/${tool%.sh}"
  echo "[install_bin_tools] linked ${tool%.sh} -> $TARGET_BIN/${tool%.sh}"
done
