#!/usr/bin/env bash
# Run the whole pipeline:   bash run_pipeline.sh [config/config.yaml]
# Resume from a step:       START_AT=4 bash run_pipeline.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

CONFIG="${1:-config/config.yaml}"
START_AT="${START_AT:-1}"
cfg() { python3 scripts/cfg_get.py "$CONFIG" "$1"; }

RESULTS="$(cfg paths.results_dir)"
mkdir -p "$RESULTS/logs"
LOG="$RESULTS/logs/pipeline.log"

step=0
run() {            # run <description> <command...>
  step=$((step + 1))
  local desc="$1"; shift
  if [ "$step" -lt "$START_AT" ]; then echo "[skip] $step $desc"; return; fi
  echo "=== STEP $step: $desc" | tee -a "$LOG"
  "$@" 2>&1 | tee -a "$LOG"
}

run "Decompress / validate IDATs"          bash scripts/01_prepare_idats.sh "$(cfg paths.idat_dir)" "$(cfg paths.work_idat_dir)"
run "QC: detection p-values (SeSAMe)"      Rscript scripts/02_qc.R "$CONFIG"
run "Preprocess (SeSAMe) + beta density"   Rscript scripts/03_preprocess.R "$CONFIG"
run "PCA + top-variable heatmap"           Rscript scripts/04_exploratory.R "$CONFIG"
run "UMAP"                                 python3 scripts/05_umap.py "$CONFIG"
run "Differential methylation (SeSAMe DML)" Rscript scripts/06_dml.R "$CONFIG"
run "Annotation + volcano plots"           Rscript scripts/07_annotate_volcano.R "$CONFIG"
run "Pathway ORA + GSEA"                   Rscript scripts/08_pathways.R "$CONFIG"
run "Lollipop plots"                       python3 scripts/09_lollipop.py "$CONFIG"

echo "Pipeline finished. Outputs are in: $RESULTS/{figures,tables}" | tee -a "$LOG"
