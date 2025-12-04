#!/usr/bin/env bash
#SBATCH --job-name=cas12_predict_all_dual
#SBATCH --mem=256G
#SBATCH --cpus-per-task=2
#SBATCH --partition=unlimited
#SBATCH -t 72:00:00
#SBATCH -o /path/to/logs/slurm_predict_all_dual_%A_%a.out
#SBATCH --array=0-24
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=email@mail.com

set -euo pipefail

# --- Paths ---
PY_BIN="path/to/mamba/envs/cas12_tf_env/bin/python"
MODELS_DIR="/path/to/trained_models"
PREDICT_SCRIPT="path/to/predict_scores.py"
FEATURES_DIR="/path/to/feature_files"

mkdir -p "${MODELS_DIR}/logs"

# Chromosome list aligned with array indices 0..23
CHROMS=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" \
        "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY" "chrM")
CHR="${CHROMS[$SLURM_ARRAY_TASK_ID]}"
FEATURES_FILE="${FEATURES_DIR}/${CHR}_final_features.csv"

ts(){ date +"%Y-%m-%d %H:%M:%S"; }

echo "[$(ts)]   Starting dual prediction for ${CHR}"
echo "[$(ts)]   File: ${FEATURES_FILE}"

# Pre-flight checks
if [ ! -f "$FEATURES_FILE" ]; then
  echo "❌ ERROR: Missing features file: $FEATURES_FILE"
  exit 1
fi
if [ ! -f "$PREDICT_SCRIPT" ]; then
  echo "❌ ERROR: Prediction script missing: $PREDICT_SCRIPT"
  exit 1
fi

# Run dual scoring (exonic + intronic/minimal). Updates file in place.
"$PY_BIN" "$PREDICT_SCRIPT" \
  --features_file "$FEATURES_FILE" \
  --models_dir "$MODELS_DIR"

echo "[$(ts)] 🐉🐉🐉 Finished ${CHR}  (updated: ${FEATURES_FILE})"

