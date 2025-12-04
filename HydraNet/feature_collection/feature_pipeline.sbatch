#!/bin/bash

#SBATCH --job-name=cas12_feature
#SBATCH --mem=256G
#SBATCH --cpus-per-task=12
#SBATCH --partition=norm
#SBATCH -t 240:00:00
#SBATCH --output=path/to/logs/slurm_FULL_%A_%a.out
#SBATCH --error=path/to/logs/slurm_FULL_%A_%a.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=email@mail.com
#SBATCH --array=0-21

set -eo pipefail

# Chromosome list: mm10 (1–19, X, Y, M) or (1-22) for hg38
CHR_LIST=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY chrM)
CHR_NAME=${CHR_LIST[${SLURM_ARRAY_TASK_ID}]}

echo "======================================================="
echo "🐉🐉🐉Initializing Cas12 FULL PIPELINE Job (MOUSE)"
echo "🐉🐉🐉Job ID: ${SLURM_JOB_ID:-unset}  Array Task: ${SLURM_ARRAY_TASK_ID}  Chr: ${CHR_NAME}"
echo "======================================================="

echo "🐉🐉🐉Activating micromamba environment..."
eval "$(micromamba shell hook --shell bash)"
micromamba activate cas12_r_env
echo "🐉🐉🐉Environment 'cas12_r_env' activated."

# --- PATHS ---

R_SCRIPT_01="path/to/01_get_amino_acid.r"
R_SCRIPT_02="path/to/02_get_domains.r"
R_SCRIPT_03="path/to/03_get_phylop.r"
PY_SCRIPT_04="path/to/04_postprocess.py"

#change as needed
DATE_STAMP="mm.dd.yy"
DATE_DIR="path/to/${DATE_STAMP}"
TMP_DIR="${DATE_DIR}/temp_files"
FINAL_OUTPUT_DIR="${DATE_DIR}/feature_files"

mkdir -p "${TMP_DIR}" "${FINAL_OUTPUT_DIR}"
PHYLOP_BW="path/to/mm10.60way.phyloP60way.bw"
#change to human as needed

INPUT_FILE="path/to/${CHR_NAME}_input.csv"

OUT_01="${TMP_DIR}/${CHR_NAME}_step1_amino.csv"
OUT_02="${TMP_DIR}/${CHR_NAME}_step2_domains.csv"
OUT_03="${TMP_DIR}/${CHR_NAME}_step3_phylop.csv"
FINAL_CHUNK_OUT="${FINAL_OUTPUT_DIR}/${CHR_NAME}_final_features.csv"

echo "Sanity checks for ${CHR_NAME}..."
if [[ ! -f "${R_SCRIPT_01}" ]]; then echo "❌ Script missing: ${R_SCRIPT_01}"; exit 1; fi
if [[ ! -f "${R_SCRIPT_02}" ]]; then echo "❌ Script missing: ${R_SCRIPT_02}"; exit 1; fi
if [[ ! -f "${R_SCRIPT_03}" ]]; then echo "❌ Script missing: ${R_SCRIPT_03}"; exit 1; fi
if [[ ! -f "${PY_SCRIPT_04}" ]]; then echo "❌ Script missing: ${PY_SCRIPT_04}"; exit 1; fi
if [[ ! -f "${PHYLOP_BW}" ]];   then echo "❌ BigWig file missing: ${PHYLOP_BW}"; exit 1; fi
if [[ ! -s "${INPUT_FILE}" ]]; then echo "❌ Input file missing or empty: ${INPUT_FILE}"; exit 1; fi
echo "🐲🐲🐲 All required files found for ${CHR_NAME}."

echo "🐉🐉🐉 Starting FULL MOUSE pipeline for: ${CHR_NAME}"
echo "🐉🐉🐉Input file: ${INPUT_FILE}"

echo "🐉🐉🐉--- [1/4] Amino Acid Annotation (${CHR_NAME}) ---"
Rscript "${R_SCRIPT_01}" "${INPUT_FILE}" "${OUT_01}"

echo "🐉🐉🐉--- [2/4] Domain Annotation (${CHR_NAME}) ---"
Rscript "${R_SCRIPT_02}" "${OUT_01}" "${OUT_02}"

echo "🐉🐉🐉--- [3/4] PhyloP Annotation (${CHR_NAME}) ---"
Rscript "${R_SCRIPT_03}" "${OUT_02}" "${PHYLOP_BW}" "${OUT_03}"

echo "🐉🐉🐉--- [4/4] Final Feature Postprocessing (${CHR_NAME}) ---"
python3 "${PY_SCRIPT_04}" --input "${OUT_03}" --output "${FINAL_CHUNK_OUT}"

echo "🐉🐉🐉 Mouse pipeline finished successfully for ${CHR_NAME}."
echo "======================================================="

