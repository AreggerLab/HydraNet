#!/bin/bash
#SBATCH --job-name=all_chr
#SBATCH --array=0-24
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=norm
#SBATCH -t 12:00:00
#SBATCH --output=/path/cleaned_all_guides/chr_%A_%a.out
#SBATCH --error=path/cleaned_all_guides/chr_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mail@mail.com

set -euo pipefail


# ============================================================
# ENVIRONMENT
# ============================================================

export MAMBA_ROOT_PREFIX="/home/nahar2/micromamba"

MICROMAMBA_PATH="/path/micromamba"
CONDA_ENV_PATH="/path//mamba/envs/cas12_r_env_ens14"

eval "$(${MICROMAMBA_PATH} shell hook --shell bash)"
micromamba activate "${CONDA_ENV_PATH}"

echo "======================================================="
echo "TRAVER ALL-GUIDES PER-CHROMOSOME PIPELINE"
echo "Job ID:    ${SLURM_JOB_ID}"
echo "Array ID:  ${SLURM_ARRAY_TASK_ID}"
echo "Started:   $(date)"
echo "======================================================="


# ============================================================
# CHROMOSOME ARRAY
# ============================================================

CHROMS=(
    "chr1"
    "chr2"
    "chr3"
    "chr4"
    "chr5"
    "chr6"
    "chr7"
    "chr8"
    "chr9"
    "chr10"
    "chr11"
    "chr12"
    "chr13"
    "chr14"
    "chr15"
    "chr16"
    "chr17"
    "chr18"
    "chr19"
    "chr20"
    "chr21"
    "chr22"
    "chrX"
    "chrY"
    "chrM"
)

CHR="${CHROMS[$SLURM_ARRAY_TASK_ID]}"

echo "Chromosome: ${CHR}"


# ============================================================
# PATHS
# ============================================================

SCRIPT_DIR="/path"

DATA_DIR="/path/data_dir"

OUTPUT_DIR="${DATA_DIR}/all_guides"

INPUT_ALL="${DATA_DIR}/all_guides_cleaned.csv"

CHR_DIR="${OUTPUT_DIR}/per_chr/${CHR}"

mkdir -p "${CHR_DIR}"


# Scripts

R_SCRIPT_01="${SCRIPT_DIR}/01_get_amino_acid.r"

R_SCRIPT_02="${SCRIPT_DIR}/02_get_domains.r"

R_SCRIPT_03="${SCRIPT_DIR}/03_get_phylop.r"

PY_SCRIPT_04="${SCRIPT_DIR}/04_postprocess.py"



# References

PHYLOP_BW="${SCRIPT_DIR}/hg38.phyloP100way.bw"

SEQINFO_RDS="${SCRIPT_DIR}/hg38_seqinfo.rds"

GTF_FILE="${SCRIPT_DIR}/Homo_sapiens.GRCh38.104.gtf"


# ============================================================
# PER-CHROMOSOME OUTPUTS
# ============================================================

CHR_INPUT="${CHR_DIR}/00_${CHR}_cleaned.csv"

OUT_01="${CHR_DIR}/01_${CHR}_amino_acid.csv"

OUT_02="${CHR_DIR}/02_${CHR}_domains.csv"

OUT_03="${CHR_DIR}/03_${CHR}_phylop.csv"

OUT_04="${CHR_DIR}/04_${CHR}_final_features.csv"



# ============================================================
# SANITY CHECKS
# ============================================================

REQUIRED_FILES=(
    "${INPUT_ALL}"
    "${R_SCRIPT_01}"
    "${R_SCRIPT_02}"
    "${R_SCRIPT_03}"
    "${PY_SCRIPT_04}"
    "${PY_SCRIPT_05}"
    "${PHYLOP_BW}"
    "${SEQINFO_RDS}"
    "${GTF_FILE}"
)

echo
echo "Checking required files..."

for f in "${REQUIRED_FILES[@]}"; do

    if [[ ! -s "${f}" ]]; then
        echo "ERROR: missing or empty:"
        echo "${f}"
        exit 1
    fi

    echo "OK: ${f}"
done


# ============================================================
# SPLIT CLEANED ALL-GUIDES FILE BY CHROMOSOME
# ============================================================

echo
echo "======================================================="
echo "SPLITTING ${CHR}"
echo "======================================================="

rm -f "${CHR_INPUT}"

python - "${INPUT_ALL}" "${CHR_INPUT}" "${CHR}" <<'PY'
import csv
import sys

src = sys.argv[1]
dst = sys.argv[2]
chrom = sys.argv[3]

with open(src, newline="") as fin:
    reader = csv.DictReader(fin)

    if reader.fieldnames is None:
        raise RuntimeError("Input CSV has no header")

    if "seqnames" not in reader.fieldnames:
        raise RuntimeError(
            "Input cleaned CSV does not contain 'seqnames'"
        )

    n = 0

    with open(dst, "w", newline="") as fout:
        writer = csv.DictWriter(
            fout,
            fieldnames=reader.fieldnames
        )

        writer.writeheader()

        for row in reader:

            if str(row.get("seqnames", "")).strip() == chrom:

                writer.writerow(row)
                n += 1

print(f"{chrom}: {n} guides")
PY


# ============================================================
# CSV ROW COUNTER
# ============================================================

csv_rows () {

    python - "$1" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="") as f:

    reader = csv.reader(f)

    try:
        next(reader)
    except StopIteration:
        print(0)
        raise SystemExit

    print(sum(1 for _ in reader))
PY
}


INPUT_ROWS=$(csv_rows "${CHR_INPUT}")

echo
echo "${CHR} input rows: ${INPUT_ROWS}"


# ============================================================
# SKIP CHROMOSOMES WITH NO GUIDES
# ============================================================

if [[ "${INPUT_ROWS}" -eq 0 ]]; then

    echo
    echo "======================================================="
    echo "${CHR}: NO GUIDES — SKIPPING"
    echo "======================================================="

    rm -f "${CHR_INPUT}"

    exit 0
fi


# ============================================================
# OUTPUT VALIDATION
# ============================================================

check_output () {

    FILE="$1"
    STEP="$2"

    if [[ ! -s "${FILE}" ]]; then

        echo
        echo "ERROR: ${STEP} did not create output:"
        echo "${FILE}"
        exit 1
    fi

    ROWS=$(csv_rows "${FILE}")

    echo "${STEP}: ${ROWS} rows"

    if [[ "${ROWS}" -ne "${INPUT_ROWS}" ]]; then

        echo
        echo "ERROR: row count changed during ${STEP}"
        echo "Input:  ${INPUT_ROWS}"
        echo "Output: ${ROWS}"
        echo "File:   ${FILE}"
        exit 1
    fi
}


# ============================================================
# STEP 1
# AMINO ACID / EXON ANNOTATION
# ============================================================

echo
echo "======================================================="
echo "[1/5] ${CHR} — AMINO ACID"
echo "======================================================="

rm -f "${OUT_01}"

Rscript \
    "${R_SCRIPT_01}" \
    "${CHR_INPUT}" \
    "${OUT_01}"

check_output \
    "${OUT_01}" \
    "Step 1"


# ============================================================
# STEP 2
# DOMAIN ANNOTATION
# ============================================================

echo
echo "======================================================="
echo "[2/5] ${CHR} — DOMAINS"
echo "======================================================="

rm -f "${OUT_02}"

Rscript \
    "${R_SCRIPT_02}" \
    "${OUT_01}" \
    "${OUT_02}"

check_output \
    "${OUT_02}" \
    "Step 2"


# ============================================================
# STEP 3
# PHYLOP
# ============================================================

echo
echo "======================================================="
echo "[3/5] ${CHR} — PHYLOP"
echo "======================================================="

rm -f "${OUT_03}"

Rscript \
    "${R_SCRIPT_03}" \
    "${OUT_02}" \
    "${PHYLOP_BW}" \
    "${OUT_03}"

check_output \
    "${OUT_03}" \
    "Step 3"


# ============================================================
# STEP 4
# POSTPROCESS FEATURES
# ============================================================

echo
echo "======================================================="
echo "[4/5] ${CHR} — POSTPROCESS"
echo "======================================================="

rm -f "${OUT_04}"

python -u \
    "${PY_SCRIPT_04}" \
    --input "${OUT_03}" \
    --output "${OUT_04}"

check_output \
    "${OUT_04}" \
    "Step 4"


# ============================================================
# STEP 5
# BACKFILL CUT POSITION + CDS NORMALIZATION
#
# Script 05 overwrites --final, therefore copy Step 04 first.
# ============================================================

echo
echo "======================================================="
echo "[5/5] ${CHR} — CUTPOS + CDS NORMALIZATION"
echo "======================================================="

rm -f "${OUT_05}"

cp \
    "${OUT_04}" \
    "${OUT_05}"

python -u \
    "${PY_SCRIPT_05}" \
    --final "${OUT_05}" \
    --orig "${CHR_INPUT}" \
    --gtf "${GTF_FILE}"

check_output \
    "${OUT_05}" \
    "Step 5"


# ============================================================
# FINAL QC
# ============================================================

echo
echo "======================================================="
echo "${CHR} COMPLETE"
echo "======================================================="
echo "Input guides: ${INPUT_ROWS}"
echo
echo "Final file:"
echo "${OUT_05}"
echo

python - "${OUT_05}" <<'PY'
import pandas as pd
import sys

f = sys.argv[1]

df = pd.read_csv(
    f,
    low_memory=False,
    na_values=["NA"]
)

print(f"Rows: {len(df):,}")

if "cut_position_normalized_CDS" in df.columns:

    exonic = (
        df["exon_idx"].notna()
        if "exon_idx" in df.columns
        else pd.Series(False, index=df.index)
    )

    print(
        "Exonic rows:",
        int(exonic.sum())
    )

    print(
        "Exonic rows with CDS normalized position:",
        int(
            df.loc[
                exonic,
                "cut_position_normalized_CDS"
            ].notna().sum()
        )
    )

if "PhyloP_Nuc_1" in df.columns:

    print(
        "Rows with nucleotide PhyloP:",
        int(df["PhyloP_Nuc_1"].notna().sum())
    )

print()
print("First guide:")

cols = [
    c for c in [
        "chromosome",
        "protospacer_sequence",
        "cut_site",
        "cut_pos",
        "PhyloP_AA_Up1",
        "PhyloP_AA_Cut",
        "PhyloP_AA_Down1",
        "PhyloP_Nuc_1",
        "cut_position_normalized_CDS"
    ]
    if c in df.columns
]

if len(df):
    print(df.loc[df.index[0], cols].to_string())
PY


echo
echo "Finished: $(date)"
echo "======================================================="
