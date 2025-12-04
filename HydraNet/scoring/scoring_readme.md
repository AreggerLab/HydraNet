# HydraNet Scoring Pipeline

This folder contains code and models for scoring guides with HydraNet using precomputed features.

There are two model families:

- **CDS models (`*_cds`)**: trained on exonic guides with full protein, conservation, and sequence features.
- **Seq models (`*_seq`)**: trained on sequence-only features and can score any guide with a valid protospacer.

Models are expected under:

```
scoring/trained_models/
    D156R_cds/
    KX10_cds/
    KX12_cds/
    D156R_seq/
    KX10_seq/
    KX12_seq/
```

Each model directory contains a TensorFlow SavedModel export (e.g. `saved_model.pb`, `variables/`) and `temperature.json` with a single key `temperature` used for calibration.

---

## 1. 🐉🐉🐉  Dependencies (cas12_tf_env)

Exact versions from inference environment:

```
python         3.11.13
tensorflow     2.18.0
scikit-learn   1.7.1
pandas         2.3.2
numpy          2.3.2
```

Standard library modules used:

```
os, json, warnings, argparse, collections
```

Install with micromamba (example):

```bash
micromamba install -n cas12_tf_env \
  tensorflow=2.18.0 scikit-learn=1.7.1 pandas=2.3.2 numpy=2.3.2
```

---

## 2. 🐉🐉🐉  Input Requirements

The scoring script expects a **feature CSV** produced by the feature_collection pipeline (e.g. `final_features.csv`).

### 2.1 Required columns for protospacer handling

The script will use `protospacer_sequence` to construct one-hot encoding and embeddings


Rules:

- For minus strand (`strand == "-"`): `protospacer_sequence = pam_sequence` (must be length 23).
- For plus strand (`strand == "+"`): protospacer is taken from `guide_sequence` after removing the PAM prefix, then truncated to 23 nt.

At least one of the following must therefore be true:

- `protospacer_sequence` exists and has length 23 for all rows, **or**
- `guide_sequence`, `pam_sequence`, and `strand` exist with lengths that allow reconstruction.

### 2.2 Required feature columns for CDS (exonic) models

CDS models are only applied to rows with a non-null `exon_idx`. For those rows, the following columns must be present (they are produced by the feature_collection pipeline):

- PhyloP AA and flank:
  - `PhyloP_AA_Up1`, `PhyloP_AA_Cut`, `PhyloP_AA_Down1`, `PhyloP_GuidePlusFlank`
- Domain and CDS position:
  - `Binary_Domain`, `cut_position_normalized_CDS`
- Sequence-level features:
  - `GC_Content`, `Melting_Temperature`, `RNAfold_MFE`, `has_TTT`,
  - `AT_skew`, `palindromic_score`, `self_complementarity_score`
- Amino-acid biochemical features (cut site, upstream, downstream, averaged):
  - `Cut_Site_*`, `Upstream1_AA_*`, `Downstream1_AA_*`, `Avg_3_*`
- Nucleotide-level PhyloP (23 positions):
  - `PhyloP_Nuc_1` ... `PhyloP_Nuc_23`

The script uses a fixed list of feature columns for the exonic branch. Any rows missing, or having NA values, these guides are dropped from CDS model evaluation and processed using Seq model.

### 2.3 Required feature columns for Seq models

Seq models use a smaller tabular feature set plus the one-hot encoded protospacer. They require:

- `GC_Content`
- `Melting_Temperature`
- `RNAfold_MFE`
- `has_TTT`
- `AT_skew`
- `palindromic_score`
- `self_complementarity_score`

These are expected to be produced by `04_postprocess.py`.

All rows with a valid 23-nt `protospacer_sequence` are eligible for Seq model scoring.

---

## 3. 🐉🐉🐉  Script: predict_scores.py

The main scoring script performs:

- Exonic scoring with 4-input HydraNet-style models (`*_cds`).
- Minimal scoring with 2-input sequence models (`*_seq`).
- Writes six new columns into the same CSV, in place.

### Usage

```bash
python3 predict_scores.py \
  --features_file path/to/final_features.csv \
  --models_dir    path/to/trained_models
```

After running, the input CSV is modified in place and gains the following columns:

- `D156R_CDS`
- `KX10_CDS`
- `KX12_CDS`
- `D156R_Seq`
- `KX10_Seq`
- `KX12_Seq`

Each column contains a probability score in [0, 1].

---

## 4. 🐉🐉🐉 SLURM Array Example

The provided SLURM script (`predict_all_dual.sbatch`) runs scoring per chromosome using a job array, updating one feature file per task. The paths and file names can be adapted

Example template:

```bash
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

PY_BIN="path/to/mamba/envs/cas12_tf_env/bin/python"
MODELS_DIR="/path/to/trained_models"
PREDICT_SCRIPT="path/to/predict_scores.py"
FEATURES_DIR="/path/to/feature_files"

CHROMS=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" \
        "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY" "chrM")
CHR="${CHROMS[$SLURM_ARRAY_TASK_ID]}"
FEATURES_FILE="${FEATURES_DIR}/${CHR}_final_features.csv"

"$PY_BIN" "$PREDICT_SCRIPT" \
  --features_file "$FEATURES_FILE" \
  --models_dir "$MODELS_DIR"
```

Edit the paths (`PY_BIN`, `MODELS_DIR`, `PREDICT_SCRIPT`, `FEATURES_DIR`) to match your cluster layout.

---

## 5. 🐉🐉🐉  Outputs

The scoring pipeline does **not** create new files by default. It overwrites the feature CSV specified by `--features_file` by adding six new **score columns** in place.

**Keep a copy of the original features file if you need an unscored version.**
