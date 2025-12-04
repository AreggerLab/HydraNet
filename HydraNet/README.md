# HydraNet: Multi-Head Deep Learning for Cas12a Guide Scoring

HydraNet is a hybrid deep-learning system for predicting Cas12a guide activity using both protein-context and sequence-only features. This repository provides:

- A feature-collection pipeline for transforming raw genomic coordinates + protospacer sequences into ML-ready feature tables.
- A scoring pipeline using pretrained HydraNet models to assign activity probabilities.
- Example inputs for both stages.

## Model Families

### 1. CDS Models (`*_cds`)
For **exonic guides**, using:
- amino-acid context
- protein domain annotations
- PhyloP conservation (codon + nucleotide)
- biochemical AA features
- sequence-derived metrics

### 2. Seq Models (`*_seq`)
For **all guides**, including non-coding:
- protospacer one-hot encoding
- GC content
- MFE (RNAfold)
- k-mer/palindromic structure features

---

## Repository Structure

```
HydraNet/
│
├── README.md
├── .gitignore
│
├── feature_collection/
│   ├── 01_get_amino_acid.r
│   ├── 02_get_domains.r
│   ├── 03_get_phylop.r
│   ├── 04_postprocess.py
│   ├── feature_pipeline.sbatch
│   ├── example_input/
│   │   └── example_guides.csv
│   └── feature_collection_readme.md
│
└── scoring/
    ├── trained_models/
    │   ├── D156R_cds/
    │   ├── KX10_cds/
    │   ├── KX12_cds/
    │   ├── D156R_seq/
    │   ├── KX10_seq/
    │   ├── KX12_seq/
    │
    ├── score_guides.py
    ├── example_input/
    │   └── example_features.csv
    └── scoring_readme.md
└──env_dependencies.md
```

---

## Workflow

### Step 1 — Feature Extraction
```
cd feature_collection/

Rscript 01_get_amino_acid.r raw.csv step1.csv
Rscript 02_get_domains.r    step1.csv step2.csv
Rscript 03_get_phylop.r     step2.csv phyloP.bw step3.csv
python3 04_postprocess.py --input step3.csv --output final_features.csv
```
Produces: **final_features.csv**

### Step 2 — Scoring with HydraNet
```
cd scoring/

python3 score_guides.py   --features_file final_features.csv   --models_dir trained_models/
```
Adds 6 columns:
- `D156R_CDS`, `KX10_CDS`, `KX12_CDS`
- `D156R_Seq`, `KX10_Seq`, `KX12_Seq`

---

## Environments

### R (`cas12_r_env`)
Used for AA projection, domain annotations, PhyloP, ViennaRNA calls.

### Python (`cas12_tf_env`)
Versions:
- Python 3.11.13  
- TensorFlow 2.18  
- Pandas 2.3.2  
- NumPy 2.3.2  
- Scikit-learn 1.7.1  
- Biopython 1.85  

---

## Citation
Please cite **Jeon et al., An Optimized Cas12a Toolkit for Scalable Combinatorial Genetic Screening and Single-Cell Transcriptomics. Submitted.** if used in your work.

