# 🐉🐉🐉 Feature Collection Pipeline

This folder contains the full HydraNet feature extraction workflow. It transforms raw Cas12a guide coordinates and protospacer sequences into a complete machine-learning feature table that is compatible with HydraNet scoring.

The pipeline consists of four scripts:

```
01_get_amino_acid.r     → Gene mapping, exon overlap, amino-acid context
02_get_domains.r        → PFAM and InterPro domain overlap
03_get_phylop.r         → PhyloP conservation scores (nucleotide + codon)
04_postprocess.py       → RNAfold MFE, biochemical AA features, GC %, k-mers
```

Each script consumes the output of the previous script.

---

# 🐉🐉🐉 Dependencies (Exact Versions)

## 🐉🐉🐉 R Environment (`cas12_r_env`)
```
r-base 4.4.3
AnnotationHub 3.14.0
ensembldb 2.30.0
GenomicRanges 1.58.0
GenomeInfoDb 1.42.0
rtracklayer 1.66.0
biomaRt 2.62.0
dplyr 1.1.4
readr 2.1.5
tidyr 1.3.1
purrr 1.1.0
stringr 1.5.2
```

## 🐉🐉🐉 Python Environment (`cas12_tf_env`)
```
python 3.11.13
pandas 2.3.2
numpy 2.3.2
biopython 1.85
tqdm 4.67.1
```

```viennarna 2.4.7```

Standard modules (no install):
```
os, re, shutil, subprocess, argparse, collections
```

## 🐉🐉🐉 External Required Files
```
gencode.vM10.annotation.gtf   OR   Homo_sapiens.GRCh38.104.gtf
phyloP bigWig (matching genome build)
```

---

# 1. 🐲🐲🐲 Required Input Format

Your input CSV must contain these mandatory columns:

| Column | Description |
|--------|-------------|
| seqnames | Chromosome name |
| cut_site | Genomic cut coordinate |
| protospacer_sequence | 23-nt protospacer |
| protospacer_start | Genomic start |
| protospacer_end | Genomic end |

Example:
```
seqnames,cut_site,protospacer_sequence,protospacer_start,protospacer_end
chr1,12345678,TTTAGGCTTAGATTGGGCGAGGC,12345656,12345678
```

Optional but recommended:
```
guide_id, strand, gene_name, exon_id
```

---

# 2. 🐲🐲🐲 Running the Pipeline

```
Rscript 01_get_amino_acid.r raw_guides.csv step1_amino.csv
Rscript 02_get_domains.r step1_amino.csv step2_domains.csv
Rscript 03_get_phylop.r step2_domains.csv phyloP.bw step3_phylo.csv
python3 04_postprocess.py --input step3_phylo.csv --output final_features.csv
```
Refer to: feature_pipeline.sbatch
---

# 3. 🐲🐲🐲 Example Input

Place inside:
```
feature_collection/example_input/example_guides.csv
```
# 4. 🐲🐲🐲 Example Execution sbatch
feature_pipeline.sbatch