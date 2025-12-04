#!/usr/bin/env python3
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import defaultdict
import subprocess
from tqdm import tqdm
import re
import shutil
import os
import argparse

# === CONFIG ===
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

input_csv = args.input
output_csv = args.output

gtf_path = "PATH/gencode.vM10.annotation.gtf"
#gtf_path = "PATH/Homo_sapiens.GRCh38.104.gtf"

# Columns from upstream R scripts
PROTOSPACER_COL = 'protospacer_sequence'
LOCATION_COL = 'exon_id'

# === RNAfold detection ===
CONDA_ENV_PATH = os.environ.get("CONDA_PREFIX")
RNAFOLD_PATH = shutil.which("RNAfold", path=os.path.join(CONDA_ENV_PATH, "bin") if CONDA_ENV_PATH else None)
if not RNAFOLD_PATH:
    print("!!! RNAfold not found. MFE column will be empty.")
else:
    print(f"🐉🐉🐉 Found RNAfold at: {RNAFOLD_PATH}")

# === Load Data ===
print(f" Loading input data from: {input_csv}")
try:
    df = pd.read_csv(input_csv, na_values="NA", low_memory=False)
    print(f"🐉🐉🐉 Loaded {len(df)} rows.")
except FileNotFoundError:
    raise FileNotFoundError(f"FATAL: The input file from the R script was not found: {input_csv}")


print(" Standardizing column names...")
rename_map = {
    'seqnames': 'chromosome',
    #'cut_site': 'cut_pos'
}
# Select only columns that actually exist in the dataframe to avoid errors
cols_to_rename = {k: v for k, v in rename_map.items() if k in df.columns}
if cols_to_rename:
    df.rename(columns=cols_to_rename, inplace=True)
    print(f"   - Renamed columns: {cols_to_rename}")

# Sanity check: require protospacer from R step
if PROTOSPACER_COL not in df.columns:
    raise ValueError(f"Missing '{PROTOSPACER_COL}' from previous R step.")

tqdm.pandas()

# === Helper: RNAfold MFE ===
def compute_mfe(seq):
    if not RNAFOLD_PATH or pd.isna(seq) or not isinstance(seq, str) or len(seq) < 4:
        return pd.NA
    seq_rna = seq.upper().replace("T", "U").replace("N", "")
    if not seq_rna:
        return pd.NA
    try:
        result = subprocess.run(
            [RNAFOLD_PATH, "--noPS"],
            input=seq_rna,
            capture_output=True,
            text=True,
            check=True,
            timeout=15
        )
        match = re.search(r'\(\s*(-?\d+\.\d+)\s*\)$', result.stdout.strip())
        return float(match.group(1)) if match else pd.NA
    except Exception:
        return pd.NA

# === Helper: Sequence Features ===
def sequence_features(seq):
    if pd.isna(seq) or not isinstance(seq, str) or len(seq) == 0:
        return pd.Series({
            'has_TTT': pd.NA,
            'AT_skew': pd.NA,
            'palindromic_score': pd.NA,
            'self_complementarity_score': pd.NA
        })
    seq = seq.upper()
    has_TTT = int('TTT' in seq)
    a, t = seq.count("A"), seq.count("T")
    at_skew = (a - t) / (a + t) if (a + t) > 0 else 0
    rev = str(Seq(seq).reverse_complement())
    pal_score = sum(seq[i] == rev[i] for i in range(len(seq))) / len(seq)
    max_comp = 0
    for k in range(4, 12):
        for i in range(len(seq) - k + 1):
            sub_seq = seq[i:i+k]
            if str(Seq(sub_seq).reverse_complement()) in seq:
                max_comp = max(max_comp, k)
    return pd.Series({
        'has_TTT': has_TTT,
        'AT_skew': at_skew,
        'palindromic_score': pal_score,
        'self_complementarity_score': max_comp
    })

# === GTF parsing for CDS-relative position ===
print("🐉🐉🐉 Computing CDS-relative cut positions...")
df["cut_position_normalized_CDS"] = pd.NA
is_exonic = df.get(LOCATION_COL).notna() & (df.get(LOCATION_COL) != 'NA') if LOCATION_COL in df.columns else pd.Series([False]*len(df))
df_exonic = df[is_exonic].copy()

if not df_exonic.empty:
    print("   Parsing GTF...")
    gene_to_cds = defaultdict(lambda: defaultdict(list))
    transcript_tags = {}
    with open(gtf_path, "r") as gtf:
        for line in gtf:
            if line.startswith("#") or "\tCDS\t" not in line:
                continue
            cols = line.strip().split("\t")
            attr = {
                k.strip(): v.strip().strip('"')
                for entry in cols[8].split(";")
                if entry.strip() and ' ' in entry
                for k, v in [entry.strip().split(" ", 1)]
            }
            gene = attr.get("gene_name")
            tid  = attr.get("transcript_id")
            if gene and tid:
                transcript_tags[tid] = cols[8]
                gene_to_cds[gene][tid].append({
                    "start": int(cols[3]),
                    "end": int(cols[4]),
                    "strand": cols[6],
                    "chrom": cols[0]
                })

    normalized_positions = []
    # Use the standardized column names 'gene_name', 'cut_pos', 'chromosome'
    for _, row in tqdm(df_exonic.iterrows(), total=len(df_exonic), desc="   RelPos"):
        gene = row["gene_name"]
        cut  = row["cut_pos"]
        chrom = str(row["chromosome"]).replace("chr", "")

        # Robust checks: make sure gene is a string and cut is a number
        if not isinstance(gene, str) or pd.isna(cut):
            normalized_positions.append(pd.NA)
            continue

        if gene not in gene_to_cds:
            normalized_positions.append(pd.NA)
            continue

        transcripts = gene_to_cds.get(gene, {})
        if not transcripts:
            normalized_positions.append(pd.NA)
            continue

        canonical_tid = next((t for t in transcripts if "appris_principal" in transcript_tags[t]), None)
        if not canonical_tid:
            ccds = [t for t in transcripts if "CCDS" in transcript_tags.get(t, "")]
            canonical_tid = ccds[0] if ccds else max(
                transcripts,
                key=lambda t: sum(e["end"] - e["start"] + 1 for e in transcripts[t])
            )

        blocks = [e for e in transcripts[canonical_tid] if e["chrom"].replace("chr", "") == chrom]
        if not blocks:
            normalized_positions.append(pd.NA)
            continue

        blocks = sorted(blocks, key=lambda x: x["start"] if blocks[0]["strand"] == "+" else -x["end"])
        transcript_pos = 0
        rel_cut = None
        for b in blocks:
            if b["start"] <= cut <= b["end"]:
                if b["strand"] == "+":
                    rel_cut = transcript_pos + (cut - b["start"] + 1)
                else:
                    rel_cut = transcript_pos + (b["end"] - cut + 1)
                break
            transcript_pos += (b["end"] - b["start"] + 1)

        total_len = sum(b["end"] - b["start"] + 1 for b in blocks)
        normalized_positions.append(
            round(rel_cut / total_len, 5) if rel_cut and total_len > 0 else pd.NA
        )


    df.loc[is_exonic, "cut_position_normalized_CDS"] = normalized_positions

# === Sequence-based features (protospacer) ===
print("🐉🐉🐉 Calculating GC content and melting temperature...")
df["GC_Content"] = df[PROTOSPACER_COL].apply(
    lambda s: round(100 * (s.count("G")+s.count("C")) / len(s), 2) if isinstance(s, str) and len(s) > 0 else pd.NA
)
df["Melting_Temperature"] = df[PROTOSPACER_COL].apply(
    lambda s: 2*(s.count("A")+s.count("T")) + 4*(s.count("G")+s.count("C")) if isinstance(s, str) else pd.NA
)

print("🐉🐉🐉 Calculating RNAfold MFE...")
df["RNAfold_MFE"] = df[PROTOSPACER_COL].progress_apply(compute_mfe)

print("🐉🐉🐉 Computing k-mer and palindromic features...")
df = pd.concat([df, df[PROTOSPACER_COL].progress_apply(sequence_features)], axis=1)

# --- biochemical properties ---
print("🐉🐉🐉 AA biochemical properties…")

# Scales
kd_hydro = {'I':4.5,'V':4.2,'L':3.8,'F':2.8,'C':2.5,'M':1.9,'A':1.8,'G':-0.4,'T':-0.7,'S':-0.8,'W':-0.9,'Y':-1.3,'P':-1.6,'H':-3.2,'E':-3.5,'Q':-3.5,'D':-3.5,'N':-3.5,'K':-3.9,'R':-4.5}
vdw_vol  = {'A':67,'R':148,'N':96,'D':91,'C':86,'E':109,'Q':114,'G':48,'H':118,'I':124,'L':124,'K':135,'M':124,'F':135,'P':90,'S':73,'T':93,'W':163,'Y':141,'V':105}
FLEX_SCALE = {
    'A':0.357,'R':0.529,'N':0.463,'D':0.511,'C':0.346,'Q':0.493,'E':0.497,'G':0.544,'H':0.323,'I':0.462,
    'L':0.365,'K':0.466,'M':0.295,'F':0.314,'P':0.509,'S':0.507,'T':0.444,'W':0.305,'Y':0.420,'V':0.386
}

def compute_props(seq: str):
    n_na = 13
    if pd.isna(seq) or not isinstance(seq, str): return [pd.NA]*n_na
    seq = seq.replace("*","")
    if not seq or not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq): return [pd.NA]*n_na
    try:
        pa = ProteinAnalysis(seq)
        hydro  = sum(kd_hydro.get(aa, 0) for aa in seq) / len(seq)
        charge = pa.charge_at_pH(7.0)
        mw     = pa.molecular_weight()
        iso    = pa.isoelectric_point()
        aroma  = pa.aromaticity()
        instab = pa.instability_index()
        flex_vals = pa.flexibility()
        flex = float(np.nanmean(flex_vals)) if len(flex_vals) > 0 else float(np.nanmean([FLEX_SCALE.get(aa, np.nan) for aa in seq]))
        helix, sheet, coil = pa.secondary_structure_fraction()
        vdw = sum(vdw_vol.get(aa, 0) for aa in seq) / len(seq)
        vec = [hydro, charge, mw, iso, aroma, instab, flex, helix, sheet, coil, pd.NA, pd.NA, vdw]
        return [round(float(x), 4) if isinstance(x, (int, float, np.floating)) and not pd.isna(x) else pd.NA for x in vec]
    except Exception: return [pd.NA]*n_na

for c in ["aa_cut", "aa_up", "aa_down"]:
    if c not in df.columns: df[c] = pd.NA

df["aa_up1"]   = df["aa_up"].astype("string").str[-1]
df["aa_down1"] = df["aa_down"].astype("string").str[0]
df["aa_avg3"]  = df["aa_cut"].fillna('') + df["aa_up1"].fillna('') + df["aa_down1"].fillna('')

prefix_map = {"Cut": "aa_cut", "Up1": "aa_up1", "Down1": "aa_down1", "Avg3": "aa_avg3"}
label_map = {"Cut": "Cut_Site", "Up1": "Upstream1_AA", "Down1": "Downstream1_AA", "Avg3": "Avg_3"}
prop_names = ["Hydrophobicity","Charge","Molecular_Weight","Isoelectric_Point","Aromaticity","Instability","Flexibility","Helix","Sheet","Coil","BetaTurn","Solvent_Accessibility","VDW"]

for prefix, col in prefix_map.items():
    props = df[col].progress_apply(compute_props)
    df_props = pd.DataFrame(props.tolist(), index=df.index, columns=[f"{label_map[prefix]}_{p}" for p in prop_names])
    df = pd.concat([df, df_props], axis=1)

PREFIXES_FOR_FILL = ["Cut_Site", "Upstream1_AA", "Downstream1_AA", "Avg_3"]
def fill_beta_solvent(df_in: pd.DataFrame, prefix: str) -> pd.DataFrame:
    hydro_col, coil_col, beta_col, solv_col = f"{prefix}_Hydrophobicity", f"{prefix}_Coil", f"{prefix}_BetaTurn", f"{prefix}_Solvent_Accessibility"
    if hydro_col not in df_in.columns or coil_col not in df_in.columns: return df_in
    hydro, coil = pd.to_numeric(df_in[hydro_col], errors="coerce"), pd.to_numeric(df_in[coil_col], errors="coerce")
    df_in[beta_col] = (coil + np.maximum(0.0, -hydro/10.0)).round(4)
    df_in[solv_col] = (np.maximum(0.0, 1.0 - hydro/4.0)).round(4)
    return df_in

for pfx in PREFIXES_FOR_FILL:
    df = fill_beta_solvent(df, pfx)

print(f"🐉🐉🐉 Saving: {output_csv}")
df.to_csv(output_csv, index=False, na_rep="NA")
print("🐉🐉🐉 Done.")
