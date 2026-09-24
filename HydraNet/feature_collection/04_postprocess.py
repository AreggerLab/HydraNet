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
import sys # Import sys to exit

# === CONFIG ===
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument(
    "--gtf",
    default="/path/Homo_sapiens.GRCh38.104.gtf",
    help="Path to the GRCh38 GTF file",
)
args = parser.parse_args()

input_csv = args.input
output_csv = args.output
gtf_path = args.gtf

# Columns from upstream R scripts
PROTOSPACER_COL = 'protospacer_sequence'
LOCATION_COL = 'exon_id'

# === RNAfold detection ===
CONDA_ENV_PATH = os.environ.get("CONDA_PREFIX")
RNAFOLD_PATH = shutil.which("RNAfold", path=os.path.join(CONDA_ENV_PATH, "bin") if CONDA_ENV_PATH else None)

if not RNAFOLD_PATH:
    print("⚠️ RNAfold not found. MFE column will be empty.")
else:
    print(f"✅ Found RNAfold at: {RNAFOLD_PATH}")

# === Load Data ===
print(f"📦 Loading input data from: {input_csv}")
try:
    df = pd.read_csv(input_csv, na_values="NA", low_memory=False)
    print(f"🔹 Loaded {len(df)} rows.")
except FileNotFoundError:
    raise FileNotFoundError(f"FATAL: The input file from the R script was not found: {input_csv}")

# =================================================================================
# FIX: Standardize and VALIDATE column names.
# =================================================================================
print("🔹 Standardizing and validating column names...")
if 'seqnames' in df.columns and 'chromosome' not in df.columns:
    df.rename(columns={'seqnames': 'chromosome'}, inplace=True)
    print("   - Renamed 'seqnames' to 'chromosome'")

if 'cut_site' in df.columns:
    # Preserve the original cut_site column and derive a numeric cut_pos from it.
    # Script 5 had to backfill cut_site because Script 4 renamed it away.
    df['cut_pos'] = pd.to_numeric(df['cut_site'], errors='coerce')
    print("   - Derived numeric 'cut_pos' from 'cut_site' (kept 'cut_site')")
elif 'cut_pos' in df.columns:
    df['cut_pos'] = pd.to_numeric(df['cut_pos'], errors='coerce')

# --- START: NEW VALIDATION BLOCK ---
required_cols = ['cut_pos', 'chromosome']
missing_cols = [col for col in required_cols if col not in df.columns]
if missing_cols:
    print(f"❌ FATAL: Input CSV is missing required columns: {missing_cols}", file=sys.stderr)
    sys.exit(1) # Exit with an error code

if 'gene_id' not in df.columns and 'gene_name' not in df.columns:
    print("❌ FATAL: Input CSV must have at least one of 'gene_id' or 'gene_name'.", file=sys.stderr)
    sys.exit(1)
# --- END: NEW VALIDATION BLOCK ---


# Sanity check: require protospacer from R step
if PROTOSPACER_COL not in df.columns:
    raise ValueError(f"Missing '{PROTOSPACER_COL}' from previous R step. Do not recompute here.")

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
        match = re.search(r'\(\s*(-?\d+(?:\.\d+)?)\s*\)\s*$', result.stdout.strip())
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

# === GTF attr parser (keeps repeated keys like tag) ===
def parse_gtf_attrs(attr_str: str):
    out = defaultdict(list)
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr_str):
        out[m.group(1)].append(m.group(2))
    return out

def chrom_norm(x):
    s = str(x)
    s = s.replace("chr", "")
    # keep Ensembl style: "MT" not "M"
    if s == "M":
        s = "MT"
    return s

def has_tag_contains(tid: str, needle: str, transcript_tags: dict) -> bool:
    tags = transcript_tags.get(tid, [])
    needle = needle.lower()
    return any(needle in t.lower() for t in tags)

# === GTF parsing for CDS-relative position ===
print("📍 Computing CDS-relative cut positions...")
df["cut_position_normalized_CDS"] = pd.NA
is_exonic = df.get(LOCATION_COL).notna() & (df.get(LOCATION_COL) != 'NA') if LOCATION_COL in df.columns else pd.Series([False]*len(df))
df_exonic = df[is_exonic].copy()

if not df_exonic.empty:
    print("   Parsing GTF (transcript + CDS features)...")
    # index CDS by both gene_name and gene_id to avoid lookup mismatch
    gene_to_cds_by_name = defaultdict(lambda: defaultdict(list))
    gene_to_cds_by_id   = defaultdict(lambda: defaultdict(list))
    # transcript_id -> list of tags (from transcript feature lines)
    transcript_tags = defaultdict(list)

    with open(gtf_path, "r") as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feature = cols[2]
            attrs = parse_gtf_attrs(cols[8])
            tid = (attrs.get("transcript_id") or [None])[0]
            gid = (attrs.get("gene_id") or [None])[0]
            gname = (attrs.get("gene_name") or [None])[0]
            tags = attrs.get("tag") or []

            # Capture transcript metadata from transcript and CDS records.
            # APPRIS tags and CCDS identifiers can appear in GTF attributes.
            if tid:
                if tags:
                    transcript_tags[tid].extend(tags)
                ccds_id = (attrs.get("ccds_id") or [None])[0]
                if ccds_id:
                    transcript_tags[tid].append(ccds_id)

            if feature == "transcript":
                continue

            if feature != "CDS" or not tid:
                continue

            # store CDS blocks
            block = {
                "start": int(cols[3]),
                "end": int(cols[4]),
                "strand": cols[6],
                "chrom": cols[0]
            }
            if gname:
                gene_to_cds_by_name[gname][tid].append(block)
            if gid:
                gene_to_cds_by_id[gid][tid].append(block)

    normalized_positions = []
    for _, row in tqdm(df_exonic.iterrows(), total=len(df_exonic), desc="   RelPos"):
        cut = row.get("cut_pos")
        chrom = chrom_norm(row.get("chromosome"))
        # gene lookup: prefer gene_id (stable), then gene_name
        gene_id = row.get("gene_id")
        gene_name = row.get("gene_name")
        transcripts = {}

        if isinstance(gene_id, str) and gene_id in gene_to_cds_by_id:
            transcripts = gene_to_cds_by_id[gene_id]
        elif isinstance(gene_name, str) and gene_name in gene_to_cds_by_name:
            transcripts = gene_to_cds_by_name[gene_name]
        else:
            normalized_positions.append(pd.NA)
            continue

        if pd.isna(cut) or not transcripts:
            normalized_positions.append(pd.NA)
            continue

        # Canonical transcript selection matched to the Script 5 repair:
        # prefer APPRIS principal, then CCDS, otherwise the longest CDS.
        canonical_tid = next(
            (t for t in transcripts if has_tag_contains(t, "appris_principal", transcript_tags)),
            None,
        )
        if not canonical_tid:
            canonical_tid = next(
                (t for t in transcripts if has_tag_contains(t, "CCDS", transcript_tags)),
                None,
            )
        if not canonical_tid and transcripts:
            canonical_tid = max(
                transcripts,
                key=lambda t: sum(
                    e["end"] - e["start"] + 1 for e in transcripts.get(t, [])
                ),
            )

        if not canonical_tid:
            normalized_positions.append(pd.NA)
            continue

        blocks = [e for e in transcripts[canonical_tid] if chrom_norm(e["chrom"]) == chrom]
        if not blocks:
            normalized_positions.append(pd.NA)
            continue

        strand = blocks[0]["strand"]
        if strand == "+":
            blocks = sorted(blocks, key=lambda x: x["start"])
        else:
            blocks = sorted(blocks, key=lambda x: x["end"], reverse=True)

        transcript_pos = 0
        rel_cut = None
        for b in blocks:
            if b["start"] <= cut <= b["end"]:
                if strand == "+":
                    rel_cut = transcript_pos + (cut - b["start"] + 1)
                else:
                    rel_cut = transcript_pos + (b["end"] - cut + 1)
                break
            transcript_pos += (b["end"] - b["start"] + 1)
        
        total_len = sum(b["end"] - b["start"] + 1 for b in blocks)

        if rel_cut is None or total_len <= 0:
            normalized_positions.append(pd.NA)
        else:
            normalized_positions.append(round(rel_cut / total_len, 5))

    df.loc[is_exonic, "cut_position_normalized_CDS"] = normalized_positions

# === Sequence-based features (protospacer) ===
print("🌡️ Calculating GC content and melting temperature...")
df["GC_Content"] = df[PROTOSPACER_COL].apply(
    lambda s: round(100 * (s.count("G")+s.count("C")) / len(s), 2) if isinstance(s, str) and len(s) > 0 else pd.NA
)
df["Melting_Temperature"] = df[PROTOSPACER_COL].apply(
    lambda s: 2*(s.count("A")+s.count("T")) + 4*(s.count("G")+s.count("C")) if isinstance(s, str) else pd.NA
)

print("🧬 Calculating RNAfold MFE...")
df["RNAfold_MFE"] = df[PROTOSPACER_COL].progress_apply(compute_mfe)

print("🔢 Computing k-mer and palindromic features...")
df = pd.concat([df, df[PROTOSPACER_COL].progress_apply(sequence_features)], axis=1)

# --- biochemical properties ---
print("🧪 AA biochemical properties…")
kd_hydro = {'I':4.5,'V':4.2,'L':3.8,'F':2.8,'C':2.5,'M':1.9,'A':1.8,'G':-0.4,'T':-0.7,'S':-0.8,'W':-0.9,'Y':-1.3,'P':-1.6,'H':-3.2,'E':-3.5,'Q':-3.5,'D':-3.5,'N':-3.5,'K':-3.9,'R':-4.5}
vdw_vol  = {'A':67,'R':148,'N':96,'D':91,'C':86,'E':109,'Q':114,'G':48,'H':118,'I':124,'L':124,'K':135,'M':124,'F':135,'P':90,'S':73,'T':93,'W':163,'Y':141,'V':105}
FLEX_SCALE = {
    'A':0.357,'R':0.529,'N':0.463,'D':0.511,'C':0.346,'Q':0.493,'E':0.497,'G':0.544,'H':0.323,'I':0.462,
    'L':0.365,'K':0.466,'M':0.295,'F':0.314,'P':0.509,'S':0.507,'T':0.444,'W':0.305,'Y':0.420,'V':0.386
}

def compute_props(seq: str):
    n_na = 13
    if pd.isna(seq) or not isinstance(seq, str):
        return [pd.NA]*n_na
    
    seq = seq.replace("*", "")
    if not seq or not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq):
        return [pd.NA]*n_na
    
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
    except Exception:
        return [pd.NA]*n_na

for c in ["aa_cut", "aa_up", "aa_down"]:
    if c not in df.columns:
        df[c] = pd.NA

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
    
    if hydro_col not in df_in.columns or coil_col not in df_in.columns:
        return df_in

    hydro = pd.to_numeric(df_in[hydro_col], errors="coerce")
    coil  = pd.to_numeric(df_in[coil_col], errors="coerce")
    
    df_in[beta_col] = (coil + np.maximum(0.0, -hydro/10.0)).round(4)
    df_in[solv_col] = (np.maximum(0.0, 1.0 - hydro/4.0)).round(4)
    return df_in

for pfx in PREFIXES_FOR_FILL:
    df = fill_beta_solvent(df, pfx)

print(f"💾 Saving: {output_csv}")
df.to_csv(output_csv, index=False, na_rep="NA")
print("✅ Done.")
