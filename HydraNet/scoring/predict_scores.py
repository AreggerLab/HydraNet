#!/usr/bin/env python3
import os, json, warnings, argparse
import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
tf.get_logger().setLevel("ERROR")

# =========================
# 🐉🐉🐉Encoding helpers
# =========================
NT2I = {'A':0,'C':1,'G':2,'T':3}

def encode_1mer(seq: str):
    s = (seq or "").upper()
    return [NT2I.get(b, 0) for b in s]

def make_kmer_mapping(k: int):
    vocab_array = np.array(np.meshgrid(*[['A','T','C','G']]*k)).T.reshape(-1, k)
    vocab = [''.join(p) for p in vocab_array]
    return {v: i for i, v in enumerate(vocab)}

KMER_MAP_3 = make_kmer_mapping(3)
KMER_MAP_5 = make_kmer_mapping(5)

def encode_kmer(seq: str, k: int, mapping: dict):
    s = (seq or "").upper()
    L = len(s)
    if L < k: return []
    return [mapping.get(s[i:i+k], 0) for i in range(L-k+1)]

def one_hot_23(seq: str):
    """(23,4) one-hot; pads/truncates to 23; unknown -> 0s."""
    TARGET_LEN = 23
    s = (seq or "").upper()
    if len(s) < TARGET_LEN:
        s = s + ("N" * (TARGET_LEN - len(s)))
    elif len(s) > TARGET_LEN:
        s = s[:TARGET_LEN]
    arr = np.zeros((TARGET_LEN, 4), dtype=np.float32)
    for i, ch in enumerate(s):
        j = NT2I.get(ch)
        if j is not None:
            arr[i, j] = 1.0
    return arr

# =========================
# Protospacer extraction:
# =========================
def ensure_protospacer(df: pd.DataFrame) -> pd.DataFrame:
    # make sure the column exists
    if "protospacer_sequence" not in df.columns:
        df["protospacer_sequence"] = pd.NA

    # required columns
    for c in ("guide_sequence", "pam_sequence", "strand"):
        if c not in df.columns:
            raise ValueError(f"Missing required column: {c}")

    # normalize to string
    g = df["guide_sequence"].astype("string").str.upper()
    p = df["pam_sequence"].astype("string").str.upper()
    s = df["strand"].astype("string")

    # which rows still need filling?
    ok = df["protospacer_sequence"].astype("string").str.len().eq(23)
    need = (~ok).fillna(True)

    # 1) minus strand: protospacer = pam_sequence (must be 23)
    m_mask = need & s.eq("-") & p.str.len().eq(23)
    if m_mask.any():
        df.loc[m_mask, "protospacer_sequence"] = p[m_mask]

    # 2) plus strand: trim left PAM from guide, then take 23
    #    require guide long enough to contain (pam + 23)
    plen = p.str.len().fillna(0)
    glen = g.str.len().fillna(0)
    plus_ok_len = glen.ge(plen + 23)
    p_mask = need & s.eq("+") & plus_ok_len
    if p_mask.any():
        df.loc[p_mask, "protospacer_sequence"] = df.loc[p_mask, ["guide_sequence","pam_sequence"]].apply(
            lambda r: str(r["guide_sequence"]).upper()[len(str(r["pam_sequence"]).upper()):][:23],
            axis=1
        )

    return df


# =========================
# 🐉🐉🐉TFSMLayer output normalization
# =========================
def to_logits_array(pred):
    if isinstance(pred, dict):
        pred = pred.get("final_prediction", next(iter(pred.values())))
    elif isinstance(pred, (list, tuple)):
        pred = pred[0]
    t = tf.convert_to_tensor(pred, dtype=tf.float32)
    return t.numpy().ravel()

# =========================
# 🐉🐉🐉Batched inference helpers
# =========================
def run_model_in_batches_4in(model, X_feat, X1, X3, X5, T=1.0, batch_size=256):
    preds = []
    N = X_feat.shape[0]
    for i in range(0, N, batch_size):
        sl = slice(i, i+batch_size)
        raw = model([
            tf.convert_to_tensor(X_feat[sl], dtype=tf.float32),
            tf.convert_to_tensor(X1[sl],   dtype=tf.float32),
            tf.convert_to_tensor(X3[sl],   dtype=tf.float32),
            tf.convert_to_tensor(X5[sl],   dtype=tf.float32),
        ])
        logits = to_logits_array(raw)
        p = 1.0 / (1.0 + np.exp(-logits / float(T)))
        preds.append(p)
    return np.concatenate(preds, axis=0) if preds else np.array([], dtype=np.float32)

def run_model_in_batches_2in(model, X_tab, X_ohe, T=1.0, batch_size=256):
    preds = []
    N = X_tab.shape[0]
    for i in range(0, N, batch_size):
        sl = slice(i, i+batch_size)
        raw = model([
            tf.convert_to_tensor(X_tab[sl], dtype=tf.float32),
            tf.convert_to_tensor(X_ohe[sl], dtype=tf.float32),
        ])
        logits = to_logits_array(raw)
        p = 1.0 / (1.0 + np.exp(-logits / float(T)))
        preds.append(p)
    return np.concatenate(preds, axis=0) if preds else np.array([], dtype=np.float32)

# =========================
# 🐉🐉🐉Main
# =========================
def main():
    ap = argparse.ArgumentParser(description="Dual scoring: exonic (4-input) + intronic/minimal (2-input); writes 6 columns in-place.")
    ap.add_argument("--features_file", required=True)
    ap.add_argument("--models_dir",    required=True)
    args = ap.parse_args()

    print(f"[INFO] Loading: {args.features_file}")
    df_raw = pd.read_csv(args.features_file, low_memory=False, na_values=["NA"])

    # Build (or fix) protospacer exactly per your CSV scheme
    df_raw = ensure_protospacer(df_raw)
    lens = df_raw["protospacer_sequence"].astype("string").str.len().value_counts(dropna=False).sort_index()
    print("[INFO] protospacer length distribution:", dict(lens))

    # Exonic mask
    exon_col = "exon_idx"
    is_exonic = df_raw[exon_col].notna() if exon_col in df_raw.columns else pd.Series(False, index=df_raw.index)

    # =========================
    # 🐉🐉🐉EXONIC models
    # =========================
    FEATURE_COLS_EXONIC = [
        'PhyloP_AA_Up1','PhyloP_AA_Cut','PhyloP_AA_Down1','PhyloP_GuidePlusFlank','Binary_Domain',
        'cut_position_normalized_CDS','GC_Content','Melting_Temperature','RNAfold_MFE','has_TTT',
        'AT_skew','palindromic_score','self_complementarity_score','Cut_Site_Hydrophobicity',
        'Cut_Site_Charge','Cut_Site_Molecular_Weight','Cut_Site_Isoelectric_Point','Cut_Site_Aromaticity',
        'Cut_Site_Instability','Cut_Site_Flexibility','Cut_Site_Helix','Cut_Site_Sheet','Cut_Site_Coil',
        'Cut_Site_BetaTurn','Cut_Site_Solvent_Accessibility','Cut_Site_VDW','Upstream1_AA_Hydrophobicity',
        'Upstream1_AA_Charge','Upstream1_AA_Molecular_Weight','Upstream1_AA_Isoelectric_Point',
        'Upstream1_AA_Aromaticity','Upstream1_AA_Instability','Upstream1_AA_Flexibility',
        'Upstream1_AA_Helix','Upstream1_AA_Sheet','Upstream1_AA_Coil','Upstream1_AA_BetaTurn',
        'Upstream1_AA_Solvent_Accessibility','Upstream1_AA_VDW','Downstream1_AA_Hydrophobicity',
        'Downstream1_AA_Charge','Downstream1_AA_Molecular_Weight','Downstream1_AA_Isoelectric_Point',
        'Downstream1_AA_Aromaticity','Downstream1_AA_Instability','Downstream1_AA_Flexibility',
        'Downstream1_AA_Helix','Downstream1_AA_Sheet','Downstream1_AA_Coil','Downstream1_AA_BetaTurn',
        'Downstream1_AA_Solvent_Accessibility','Downstream1_AA_VDW','Avg_3_Hydrophobicity','Avg_3_Charge',
        'Avg_3_Molecular_Weight','Avg_3_Isoelectric_Point','Avg_3_Aromaticity','Avg_3_Instability',
        'Avg_3_Flexibility','Avg_3_Helix','Avg_3_Sheet','Avg_3_Coil','Avg_3_BetaTurn',
        'Avg_3_Solvent_Accessibility','Avg_3_VDW','PhyloP_Nuc_1','PhyloP_Nuc_2','PhyloP_Nuc_3',
        'PhyloP_Nuc_4','PhyloP_Nuc_5','PhyloP_Nuc_6','PhyloP_Nuc_7','PhyloP_Nuc_8','PhyloP_Nuc_9',
        'PhyloP_Nuc_10','PhyloP_Nuc_11','PhyloP_Nuc_12','PhyloP_Nuc_13','PhyloP_Nuc_14',
        'PhyloP_Nuc_15','PhyloP_Nuc_16','PhyloP_Nuc_17','PhyloP_Nuc_18','PhyloP_Nuc_19',
        'PhyloP_Nuc_20','PhyloP_Nuc_21','PhyloP_Nuc_22','PhyloP_Nuc_23'
    ]
    exonic_models = [
        {'name':'D156R_CDS', 'path':'D156R_cds'},
        {'name':'KX10_CDS',  'path':'KX10_cds'},
        {'name':'KX12_CDS',  'path':'KX12_cds'},
    ]

    avail_ex = [c for c in FEATURE_COLS_EXONIC if c in df_raw.columns]
    work_ex = df_raw.loc[is_exonic, ["protospacer_sequence"] + avail_ex].copy()
    for c in avail_ex:
        work_ex[c] = pd.to_numeric(work_ex[c], errors="coerce")

    mask_len23_ex = work_ex["protospacer_sequence"].astype("string").str.len().eq(23)
    work_ex = work_ex.loc[mask_len23_ex]
    work_ex = work_ex.dropna(subset=avail_ex)  # strict for exonic
    idx_ex = work_ex.index.to_numpy()
    print(f"[INFO] Exonic rows eligible: {len(idx_ex)}")

    if len(idx_ex) > 0:
        X_feat_ex = StandardScaler().fit_transform(work_ex[avail_ex].values)
        seqs_ex = work_ex["protospacer_sequence"].astype(str)
        X1_ex = np.stack(seqs_ex.apply(encode_1mer).values).astype(np.float32)                    # (N,23)
        X3_ex = np.stack(seqs_ex.apply(lambda s: encode_kmer(s,3,KMER_MAP_3)).values).astype(np.float32)  # (N,21)
        X5_ex = np.stack(seqs_ex.apply(lambda s: encode_kmer(s,5,KMER_MAP_5)).values).astype(np.float32)  # (N,19)

        for cfg in exonic_models:
            name = cfg["name"]
            model_path = os.path.join(args.models_dir, cfg["path"])
            print(f"\n[INFO] Exonic model: {name}")
            if not os.path.exists(model_path):
                print(f"[WARN] Missing model path, skipping: {model_path}")
                if name not in df_raw.columns: df_raw[name] = np.nan
                continue
            model = tf.keras.layers.TFSMLayer(model_path, call_endpoint="serve")
            T = 1.0
            try:
                with open(os.path.join(model_path, "temperature.json"), "r") as f:
                    T = float(json.load(f)["temperature"])
            except FileNotFoundError:
                pass
            probs = run_model_in_batches_4in(model, X_feat_ex, X1_ex, X3_ex, X5_ex, T=T, batch_size=256)
            if name not in df_raw.columns: df_raw[name] = np.nan
            df_raw.loc[idx_ex, name] = probs
            print(f"[OK] Wrote '{name}' for {len(idx_ex)} rows.")

    # =========================
    # 🐉🐉🐉INTRONIC / INTERGENIC minimal (2-input; EXACT 8 features)
    # =========================
    FEATURE_COLS_INTRONIC = [
        'GC_Content','Melting_Temperature','RNAfold_MFE','has_TTT',
        'AT_skew','palindromic_score','self_complementarity_score'
    ]
    intron_models = [
        {'name':'D156R_Seq', 'path':'D156R_seq'},
        {'name':'KX10_Seq',  'path':'KX10_seq'},
        {'name':'KX12_Seq',  'path':'KX12_seq'},
    ]

    for c in FEATURE_COLS_INTRONIC:
        if c not in df_raw.columns:
            df_raw[c] = np.nan

    work_in = df_raw.loc[:, ["protospacer_sequence"] + FEATURE_COLS_INTRONIC].copy()
    for c in FEATURE_COLS_INTRONIC:
        work_in[c] = pd.to_numeric(work_in[c], errors="coerce")

    mask_len23_all = work_in["protospacer_sequence"].astype("string").str.len().eq(23)
    work_in = work_in.loc[mask_len23_all]
    idx_all = work_in.index.to_numpy()
    print(f"[INFO] Minimal (all guides) rows encodable: {len(idx_all)}")

    if len(idx_all) > 0:
        # IMPORTANT: no scaling, NaNs preserved
        X_tab = work_in[FEATURE_COLS_INTRONIC].astype(np.float32).values
        seqs_all = work_in["protospacer_sequence"].astype(str)
        X_ohe = np.stack(seqs_all.apply(one_hot_23).values).astype(np.float32)

        for cfg in intron_models:
            name = cfg["name"]
            model_path = os.path.join(args.models_dir, cfg["path"])
            print(f"\n[INFO] Minimal model: {name}")
            if not os.path.exists(model_path):
                print(f"[WARN] Missing model path, skipping: {model_path}")
                if name not in df_raw.columns: df_raw[name] = np.nan
                continue
            model = tf.keras.layers.TFSMLayer(model_path, call_endpoint="serve")
            T = 1.0
            try:
                with open(os.path.join(model_path, "temperature.json"), "r") as f:
                    T = float(json.load(f)["temperature"])
            except FileNotFoundError:
                pass
            probs = run_model_in_batches_2in(model, X_tab, X_ohe, T=T, batch_size=256)
            if name not in df_raw.columns: df_raw[name] = np.nan
            df_raw.loc[idx_all, name] = probs
            print(f"[OK] Wrote '{name}' for {len(idx_all)} rows.")

    # =========================
    # 🐉🐉🐉Save back in-place (The scores will be saved in the same csv as feature csv i.e. the csvs would be over-written with score cols)
    # =========================
    print(f"\n[INFO] Saving: {args.features_file}")
    df_raw.to_csv(args.features_file, index=False, na_rep="NA")
    print("[DONE] Dual scoring complete.")

if __name__ == "__main__":
    main()

