#!/usr/bin/env python3
"""
Predict with a pre-trained LightGBM model on the condensed .npz file,
and save predictions together with the original feature columns.
"""

import argparse, os, pickle as pkl
from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("-n", "--npz", required=True,
                   help="Condensed feature file from condense_training_data.py")
    p.add_argument("-m", "--model", required=True,
                   help="Pickled LightGBM model (.pkl)")
    p.add_argument("-o", "--output-dir", required=True,
                   help="Directory for predictions_with_features.csv")
    p.add_argument("-t", "--threads", type=int, default=os.cpu_count(),
                   help="Threads for LightGBM predict (default: all cores)")
    return p.parse_args()

def _best_feature_names(npz_obj, n_features: int):
    # Try common keys written by pipelines; otherwise fallback
    candidate_keys = ["feature_names", "features", "feat_names", "columns"]
    for k in candidate_keys:
        if k in npz_obj.files:
            names = npz_obj[k]
            # ensure 1D of strings
            try:
                names = names.astype(str).tolist()
            except Exception:
                names = [str(x) for x in names]
            if len(names) == n_features:
                return names
    return [f"feature_{i}" for i in range(n_features)]

def main():
    args = parse_args()

    # 1) Load model
    with open(args.model, "rb") as f:
        model = pkl.load(f)

    # Set threads for prediction
    if hasattr(model, "n_jobs"):
        model.n_jobs = args.threads
    elif hasattr(model, "set_param"):
        model.set_param({"num_threads": args.threads})

    # 2) Load data
    with np.load(args.npz, mmap_mode="r") as z:
        if "X" not in z.files or "files" not in z.files or "labels" not in z.files:
            out = Path(args.output_dir) / "predictions_with_features.csv"
            pd.DataFrame(columns=["file", "true", "pred", "cntrl_score"]).to_csv(out, index=False)
            print(f"No reference TEs found in {args.npz} – wrote empty predictions_with_features.csv")
            return

        X      = z["X"]
        files  = z["files"].astype(str)
        labels = z["labels"].astype(str)
        feat_names = _best_feature_names(z, X.shape[1])

    if X.size == 0:
        out = Path(args.output_dir) / "predictions_with_features.csv"
        pd.DataFrame(columns=["file", "true", "pred", "cntrl_score"]).to_csv(out, index=False)
        print("No samples found – wrote empty predictions_with_features.csv")
        return

    # 3) Predict in chunks (keeps memory steady for huge X)
    CHUNK = 250_000
    preds = np.empty(X.shape[0], dtype=object)  # object to handle string/float/int labels robustly
    for i in tqdm(range(0, X.shape[0], CHUNK), desc="LightGBM predict"):
        preds[i:i+CHUNK] = model.predict(X[i:i+CHUNK])

    # 4) Build output with features
    os.makedirs(args.output_dir, exist_ok=True)
    base_df = pd.DataFrame({
        "file": files,
        "true": labels,
        "pred": preds,
        "cntrl_score": 0,       # keep column for parity with your existing outputs
    })

    # X may be a memmap; build features DataFrame without copying when possible
    features_df = pd.DataFrame(X, columns=feat_names)

    full_df = pd.concat([base_df, features_df], axis=1)

    out_csv = Path(args.output_dir) / "predictions_with_features.csv"
    full_df.to_csv(out_csv, index=False)
    print(f"Wrote predictions with features for {len(full_df):,} samples → {out_csv}")

if __name__ == "__main__":
    main()
