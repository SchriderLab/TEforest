#!/usr/bin/env python3
"""
Predict with a pre-trained LightGBM model on the condensed .npz file.
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
                   help="Directory for predictions.csv")
    p.add_argument("-t", "--threads", type=int, default=os.cpu_count(),
                   help="Threads for LightGBM predict (default: all cores)")
    return p.parse_args()

def main():
    args = parse_args()

    # 1 . load model
    with open(args.model, "rb") as f:
        model = pkl.load(f)

    # LightGBM sklearn wrapper: honour n_jobs
    if hasattr(model, "n_jobs"):
        model.n_jobs = args.threads
    # Native Booster: set OMP threads
    elif hasattr(model, "set_param"):
        model.set_param({"num_threads": args.threads})

    with np.load(args.npz, mmap_mode="r") as z:
        # if condense wrote no arrays, bail out gracefully
        if "X" not in z.files or "files" not in z.files or "labels" not in z.files:
            out = Path(args.output_dir) / "predictions.csv"
            pd.DataFrame(columns=["file", "true", "pred", "cntrl_score"]) \
              .to_csv(out, index=False)
            print(f"No reference TEs found in {args.npz} – wrote empty predictions.csv")
            return
        X      = z["X"]
        files  = z["files"]
        labels = z["labels"]

    if X.size == 0:
        out = Path(args.output_dir) / "predictions.csv"
        pd.DataFrame(columns=["file", "true", "pred", "cntrl_score"]).to_csv(out, index=False)
        print("No samples found – wrote empty predictions.csv")
        return

    # 3 . predict in chunks if RAM limited (usually not needed now)
    CHUNK = 250_000
    preds = np.empty(X.shape[0], dtype=np.int16)
    for i in tqdm(range(0, X.shape[0], CHUNK), desc="LightGBM predict"):
        preds[i:i+CHUNK] = model.predict(X[i:i+CHUNK])

    # 4 . write output
    df = pd.DataFrame({"file": files, "true": labels,
                       "pred": preds, "cntrl_score": 0})
    os.makedirs(args.output_dir, exist_ok=True)
    df.to_csv(Path(args.output_dir) / "predictions.csv", index=False)
    print(f"Wrote predictions for {len(df):,} samples → {args.output_dir}/predictions.csv")

if __name__ == "__main__":
    main()
