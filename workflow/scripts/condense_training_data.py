#!/usr/bin/env python3
"""
Condense a directory of *.npy feature vectors into a single .npz

* Skips any file whose basename ends with 'readcount.npy'.
* Writes:
      X      – 2-D float32 (n_samples × n_features)
      files  – UTF-8 array of original basenames
      labels – UTF-8 array of true labels (last-but-one field in name)
"""

import argparse, os, re
from glob import glob
from pathlib import Path
import numpy as np
from tqdm import tqdm

READCOUNT_RE = re.compile(r"readcount\.npy$")

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input_dir", required=True,
                   help="Directory containing .npy files (recursed).")
    p.add_argument("-o", "--output-npz", required=True,
                   help="Output .npz path.")
    return p.parse_args()

def main() -> None:
    args = parse_args()

    files, data = [], []
    for f in tqdm(glob(f"{args.input_dir}/**/*.npy", recursive=True),
                  desc="Reading feature vectors"):
        if READCOUNT_RE.search(f):
            continue                       # ← skip read-count files
        arr = np.load(f, mmap_mode="r")    # no copy, fast
        data.append(arr.astype(np.float32, copy=False))
        files.append(Path(f).stem)         # store basename without ".npy"

    if not data:
        raise RuntimeError("No feature vectors found (all were readcount?)")

    # sanity-check & stack
    feat_len = {x.shape for x in data}
    if len(feat_len) != 1:
        raise ValueError(f"Inconsistent vector lengths: {feat_len}")
    X = np.stack(data, axis=0)             # shape (n_samples, n_features)

    labels = np.array([n.split("-")[-2] for n in files], dtype="U")

    np.savez_compressed(args.output_npz, X=X, files=np.array(files, dtype="U"),
                        labels=labels)
    print(f"Wrote {X.shape[0]} samples × {X.shape[1]} features → {args.output_npz}")

if __name__ == "__main__":
    main()
