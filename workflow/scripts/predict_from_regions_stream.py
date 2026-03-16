#!/usr/bin/env python3
"""
Build feature vectors directly from candidate regions and predict on the fly.

This avoids writing per-region .npy feature vectors and condensed .npz files.
"""

import argparse
import multiprocessing as mp
import os
import pickle as pkl
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import bam_to_fvec_cov_normalized_rescale_regions_majorchanges as nonref_fvec
import bam_to_fvec_fullyupdated_forreference as ref_fvec


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stream feature-vector generation and inference."
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        help="Input candidate CSV (Sample,Chrom,Ref_begin,Ref_end,Class,TE).",
    )
    parser.add_argument(
        "-bd",
        "--bam_dir",
        required=True,
        help="Directory containing {sample}.bam files.",
    )
    parser.add_argument(
        "-tebd",
        "--te_bam_dir",
        required=True,
        help="Directory containing {sample}/{TE}_to_ISO1.bam files.",
    )
    parser.add_argument(
        "-m",
        "--model",
        required=True,
        help="Pickled model path (.pkl).",
    )
    parser.add_argument(
        "-o",
        "--output_csv",
        required=True,
        help="Output predictions CSV path.",
    )
    parser.add_argument(
        "--mode",
        choices=("nonreference", "reference"),
        default="nonreference",
        help="Feature-vector backend to use.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=os.cpu_count() or 1,
        help="Total prediction threads budget.",
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=1,
        help="Region-worker processes (1 disables process parallelism).",
    )
    parser.add_argument(
        "--batch_size",
        type=int,
        default=256,
        help="Regions to accumulate before predict+flush.",
    )
    return parser.parse_args()


def load_model(path: str, threads: int):
    with open(path, "rb") as f:
        model = pkl.load(f)

    if hasattr(model, "n_jobs"):
        model.n_jobs = threads
    elif hasattr(model, "set_param"):
        model.set_param({"num_threads": threads})
    return model


def get_backend(mode: str):
    return nonref_fvec if mode == "nonreference" else ref_fvec


def ensure_parent(path: str) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)


def init_output(path: str) -> None:
    pd.DataFrame(columns=["file", "true", "pred", "cntrl_score"]).to_csv(path, index=False)


def flush_batch(model, out_csv, files, labels, vectors) -> int:
    if not vectors:
        return 0

    X = np.stack(vectors, axis=0)
    preds = model.predict(X)
    out_df = pd.DataFrame(
        {
            "file": files,
            "true": labels,
            "pred": preds,
            "cntrl_score": 0,
        }
    )
    out_df.to_csv(out_csv, mode="a", header=False, index=False)

    files.clear()
    labels.clear()
    vectors.clear()
    return len(out_df)


def stream_regions(
    regions: pd.DataFrame,
    mode: str,
    bam_dir: str,
    te_bam_dir: str,
    model,
    out_csv: str,
    batch_size: int,
):
    backend = get_backend(mode)

    files = []
    labels = []
    vectors = []
    written = 0
    processed = 0
    svm_cache = {}

    for idx, region in regions.iterrows():
        processed += 1
        sample = str(region["Sample"])
        main_bam = os.path.join(bam_dir, sample + ".bam")
        te_bam = os.path.join(te_bam_dir, sample, f"{region['TE']}_to_ISO1.bam")

        if sample not in svm_cache:
            svm_cache[sample] = backend.SVMaker(main_bam, output_dir=None)
        svm = svm_cache[sample]

        base_result = svm.create_region(region, main_bam, savefile=False)
        if base_result is None:
            continue
        try:
            base_vec, _read_counts = base_result
        except Exception:
            continue

        te_vec = svm.update_region(region, te_bam, savefile=False)
        if te_vec is None:
            continue

        full_vec = np.concatenate((base_vec, te_vec), axis=0).astype(np.float32, copy=False)
        if full_vec.size == 0:
            continue

        file_id = svm._set_outfile_name(region)
        label = file_id.split("-")[-2]
        files.append(file_id)
        labels.append(label)
        vectors.append(full_vec)

        if len(vectors) >= batch_size:
            written += flush_batch(model, out_csv, files, labels, vectors)

        if idx % 500 == 0:
            import gc

            gc.collect()

    written += flush_batch(model, out_csv, files, labels, vectors)
    return processed, written


def run_worker(
    chunk: pd.DataFrame,
    mode: str,
    bam_dir: str,
    te_bam_dir: str,
    model_path: str,
    model_threads: int,
    batch_size: int,
    out_csv: str,
):
    model = load_model(model_path, model_threads)
    init_output(out_csv)
    processed, written = stream_regions(
        regions=chunk,
        mode=mode,
        bam_dir=bam_dir,
        te_bam_dir=te_bam_dir,
        model=model,
        out_csv=out_csv,
        batch_size=batch_size,
    )
    return processed, written, out_csv


def main() -> None:
    args = parse_args()

    regions = pd.read_csv(args.input_file)
    required = {"Sample", "Chrom", "Ref_begin", "Ref_end", "Class", "TE"}
    missing = required.difference(regions.columns)
    if missing:
        raise ValueError(f"Input CSV missing required columns: {sorted(missing)}")

    ensure_parent(args.output_csv)

    n_proc = max(1, int(args.processes))
    n_proc = min(n_proc, len(regions)) if len(regions) else 1

    if n_proc == 1:
        model = load_model(args.model, args.threads)
        init_output(args.output_csv)
        processed, written = stream_regions(
            regions=regions,
            mode=args.mode,
            bam_dir=args.bam_dir,
            te_bam_dir=args.te_bam_dir,
            model=model,
            out_csv=args.output_csv,
            batch_size=args.batch_size,
        )
        print(
            "Stream inference complete: "
            f"processed={processed}, predicted={written}, output={args.output_csv}, "
            f"processes=1, model_threads={max(1, int(args.threads))}"
        )
        return

    split_indices = [idx for idx in np.array_split(regions.index.to_numpy(), n_proc) if len(idx)]
    chunks = [regions.loc[idx].copy() for idx in split_indices]
    n_proc = max(1, len(chunks))
    model_threads_per_proc = max(1, int(args.threads) // n_proc)

    out_path = Path(args.output_csv).resolve()
    tmp_dir = out_path.parent / f".stream_predict_tmp_{out_path.stem}_{os.getpid()}"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    jobs = []
    for i, chunk in enumerate(chunks):
        chunk_out = tmp_dir / f"chunk_{i}.csv"
        jobs.append(
            (
                chunk,
                args.mode,
                args.bam_dir,
                args.te_bam_dir,
                args.model,
                model_threads_per_proc,
                args.batch_size,
                str(chunk_out),
            )
        )

    total_processed = 0
    total_written = 0
    init_output(args.output_csv)

    try:
        with mp.Pool(n_proc) as pool:
            results = pool.starmap(run_worker, jobs)

        for processed, written, part_csv in results:
            total_processed += processed
            total_written += written
            part_df = pd.read_csv(part_csv)
            if not part_df.empty:
                part_df.to_csv(args.output_csv, mode="a", header=False, index=False)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    print(
        "Stream inference complete: "
        f"processed={total_processed}, predicted={total_written}, output={args.output_csv}, "
        f"processes={n_proc}, model_threads_per_process={model_threads_per_proc}"
    )


if __name__ == "__main__":
    main()
