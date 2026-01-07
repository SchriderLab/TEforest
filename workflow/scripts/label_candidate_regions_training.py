#!/usr/bin/env python3

import argparse
import sys
from typing import Dict, List, Tuple, Optional

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Label TEforest candidate regions using a truth BED."
    )
    parser.add_argument(
        "--candidates",
        required=True,
        help="Candidate regions CSV from process_candidate_regions.r",
    )
    parser.add_argument(
        "--truth",
        required=True,
        help="Truth BED with columns: seqnames, start, end, zygosity, sample[, TE].",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV with Class labels updated.",
    )
    parser.add_argument(
        "--label-mode",
        choices=["classifier", "regressor"],
        default="classifier",
        help="Classifier maps 1/0.5/0 to 2/1/0; regressor keeps numeric zygosity.",
    )
    return parser.parse_args()


def has_header(first_row: List[str]) -> bool:
    expected = {"seqnames", "chrom", "chr", "start", "end", "zygosity", "sample", "te"}
    return any(cell.strip().lower() in expected for cell in first_row)


def read_truth_bed(path: str) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8") as fh:
        first = fh.readline().rstrip("\n").split("\t")

    if has_header(first):
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path, sep="\t", header=None)
        cols = ["seqnames", "start", "end", "zygosity", "sample", "te"]
        df.columns = cols[: df.shape[1]]

    col_map = {}
    for col in df.columns:
        lower = col.strip().lower()
        if lower in {"seqnames", "chrom", "chr"}:
            col_map[col] = "seqnames"
        elif lower == "start":
            col_map[col] = "start"
        elif lower == "end":
            col_map[col] = "end"
        elif lower in {"zyg", "zygosity"}:
            col_map[col] = "zygosity"
        elif lower == "sample":
            col_map[col] = "sample"
        elif lower in {"te", "family"}:
            col_map[col] = "te"

    df = df.rename(columns=col_map)
    return df


def normalize_zyg(val: float) -> float:
    try:
        z = float(val)
    except (TypeError, ValueError):
        return 0.0
    return z


def map_to_classifier(zyg: float) -> int:
    if zyg >= 0.75:
        return 2
    if zyg >= 0.25:
        return 1
    return 0


def group_truth(
    truth_df: pd.DataFrame, match_te: bool
) -> Dict[Tuple[str, Optional[str]], List[Tuple[int, int, float]]]:
    truth_groups: Dict[Tuple[str, Optional[str]], List[Tuple[int, int, float]]] = {}
    if match_te:
        grouped = truth_df.groupby(["seqnames", "te"], dropna=False)
    else:
        grouped = truth_df.groupby(["seqnames"], dropna=False)

    for key, group in grouped:
        if isinstance(key, tuple):
            chrom = key[0]
            te = key[1] if len(key) > 1 else None
        else:
            chrom = key
            te = None
        intervals = (
            group[["start", "end", "zygosity"]]
            .sort_values(["start", "end"])
            .itertuples(index=False, name=None)
        )
        truth_groups[(chrom, te)] = [
            (int(s), int(e), normalize_zyg(z)) for s, e, z in intervals
        ]
    return truth_groups


def label_candidates(
    candidates: pd.DataFrame,
    truth_groups: Dict[Tuple[str, Optional[str]], List[Tuple[int, int, float]]],
    match_te: bool,
) -> pd.Series:
    labels = pd.Series(0.0, index=candidates.index, dtype=float)
    key_cols = ["Chrom"] + (["TE"] if match_te else [])

    for key, group in candidates.groupby(key_cols, dropna=False):
        if match_te:
            if isinstance(key, tuple):
                chrom = key[0]
                te = key[1] if len(key) > 1 else None
            else:
                chrom = key
                te = None
        else:
            chrom = key[0] if isinstance(key, tuple) else key
            te = None
        truth_key = (chrom, te)
        truth_intervals = truth_groups.get(truth_key, [])
        if not truth_intervals:
            continue

        truth_intervals = sorted(truth_intervals, key=lambda x: (x[0], x[1]))
        group_sorted = group.sort_values(["Ref_begin", "Ref_end"])
        t_idx = 0
        for idx, row in group_sorted.iterrows():
            c_start = int(row["Ref_begin"])
            c_end = int(row["Ref_end"])
            max_zyg = 0.0
            while t_idx < len(truth_intervals) and truth_intervals[t_idx][1] <= c_start:
                t_idx += 1
            scan_idx = t_idx
            while scan_idx < len(truth_intervals) and truth_intervals[scan_idx][0] < c_end:
                t_start, t_end, t_zyg = truth_intervals[scan_idx]
                if t_end > c_start:
                    if t_zyg > max_zyg:
                        max_zyg = t_zyg
                scan_idx += 1
            labels.at[idx] = max_zyg
    return labels


def main() -> None:
    args = parse_args()

    candidates = pd.read_csv(args.candidates)
    required = {"Sample", "Chrom", "Ref_begin", "Ref_end", "Class", "TE"}
    missing = required.difference(candidates.columns)
    if missing:
        raise ValueError(f"Candidates CSV missing columns: {sorted(missing)}")

    truth_df = read_truth_bed(args.truth)
    if "seqnames" not in truth_df.columns or "start" not in truth_df.columns or "end" not in truth_df.columns:
        raise ValueError("Truth BED must include seqnames, start, end columns.")
    if "zygosity" not in truth_df.columns:
        raise ValueError("Truth BED must include zygosity column.")

    sample_values = candidates["Sample"].unique()
    if len(sample_values) == 1:
        sample_name = sample_values[0]
    else:
        sample_name = None

    if "sample" in truth_df.columns and sample_name is not None:
        truth_df = truth_df[truth_df["sample"] == sample_name]

    match_te = "te" in truth_df.columns
    truth_groups = group_truth(truth_df, match_te=match_te)
    labels = label_candidates(candidates, truth_groups, match_te=match_te)

    if args.label_mode == "classifier":
        candidates["Class"] = labels.apply(map_to_classifier)
    else:
        candidates["Class"] = labels

    candidates.to_csv(args.output, index=False)
    class_counts = candidates["Class"].value_counts().to_dict()
    print(f"[label_candidate_regions_training] Labeled {len(candidates)} regions: {class_counts}")


if __name__ == "__main__":
    main()
