#!/usr/bin/env python3

import argparse
import os
import pickle as pkl

import numpy as np
from lightgbm import LGBMClassifier, LGBMRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, mean_squared_error


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train a LightGBM TEforest model from condensed feature vectors."
    )
    parser.add_argument(
        "--input-npz",
        required=True,
        help="Condensed .npz produced by condense_training_data.py",
    )
    parser.add_argument(
        "--model-out",
        required=True,
        help="Output path for the trained model (.pkl)",
    )
    parser.add_argument(
        "--mode",
        choices=["classifier", "regressor"],
        default="classifier",
        help="Train classifier (default) or regressor.",
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=0.0,
        help="Fraction to hold out for testing (0 disables split).",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed for train/test split.",
    )
    parser.add_argument(
        "--n-estimators",
        type=int,
        default=500,
        help="Number of boosting rounds.",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=8,
        help="Number of threads for LightGBM.",
    )
    return parser.parse_args()


def load_npz(path: str):
    with np.load(path, mmap_mode="r") as z:
        if "X" not in z.files or "labels" not in z.files:
            raise ValueError(f"Missing X/labels in {path}")
        X = z["X"]
        labels = z["labels"]
    return X, labels


def main() -> None:
    args = parse_args()
    X, labels = load_npz(args.input_npz)
    if X.size == 0:
        raise ValueError(f"No samples found in {args.input_npz}")

    if args.mode == "classifier":
        y = labels.astype(float).round().astype(int)
        classes = np.unique(y)
    else:
        y = labels.astype(float)
        classes = None

    if args.test_size and args.test_size > 0:
        stratify = y if args.mode == "classifier" and len(set(y)) > 1 else None
        X_train, X_test, y_train, y_test = train_test_split(
            X,
            y,
            test_size=args.test_size,
            random_state=args.random_state,
            stratify=stratify,
        )
    else:
        X_train, y_train = X, y
        X_test, y_test = None, None

    if args.mode == "classifier":
        objective = "binary" if len(classes) <= 2 else "multiclass"
        model = LGBMClassifier(
            n_estimators=args.n_estimators,
            random_state=args.random_state,
            class_weight="balanced",
            n_jobs=args.n_jobs,
            objective=objective,
            verbose=-1,
        )
    else:
        model = LGBMRegressor(
            n_estimators=args.n_estimators,
            random_state=args.random_state,
            n_jobs=args.n_jobs,
            objective="regression",
            verbose=-1,
        )

    model.fit(X_train, y_train)

    if X_test is not None:
        preds = model.predict(X_test)
        if args.mode == "classifier":
            acc = accuracy_score(y_test, preds)
            print(f"[train_lightgbm] Accuracy: {acc:.4f}")
        else:
            rmse = mean_squared_error(y_test, preds, squared=False)
            print(f"[train_lightgbm] RMSE: {rmse:.6f}")

    os.makedirs(os.path.dirname(args.model_out), exist_ok=True)
    with open(args.model_out, "wb") as fh:
        pkl.dump(model, fh)
    print(f"[train_lightgbm] Wrote model to {args.model_out}")


if __name__ == "__main__":
    main()
