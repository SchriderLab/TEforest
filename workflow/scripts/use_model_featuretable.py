import argparse
import os
import pickle as pkl
import numpy as np
import pandas as pd
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n",
        "--npz",
        required=True,
        dest="npz",
        help="Path to condensed npz data.",
    )
    parser.add_argument(
        "-m",
        "--model",
        required=True,
        dest="model",
        help="Path to model that will make predictions using the data.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        dest="outdir",
        help="Directory to save model outputs to.",
    )
    return parser.parse_args()

def main():
    ap = parse_args()

    # Load the model from file
    with open(ap.model, "rb") as f:
        loaded_model = pkl.load(f)

    # Load the data
    sv_data = np.load(ap.npz)
    samps = {"data": [], "labs": [], "files": []}

    for id in tqdm(sv_data.files):
        samps["data"].append(sv_data[id])
        samps["labs"].append(id.split("-")[-2])  # Extract true label
        samps["files"].append(id)  # File identifier

    # Convert data to NumPy array
    X_data = np.array(samps["data"])
    y_true = samps["labs"]
    file_ids = samps["files"]

    # Predict using the loaded model
    preds = loaded_model.predict(X_data)

    # Create a DataFrame to store predictions, true labels, and features
    features_df = pd.DataFrame(X_data, columns=[f"feature_{i}" for i in range(X_data.shape[1])])
    results_df = pd.DataFrame({
        "file": file_ids,
        "true": y_true,
        "pred": preds,
    })
    full_df = pd.concat([results_df, features_df], axis=1)

    # Save the resulting DataFrame to a CSV file
    os.makedirs(ap.outdir, exist_ok=True)
    output_path = os.path.join(ap.outdir, "predictions_with_features.csv")
    full_df.to_csv(output_path, index=False)

    print(f"Predictions with features saved to {output_path}")

if __name__ == "__main__":
    main()
