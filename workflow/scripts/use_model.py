import argparse
import os
import pickle as pkl
import argparse

import numpy as np
import pandas as pd
import sklearn as sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n",
        "--npz",
        # metavar="INDIR",
        required=True,
        dest="npz",
        help="Path to condensed npz data.",
    )

    # parser.add_argument(
    #    "-c",
    #    "--csv",
    #    #metavar="INDIR",
    #    required=True,
    #    dest="csv",
    #    help="Path to candidate regions csv."
    # )

    parser.add_argument(
        "-m",
        "--model",
        # metavar="INDIR",
        required=True,
        dest="model",
        help="Path to model that will make predictions using the data.",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        metavar="OUTDIR",
        required=True,
        dest="outdir",
        help="Directory to save model outputs to.",
    )

    return parser.parse_args()


def main():
    ap = parse_args()

    # Load model from file
    with open(ap.model, "rb") as f:
        loaded_model = pkl.load(f)

    infile = ap.npz

    sv_data = np.load(infile)
    samps = {"data": [], "labs": [], "files": []}

    for id in tqdm(sv_data.files):
        samps["data"].append(sv_data[id])
        samps["labs"].append(id.split("-")[-2])
        samps["files"].append(id)

    preds = loaded_model.predict(samps["data"])
    #probs = loaded_model.predict_proba(samps["data"])
    # pd.DataFrame(
    #    {
    #        "file": samps["files"],
    #        "true": samps["labs"],
    #        "pred": preds,
    #        "cntrl_score": probs[:, 0],
    #    }
    # ).to_csv(os.path.join(ap.csv), index=False)

    # Create the DataFrame
    data = {
        "file": samps["files"],
        "true": samps["labs"],
        "pred": preds,
        "cntrl_score": [0] * len(samps["files"]),
    }

    full_df = pd.DataFrame(data)

    # note; this doesn't work on genomes with a hyphen in their name.
    # need to work on this before this step to avoid bugs
    def custom_grouping(file_name):
        return file_name.split("-")[0]

    # Group the DataFrame using the custom grouping function
    grouped = full_df.groupby(full_df["file"].apply(custom_grouping))

    # Define the output directory
    output_directory = ap.outdir
    output_path = os.path.join(output_directory, f"predictions.csv")
    # Save the group DataFrame to a CSV file
    full_df.to_csv(output_path, index=False)
    # Iterate over the groups and save each group to a CSV file
    #for group_name, group_df in grouped:
    #    # Define the full path to the output CSV file
    #    print(group_name)
    #    # output_path = os.path.join(output_directory, f"{group_name}.csv")
    #    output_path = os.path.join(output_directory, f"predictions.csv")
    #    # Save the group DataFrame to a CSV file
    #    group_df.to_csv(output_path, index=False)


if __name__ == "__main__":
    main()
