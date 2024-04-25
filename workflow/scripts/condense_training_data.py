import numpy as np
from glob import glob
from tqdm import tqdm
import os
import argparse


def load_npy(npfile):
    return np.load(npfile)


def extract_name(filename):
    bn = os.path.basename(filename)
    return bn.split(".npy")[0]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_dir",
        metavar="INDIR",
        required=True,
        dest="indir",
        help="Top-level directory containing subdirs/npy files to condense.",
    )

    parser.add_argument(
        "-o",
        "--output-npz",
        metavar="OUTFILE",
        required=True,
        dest="outfile",
        help="NPZ file to write condensed data to.",
    )

    return parser.parse_args()


def main():
    ap = parse_args()

    data_dict = {}
    for nfile in tqdm(
        glob(f"{ap.indir}/**/*.npy", recursive=True), desc="Condensing data"
    ):
        data_dict[extract_name(nfile)] = load_npy(nfile)

    np.savez(ap.outfile, **data_dict)


if __name__ == "__main__":
    main()
