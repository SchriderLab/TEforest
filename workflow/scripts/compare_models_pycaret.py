#!/usr/bin/env python3
import os
import glob
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pycaret.classification import setup, create_model, predict_model, plot_model
from sklearn.metrics import classification_report, f1_score

# ------------------------------
# Parameters and Paths
# ------------------------------

TRAIN_DIR = "/nas/longleaf/home/adaigle/work/test_TEforest/fullnorm_feats_50X/3L3RX"
VALID_DIR = "/nas/longleaf/home/adaigle/work/test_TEforest/fullnorm_feats_50X/2L2R"
PLOTS_DIR = "/nas/longleaf/home/adaigle/TEforest/plots/pycaret_comparisons"
os.makedirs(PLOTS_DIR, exist_ok=True)

# Cache file paths for faster loading
TRAIN_CACHE_FILE = "/nas/longleaf/home/adaigle/TEforest/train_df.pkl.gz"
VALID_CACHE_FILE = "/nas/longleaf/home/adaigle/TEforest/valid_df.pkl.gz"

# Specify the list of model IDs to test
models_to_test = ['lr', 'knn', 'nb', 'dt', 'svm', 'mlp', 
                  'ridge', 'rf', 'qda', 'ada', 'et', 'xgboost', 
                  'lightgbm', 'catboost', 'dummy']

# ------------------------------
# Feature Definitions
# ------------------------------

feature_names = [
    "Cigar1", "Cigar2", "Cigar3", "Cigar4", "Cigar5",
    "Paired", "Proper_Pair", "Is_Read1_Unmapped", "Is_Read2_Unmapped",
    "Is_Read1_Rev_Comp", "Is_Read2_Rev_Comp", "Is_First_Read", "Is_Second_Read",
    "Split", "Long_Insert", "Short_Insert", "Parallel_Read",
    "Everted_Read", "Orphan_Read", "Insert_Size", "Quality",
]

base_features = []
for f in feature_names:
    base_features.extend([f"{f}_mean", f"{f}_median", f"{f}_sd", f"{f}_IQR"])
feature_list_extended = base_features + [f"TE_specific_{feat}" for feat in base_features]

# ------------------------------
# Helper Functions
# ------------------------------

def parse_filename(filepath):
    """
    Parse the filename expected to be in the format:
    sample-chrom-start-stop-genotype-TE.npy
    Uses rsplit with maxsplit=5 to properly handle sample names that include hyphens.
    """
    base = os.path.basename(filepath).replace('.npy', '')
    parts = base.rsplit('-', maxsplit=5)
    if len(parts) != 6:
        raise ValueError(f"Filename {base} does not conform to expected format.")
    sample = parts[0]
    chrom = parts[1]
    try:
        start = int(parts[2])
        stop = int(parts[3])
        genotype = int(parts[4])
    except ValueError as ve:
        raise ValueError(f"Error parsing numerical values from {base}: {ve}")
    region_length = stop - start
    return {"sample": sample, "chrom": chrom, "start": start, "stop": stop,
            "genotype": genotype, "region_length": region_length}

def load_data_from_dir(directory):
    """
    Loads all .npy files from a directory, parses the filename, and constructs a DataFrame.
    """
    data_records = []
    file_pattern = os.path.join(directory, "*.npy")
    for filepath in glob.glob(file_pattern):
        try:
            meta = parse_filename(filepath)
            features = np.load(filepath)
            if len(features) != len(feature_list_extended):
                raise ValueError(f"Feature length mismatch in file {filepath}.")
            record = meta.copy()
            for fname, value in zip(feature_list_extended, features):
                record[fname] = value
            data_records.append(record)
        except Exception as e:
            print(f"Error processing {filepath}: {e}")
    return pd.DataFrame(data_records)

# ------------------------------
# Data Loading with Caching
# ------------------------------

if os.path.exists(TRAIN_CACHE_FILE):
    print("Loading cached training data from:", TRAIN_CACHE_FILE)
    train_df = pd.read_pickle(TRAIN_CACHE_FILE, compression="gzip")
else:
    print("Loading training data from:", TRAIN_DIR)
    train_df = load_data_from_dir(TRAIN_DIR)
    print("Saving training data to cache...")
    train_df.to_pickle(TRAIN_CACHE_FILE, compression="gzip")
print(f"Training data loaded: {train_df.shape[0]} records.")
print("Unique samples in training data:", train_df["sample"].unique())

if os.path.exists(VALID_CACHE_FILE):
    print("Loading cached validation data from:", VALID_CACHE_FILE)
    valid_df = pd.read_pickle(VALID_CACHE_FILE, compression="gzip")
else:
    print("Loading validation data from:", VALID_DIR)
    valid_df = load_data_from_dir(VALID_DIR)
    print("Saving validation data to cache...")
    valid_df.to_pickle(VALID_CACHE_FILE, compression="gzip")
print(f"Validation data loaded: {valid_df.shape[0]} records.")

# ------------------------------
# PyCaret Setup for Classification
# ------------------------------

features_for_training = feature_list_extended + ['region_length']
print("Setting up PyCaret experiment ...")
# Use fold=2 (minimal CV) so that internal CV overhead is low
exp_clf = setup(
    data=train_df,
    target='genotype',
    ignore_features=['sample', 'chrom', 'start', 'stop'],
    numeric_features=features_for_training,
    fix_imbalance=True,
    session_id=42,
    fold=2,
    verbose=False,
    use_gpu=True
)

print("Models to test:", models_to_test)

# Update sample list (correcting JUT-009 to JUT-008)
sample_list = ["A2_A3", "AKA-017_GIM-024", "JUT-008_MUN-009"]

# Dictionaries to store results and sample-specific F1 scores
results_summary = {}
sample_scores = {}  # sample_scores[model_id][sample] = f1 score

# ------------------------------
# Benchmarking Loop Over Specified Models
# ------------------------------

for model_id in models_to_test:
    print("\n====================================")
    print(f"Evaluating model: {model_id}")
    try:
        try:
            model = create_model(model_id, use_gpu=True)
        except TypeError:
            model = create_model(model_id)
        
        # Evaluate on validation data
        valid_preds = predict_model(model, data=valid_df)
        if "Label" in valid_preds.columns:
            pred_col_val = "Label"
        elif "prediction" in valid_preds.columns:
            pred_col_val = "prediction"
        elif "prediction_label" in valid_preds.columns:
            pred_col_val = "prediction_label"
        else:
            print("Debug: Validation predictions output:")
            print(valid_preds.head())
            print("Columns available:", valid_preds.columns)
            raise ValueError("Prediction column not found in validation predictions.")
        
        y_true_valid = valid_preds["genotype"].astype(int)
        y_pred_valid = valid_preds[pred_col_val].astype(int)
        overall_f1 = f1_score(y_true_valid, y_pred_valid, average='macro')
        overall_report = classification_report(y_true_valid, y_pred_valid, digits=3)
        
        print("Overall F1 on validation set: {:.3f}".format(overall_f1))
        print("Classification report:\n", overall_report)
        
        # Subset metric: cases where true in [1, 2] and prediction != 0
        mask = y_true_valid.isin([1, 2]) & (y_pred_valid != 0)
        if mask.sum() > 0:
            subset_true = y_true_valid[mask]
            subset_pred = y_pred_valid[mask]
            subset_f1 = f1_score(subset_true, subset_pred, average='macro')
            subset_report = classification_report(subset_true, subset_pred, digits=3)
            print("Subset F1 (classes 1 & 2, pred ≠ 0): {:.3f}".format(subset_f1))
            print("Subset classification report:\n", subset_report)
        else:
            subset_f1 = None
            subset_report = "N/A"
            print("No subset cases; skipping subset metric.")
        
        # Evaluate per-sample performance on training data
        train_preds = predict_model(model, data=train_df)
        if "Label" in train_preds.columns:
            pred_col_train = "Label"
        elif "prediction" in train_preds.columns:
            pred_col_train = "prediction"
        elif "prediction_label" in train_preds.columns:
            pred_col_train = "prediction_label"
        else:
            print("Debug: Training predictions output:")
            print(train_preds.head())
            print("Columns available:", train_preds.columns)
            raise ValueError("Prediction column not found in training predictions.")
        
        print("\nPerformance on specified training samples:")
        sample_scores[model_id] = {}
        for sample in sample_list:
            sample_mask = train_preds["sample"] == sample
            if sample_mask.sum() > 0:
                sample_true = train_preds.loc[sample_mask, "genotype"].astype(int)
                sample_pred = train_preds.loc[sample_mask, pred_col_train].astype(int)
                sample_f1 = f1_score(sample_true, sample_pred, average='macro')
                sample_scores[model_id][sample] = sample_f1
                sample_report = classification_report(sample_true, sample_pred, digits=3)
                print(f"Sample: {sample} | F1: {sample_f1:.3f}")
                print(sample_report)
            else:
                print(f"Sample {sample} not found in training data.")
        
        # Save the confusion matrix plot with the model ID in its filename.
        old_dir = os.getcwd()
        os.chdir(PLOTS_DIR)
        plot_model(model, plot='confusion_matrix', save=True)
        time.sleep(1)  # Wait for file to be written
        os.chdir(old_dir)
        default_filename = os.path.join(PLOTS_DIR, "Confusion_Matrix.png")
        new_filename = os.path.join(PLOTS_DIR, f"{model_id}_Confusion_Matrix.png")
        if os.path.exists(default_filename):
            os.rename(default_filename, new_filename)
        print(f"Confusion matrix plot saved to: {new_filename}")
        
        results_summary[model_id] = {
            "overall_f1": overall_f1,
            "overall_report": overall_report,
            "subset_f1": subset_f1,
            "subset_report": subset_report,
        }
    except Exception as e:
        print(f"Error evaluating model {model_id}: {e}")

# ------------------------------
# Summary Plots for F1 Scores with Annotations and Ranking
# ------------------------------

# Overall F1 summary plot (sorted descending)
sorted_models = sorted(results_summary.keys(), key=lambda m: results_summary[m]["overall_f1"], reverse=True)
overall_f1_values = [results_summary[m]["overall_f1"] for m in sorted_models]

plt.figure(figsize=(10,6))
bars = plt.bar(sorted_models, overall_f1_values, color="skyblue")
plt.xticks(rotation=45, ha="right")
plt.ylabel("Overall F1 (Macro)")
plt.title("Overall F1 Comparison by Model (Best First)")
for bar in bars:
    height = bar.get_height()
    plt.annotate(f"{height:.3f}",
                 xy=(bar.get_x() + bar.get_width()/2, height),
                 xytext=(0, 3),
                 textcoords="offset points",
                 ha="center", va="bottom")
plt.tight_layout()
overall_plot_path = os.path.join(PLOTS_DIR, "Overall_F1_Comparison.png")
plt.savefig(overall_plot_path)
plt.close()
print(f"Overall F1 comparison plot saved to: {overall_plot_path}")

# Subset F1 summary plot (only for models with subset data)
models_with_subset = [m for m in sorted_models if results_summary[m]["subset_f1"] is not None]
if models_with_subset:
    subset_f1_values = [results_summary[m]["subset_f1"] for m in models_with_subset]
    plt.figure(figsize=(10,6))
    bars = plt.bar(models_with_subset, subset_f1_values, color="lightgreen")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Subset F1 (Macro)")
    plt.title("Subset F1 Comparison by Model (Best First)")
    for bar in bars:
        height = bar.get_height()
        plt.annotate(f"{height:.3f}",
                     xy=(bar.get_x() + bar.get_width()/2, height),
                     xytext=(0, 3),
                     textcoords="offset points",
                     ha="center", va="bottom")
    plt.tight_layout()
    subset_plot_path = os.path.join(PLOTS_DIR, "Subset_F1_Comparison.png")
    plt.savefig(subset_plot_path)
    plt.close()
    print(f"Subset F1 comparison plot saved to: {subset_plot_path}")
else:
    print("No subset F1 data available for plotting.")

# Grouped bar chart for sample-specific F1 (all samples together)
# First, compute the average sample F1 for each model to sort models.
avg_sample_f1 = {}
for m in sample_scores:
    values = [sample_scores[m][s] for s in sample_list if s in sample_scores[m]]
    if values:
        avg_sample_f1[m] = np.mean(values)
    else:
        avg_sample_f1[m] = np.nan

# Sort models by average sample F1 descending
sorted_models_sample = sorted(avg_sample_f1.keys(), key=lambda m: avg_sample_f1[m], reverse=True)

x = np.arange(len(sorted_models_sample))
width = 0.2  # width for each sample's bar
plt.figure(figsize=(12,8))
for i, sample in enumerate(sample_list):
    sample_f1_vals = [sample_scores[m].get(sample, np.nan) for m in sorted_models_sample]
    bars = plt.bar(x + i*width, sample_f1_vals, width, label=sample)
    # Annotate each bar with its F1 value
    for bar in bars:
        height = bar.get_height()
        plt.annotate(f"{height:.3f}",
                     xy=(bar.get_x() + bar.get_width()/2, height),
                     xytext=(0, 3),
                     textcoords="offset points",
                     ha="center", va="bottom")
plt.xticks(x + width, sorted_models_sample, rotation=45, ha="right")
plt.ylabel("Sample-specific F1 (Macro)")
plt.title("Sample-specific F1 Comparison by Model (Best Average First)")
plt.legend()
plt.tight_layout()
sample_plot_path = os.path.join(PLOTS_DIR, "Grouped_Sample_F1_Comparison.png")
plt.savefig(sample_plot_path)
plt.close()
print(f"Grouped sample-specific F1 comparison plot saved to: {sample_plot_path}")
