#!/usr/bin/env python3
import os
import glob
import time
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shap

from pycaret.classification import setup, create_model, tune_model, finalize_model, predict_model, plot_model
from sklearn.metrics import classification_report, f1_score
from lightgbm import LGBMClassifier
import catboost

# ------------------------------
# Parameters and Paths
# ------------------------------

# Data directories (training/validation from previous experiments)
TRAIN_DIR = "/nas/longleaf/home/adaigle/work/test_TEforest/fullnorm_feats_50X/3L3RX"
VALID_DIR = "/nas/longleaf/home/adaigle/work/test_TEforest/fullnorm_feats_50X/2L2RX"  # adjust if needed

# Directory for saving plots
PLOTS_DIR = "/nas/longleaf/home/adaigle/TEforest/plots/pycaret_comparisons"
os.makedirs(PLOTS_DIR, exist_ok=True)

# Directory for saving models
MODELS_DIR = "/nas/longleaf/home/adaigle/TEforest/models"
os.makedirs(MODELS_DIR, exist_ok=True)
#BEST_MODEL_PATH = os.path.join(MODELS_DIR, "best_lightgbm_model.pkl")
#BEST_MODEL_PATH = os.path.join(MODELS_DIR, "best_rf_model.pkl")
BEST_MODEL_PATH = os.path.join(MODELS_DIR, "best_catboost_model.pkl")

# Directory for unlabeled prediction data
#UNLABELED_DIR = "/nas/longleaf/home/adaigle/work/test_TEforest/test_celegans_fullnorm_newrefte/feature_vectors/RW7000"
UNLABELED_DIR = "/nas/longleaf/home/adaigle/users/rob_flies/teforest_runs/fullMA_fullnorm/feature_vectors/MA19F"

# Cache file paths for faster data loading
TRAIN_CACHE_FILE = "/nas/longleaf/home/adaigle/TEforest/train_df.pkl.gz"
VALID_CACHE_FILE = "/nas/longleaf/home/adaigle/TEforest/valid_df.pkl.gz"

# For SHAP force/waterfall plots: list of example filenames (from the unlabeled set)
#selected_examples = ["RW7000-V-17711022-17711990-0-Tc1.npy"]
selected_examples = ["MA19F-2R-17442101-17442703-0-297.npy","MA19F-2L-12229266-12230042-0-BS.npy","MA19F-3R-12786866-12787575-0-297.npy","MA19F-X-10130584-10131127-0-roo.npy"]

# ------------------------------
# Feature Definitions (same as before)
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

# The features to be used for training (excluding metadata)
features_for_training = feature_list_extended + ["region_length"]

# ------------------------------
# Helper Functions
# ------------------------------

def parse_filename(filepath):
    """
    Parse the filename expected to be in the format:
    sample-chrom-start-stop-genotype-TE.npy
    Uses rsplit with maxsplit=5 to handle sample names with hyphens.
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
    Load all .npy files from a directory, parse filenames, and build a DataFrame.
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
# Baseline Model (Default Parameters)
# ------------------------------

print("\n=== Training Baseline LightGBM Model (Default Parameters) ===")
exp_baseline = setup(
    data=train_df,
    target="genotype",
    ignore_features=["sample", "chrom", "start", "stop"],
    numeric_features=features_for_training,
    fix_imbalance=True,
    session_id=42,
    fold=2,
    verbose=False#,
    #use_gpu=True
)

baseline_model = create_model("catboost")#, use_gpu=True)
baseline_model = finalize_model(baseline_model)
baseline_preds = predict_model(baseline_model, data=valid_df)
baseline_y_true = baseline_preds["genotype"].astype(int)
# Determine prediction column (usually "Label")
if "Label" in baseline_preds.columns:
    baseline_y_pred = baseline_preds["Label"].astype(int)
elif "prediction" in baseline_preds.columns:
    baseline_y_pred = baseline_preds["prediction"].astype(int)
else:
    baseline_y_pred = baseline_preds["prediction_label"].astype(int)
baseline_f1 = f1_score(baseline_y_true, baseline_y_pred, average="macro")
print(f"Baseline Model Validation Macro F1: {baseline_f1:.3f}")

# ------------------------------
# Hyperparameter Tuning for LightGBM via PyCaret
# ------------------------------
print("\n=== Tuning LightGBM Model for Best Performance ===")
exp_tune = setup(
    data=train_df,
    target="genotype",
    ignore_features=["sample", "chrom", "start", "stop"],
    numeric_features=features_for_training,
    fix_imbalance=True,
    session_id=42,
    fold=10,
    verbose=False#,
  #  use_gpu=True  # Keep GPU enabled for speed (warnings about GPU and sparse data may still occur)
)

if os.path.exists(BEST_MODEL_PATH):
    print("Loading saved best LightGBM model from:", BEST_MODEL_PATH)
    with open(BEST_MODEL_PATH, "rb") as f:
        best_model = pickle.load(f)
else:
    print("Tuning LightGBM model ...")
    model_lgb = create_model("catboost")#, use_gpu=True)
    
    # Define a custom grid that expands the search space:
    custom_grid = {
        "num_leaves": [31, 63, 127],
        "min_data_in_leaf": [10, 20, 50],  # Lower than default to address "no meaningful features" warning
        "learning_rate": [0.01, 0.05, 0.1],
        "feature_fraction": [0.6, 0.7, 0.8],
        "bagging_fraction": [0.5, 0.6, 0.7],
        "bagging_freq": [3, 5],
        # Note: do not include use_gpu in the grid; GPU usage is already set in setup/create_model
    }
    
    tuned_model = tune_model(model_lgb, optimize="F1")# custom_grid=custom_grid, optimize="F1")
    best_model = finalize_model(tuned_model)
    with open(BEST_MODEL_PATH, "wb") as f:
        pickle.dump(best_model, f)
    print("Best model saved to:", BEST_MODEL_PATH)

# Evaluate tuned model on validation set
tuned_preds = predict_model(best_model, data=valid_df)
tuned_y_true = tuned_preds["genotype"].astype(int)
if "Label" in tuned_preds.columns:
    tuned_y_pred = tuned_preds["Label"].astype(int)
elif "prediction" in tuned_preds.columns:
    tuned_y_pred = tuned_preds["prediction"].astype(int)
else:
    tuned_y_pred = tuned_preds["prediction_label"].astype(int)
tuned_f1 = f1_score(tuned_y_true, tuned_y_pred, average="macro")
print(f"Tuned Model Validation Macro F1: {tuned_f1:.3f}")
print(f"Improvement from tuning: {tuned_f1 - baseline_f1:.3f}")

# ------------------------------
# Training Data Size Experiment (Expanded)
# ------------------------------

# Define training and validation features and labels
X_train = train_df[features_for_training]
y_train = train_df["genotype"]
X_val = valid_df[features_for_training]
y_val = valid_df["genotype"]

# Retrieve best model hyperparameters and filter out problematic ones
best_params = best_model.get_params()
allowed_types = (int, float, str, bool, type(None))
filtered_params = {k: v for k, v in best_params.items() if isinstance(v, allowed_types)}
if "verbose" in filtered_params and isinstance(filtered_params["verbose"], bool):
    filtered_params["verbose"] = int(filtered_params["verbose"])

# Use more fractions: 10%, 25%, 50%, 75%, 100%
fractions = [0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.98, 0.99, 1.0]
f1_scores_by_fraction = []

print("\n=== Running Training Data Size Experiment ===")
for frac in fractions:
    sample_df = train_df.sample(frac=frac)
    X_sample = sample_df[features_for_training]
    y_sample = sample_df["genotype"]
    
    model_instance = LGBMClassifier(**filtered_params)
    model_instance.fit(X_sample, y_sample)
    
    y_val_pred = model_instance.predict(X_val)
    f1_val = f1_score(y_val, y_val_pred, average="macro")
    f1_scores_by_fraction.append(f1_val)
    print(f"Training fraction {frac*100:.0f}%: Validation Macro F1 = {f1_val:.3f}")

plt.figure(figsize=(8,6))
plt.plot([frac*100 for frac in fractions], f1_scores_by_fraction, marker="o", linestyle="--", color="blue")
plt.xlabel("Percentage of Training Data")
plt.ylabel("Validation Macro F1")
plt.title("Effect of Training Data Size on Model F1 Score")
plt.grid(True)
plt.tight_layout()
data_size_plot_path = os.path.join(PLOTS_DIR, "Training_Data_Size_vs_F1catboost.png")
plt.savefig(data_size_plot_path)
plt.close()
print("Training data size experiment plot saved to:", data_size_plot_path)

# ------------------------------
# Predictions on Unlabeled Data and SHAP Analysis
# ------------------------------

print("\n=== Predictions on Unlabeled Data and SHAP Analysis ===")
print("Loading unlabeled data from:", UNLABELED_DIR)
unlabeled_df = load_data_from_dir(UNLABELED_DIR)
print(f"Unlabeled data loaded: {unlabeled_df.shape[0]} records.")

# Make predictions using the tuned best model
unlabeled_preds = best_model.predict(unlabeled_df[features_for_training])
unlabeled_df["predicted_genotype"] = unlabeled_preds
predictions_csv_path = os.path.join(PLOTS_DIR, "Unlabeled_Predictions.csv")
unlabeled_df.to_csv(predictions_csv_path, index=False)
print("Unlabeled predictions saved to:", predictions_csv_path)

print("\nComputing SHAP values for overall summary plot ...")
# Extract raw model from PyCaret pipeline if needed
if hasattr(best_model, "steps"):
    raw_model = best_model.steps[-1][1]
else:
    raw_model = best_model

# Subsample training data if large
X_train_sample = X_train.sample(n=min(1000, len(X_train)), random_state=42)
explainer = shap.Explainer(raw_model, X_train_sample)
shap_values = explainer(X_train_sample)

# Generate per-class SHAP summary plots if multiclass
if len(shap_values.values.shape) == 3:
    num_classes = shap_values.values.shape[2]
    for class_idx in range(num_classes):
        plt.figure()
        shap.summary_plot(shap_values[:, :, class_idx], X_train_sample, show=False)
        plot_path = os.path.join(PLOTS_DIR, f"SHAP_Summary_Class{class_idx}catboost.png")
        plt.savefig(plot_path, bbox_inches="tight")
        plt.close()
        print(f"SHAP summary plot for class {class_idx} saved to: {plot_path}")
else:
    print("SHAP output is not multiclass. Skipping per-class plots.")

# Global SHAP summary plot
plt.figure()
shap.summary_plot(shap_values, X_train_sample, show=False)
summary_plot_path = os.path.join(PLOTS_DIR, "SHAP_Overall_Summarycatboost.png")
plt.savefig(summary_plot_path, bbox_inches="tight")
plt.close()
print(f"Global SHAP summary plot saved to: {summary_plot_path}")

# ------------------------------
# SHAP Waterfall Plots for Selected Unlabeled Examples
# ------------------------------
print("\nGenerating SHAP waterfall plots for selected unlabeled examples ...")
for ex in selected_examples:
    ex_path = os.path.join(UNLABELED_DIR, ex)
    if os.path.exists(ex_path):
        try:
            meta = parse_filename(ex_path)
            features = np.load(ex_path)
            if len(features) != len(feature_list_extended):
                print(f"Feature length mismatch in file {ex}")
                continue
            record = meta.copy()
            for fname, value in zip(feature_list_extended, features):
                record[fname] = value
            ex_df = pd.DataFrame([record])
            X_ex = ex_df[features_for_training]

            # Predict class using raw_model (extracted earlier)
            predicted_class = raw_model.predict(X_ex)[0]
            print(f"Example: {ex} | Predicted class: {predicted_class}")

            # Explain the prediction
            shap_ex = explainer(X_ex)
            if len(shap_ex.values.shape) == 3:
                # For multiclass, select values for the predicted class
                shap_values_instance = shap_ex.values[0, :, predicted_class]
                base_value = explainer.expected_value[predicted_class]
            else:
                shap_values_instance = shap_ex.values[0]
                base_value = explainer.expected_value

            # Generate a waterfall plot for the prediction
            plt.figure()
            shap.plots.waterfall(
                shap.Explanation(
                    values=shap_values_instance,
                    base_values=base_value,
                    data=X_ex.iloc[0],
                    feature_names=X_ex.columns
                ),
                max_display=20
            )
            waterfall_plot_path = os.path.join(PLOTS_DIR, f"SHAP_Waterfall_{os.path.splitext(ex)[0]}catboost.png")
            plt.savefig(waterfall_plot_path, bbox_inches="tight")
            plt.close()
            print(f"Waterfall plot for {ex} saved to: {waterfall_plot_path}")

        except Exception as e:
            print(f"Failed to generate waterfall plot for {ex}: {e}")
    else:
        print(f"Example file not found: {ex_path}")

print("\nExperiment completed. All plots and models have been saved.")
