import os
import pickle as pkl

import numpy as np
import pandas as pd

# import sklearn as sklearn
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.model_selection import train_test_split
from tqdm import tqdm
import pycaret

pycaret.__version__

# infile = "/nas/longleaf/home/adaigle/work/test_TEforest/het_experiments_generate_hetdata/3L3RX/condense.npz"
# infile = "/nas/longleaf/home/adaigle/work/mcclintock_stuff/find_candidate_regions/heterozygosity_experiments/outputs/feature_vectors/condense.npz"
# infile = "/nas/longleaf/home/adaigle/work/test_TEforest/het_experiments_A2A3oldreads/3L3RX/condense.npz"
# infile = "/nas/longleaf/home/adaigle/work/test_TEforest/het_experiments_50bp/feature_vectors_het/A2A3/condense.npz"
# infile = "/nas/longleaf/home/adaigle/work/mcclintock_stuff/find_candidate_regions/heterozygosity_experiments/outputs/feature_vectors/training_not2L/condense.npz"
infile = "/nas/longleaf/home/adaigle/work/test_TEforest/het_experiments_generate_hetdata/condense_alldata/condense.npz"
infile = "/nas/longleaf/home/adaigle/work/test_TEforest/het_experiments_full50bp/3L3RX/condense.npz"
infile = "/nas/longleaf/home/adaigle/work/test_TEforest/het_experiments_full50bp_exp2/3L3RX/condense.npz"

sv_data = np.load(infile)

samps = {"data": [], "labs": [], "files": []}
for id in tqdm(sv_data.files):
    samps["data"].append(sv_data[id])
    samps["labs"].append(id.split("-")[-2])
    samps["files"].append(id)

# preds = loaded_model.predict(samps["data"])

data = {
    "file": samps["files"],
    "true": samps["labs"],
    "data": samps["data"]
    # "pred": preds,
    # "cntrl_score": [0] * len(samps["files"]),
}

df = pd.DataFrame.from_dict(samps)
# df['read_length'] = df['files'].apply(lambda x: 54 if str(x).startswith('A2') else None)

df["read_length"] = df["files"].apply(
    lambda x: 151
    if str(x).startswith("JUT")
    else (
        75
        if str(x).startswith("RAL")
        else (
            125
            if str(x).startswith("AKA")
            else (54 if str(x).startswith("A2") else None)
        )
    )
)
# df['insert_size'] = df['files'].apply(lambda x: 415 if str(x).startswith('JUT') else (220 if str(x).startswith('RAL') else (208 if str(x).startswith('AKA') else (287 if str(x).startswith('A2_A3') else None))))
df["insert_size"] = df["files"].apply(
    lambda x: 415
    if str(x).startswith("JUT")
    else (
        220
        if str(x).startswith("RAL")
        else (
            208
            if str(x).startswith("AKA")
            else (287 if str(x).startswith("A2") else None)
        )
    )
)
# df = df[df['files'].str.startswith('A2')]


full_df = pd.DataFrame(data)
# print(full_df.head())

feature_names = [
    "Cigar1",
    "Cigar2",
    "Cigar3",
    "Cigar4",
    "Cigar5",
    "Paired",
    "Proper_Pair",
    "Is_Read1_Unmapped",
    "Is_Read2_Unmapped",
    "Is_Read1_Rev_Comp",
    "Is_Read2_Rev_Comp",
    "Is_First_Read",
    "Is_Second_Read",
    "Split",
    "Long_Insert",
    "Short_Insert",
    "Parallel_Read",
    "Everted_Read",
    "Orphan_Read",
    "Insert_Size",
    "Quality",
]

feature_list = []
for x in feature_names:
    feature_list.extend([f"{x}_mean", f"{x}_median", f"{x}_sd", f"{x}_IQR"])

feature_list_extended = feature_list + [
    "TE_specific_" + feature for feature in feature_list
]

expanded_data = full_df["data"].apply(pd.Series)

full_df = pd.concat([full_df, expanded_data], axis=1)
full_df.drop(columns=["data", "file"], inplace=True)
full_df.columns = ["true"] + feature_list_extended
full_df = full_df.merge(
    df[["read_length"]], left_index=True, right_index=True, how="left"
)
full_df = full_df.merge(
    df[["insert_size"]], left_index=True, right_index=True, how="left"
)
# full_df = full_df[full_df['read_length'] != 75]
# full_df = full_df[full_df['read_length'] != 151]
# full_df = full_df[full_df['read_length'] != 125]

# full_df = full_df[full_df['read_length'] == 75]

print(full_df.head())
full_df["true"] = pd.to_numeric(full_df["true"], errors="coerce")

# full_df = full_df[full_df['true'] != 0]
# full_df.loc[full_df['true'] == 2, 'true'] = 1

# full_df.to_csv('/nas/longleaf/home/adaigle/TEforest/workflow/scripts/feature_data.csv')
# working through pycaret regression tutorial
# import pycaret regression and init setup

# from pycaret.regression import *
#
# s = setup(full_df, target = 'true', session_id = 123,html=False, use_gpu=True)
#
# best = compare_models(include=['rf', 'lightgbm'])
# print("done training")
# print(best)
#
## Pull the performance metrics
# results = pull()
#
## Print the metrics
# print(results)
from pycaret.classification import *

# full_df = full_df.loc[full_df['true'] != '0']
s = setup(
    full_df, target="true", session_id=123, html=False, use_gpu=True
)  # , fix_imbalance=True)#, normalize = True)

# best = compare_models(exclude = ['gbc'],)
best = compare_models(include=["rf", "lightgbm"])

print("done training")
print(best)

# Pull the performance metrics
results = pull()

# Print the metrics
# print(results)

plot_model(best, plot="class_report", save=True)
# interpret_model(best, plot = 'summary_classification', save=True)

predicted_values = predict_model(best, data=full_df)  # Predicted labels
print(predicted_values.head())
print(type(predicted_values))
# actual_values = predicted_values['true']  # Ground truth
# predicted_values = predicted_values['prediction_label']
## Create a DataFrame to compare actual vs. predicted
# comparison_df = pd.DataFrame({
#    'Actual': actual_values,
#    'Predicted': predicted_values
# })

# Add a column to indicate whether the prediction was correct
predicted_values["Correct"] = (
    predicted_values["true"] == predicted_values["prediction_label"]
)

# Rows where the model got it wrong
incorrect_predictions = predicted_values[predicted_values["Correct"] == False]

# Rows where the model got it right
correct_predictions = predicted_values[predicted_values["Correct"] == True]

# Print the first few rows of each
print("Incorrect Predictions:")
print(incorrect_predictions.head())
print(incorrect_predictions.shape)

print("\nCorrect Predictions:")
print(correct_predictions.head())
print(correct_predictions.shape)
