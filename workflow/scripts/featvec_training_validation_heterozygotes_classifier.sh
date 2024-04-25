#!/bin/bash

# A pipeline to take complete synthetic heterozygotes, train a random forest with them, 
# and validate on specific chromosomes (2L and 2R)
#make sure allseeingeye is conda env.
while getopts t:c:s:v: flag
do
    case "${flag}" in
        t) training_script=${OPTARG};;
        c) condense_script=${OPTARG};;
        s) sample=${OPTARG};;
        v) validation_script=${OPTARG};;
    esac
done

echo "Current directory: $(pwd)"
featvec_dir="$(pwd)/feature_vectors_het"
training_folder=$(pwd)/2L2R
validation_folder=$(pwd)/3L3RX
training_folder_classifer=$(pwd)/2L2R_classifer
validation_folder_classifer=$(pwd)/3L3RX_classifer
#cd $featvec_dir
[ -d $training_folder_classifer ] && echo "Output directory exists" || mkdir $training_folder_classifer
[ -d $validation_folder_classifer ] && echo "Output directory exists" || mkdir $validation_folder_classifer

#train model

python $training_script -i $training_folder/condense.npz -o $training_folder
python $training_script -i $validation_folder/condense.npz -o $validation_folder

#validate
python $validation_script -n $validation_folder/condense.npz -m $training_folder/svrf_classifier_all.pkl -o $training_folder_classifer
python $validation_script -n $training_folder/condense.npz -m $validation_folder/svrf_classifier_all.pkl -o $validation_folder_classifer

for dir in "$featvec_dir/*" #[a-z][0-9][a-z][0-9]
    do echo $dir
    touch $validation_folder_classifer/$sample.txt
    touch $training_folder_classifer/$sample.txt
done
