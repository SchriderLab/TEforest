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
#cd $featvec_dir
[ -d $training_folder ] && echo "Output directory exists" || mkdir $training_folder
[ -d $validation_folder ] && echo "Output directory exists" || mkdir $validation_folder

for dir in "$featvec_dir/*" #[a-z][0-9][a-z][0-9]
    do echo $dir
    find $dir -type f -name "*2L*" -exec cp {} $training_folder/ \; 
    find $dir -type f -name "*2R*" -exec cp {} $training_folder/ \; 
    touch $training_folder/$sample.txt
    find $dir -type f -name "*3L*" -exec cp {} $validation_folder/ \; 
    find $dir -type f -name "*3R*" -exec cp {} $validation_folder/ \; 
    find $dir -type f -name "*X*" -exec cp {} $validation_folder/ \; 
    touch $validation_folder/$sample.txt
done

#train model
python $condense_script -i $training_folder -o $training_folder/condense.npz
python $condense_script -i $validation_folder -o $validation_folder/condense.npz

python $training_script -i $training_folder/condense.npz -o $training_folder
python $training_script -i $validation_folder/condense.npz -o $validation_folder

#validate
python $validation_script -n $validation_folder/condense.npz -m $training_folder/svrf_all.pkl -o $training_folder
python $validation_script -n $training_folder/condense.npz -m $validation_folder/svrf_all.pkl -o $validation_folder
