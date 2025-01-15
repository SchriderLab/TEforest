#!/bin/bash

# this script downsamples training data based on class.
training_folder=$(pwd)/2L2R_reference
validation_folder=$(pwd)/3L3RX_reference

downsample_files() {
    folder=$1
    # Get the minimum number of files matching each pattern
    min_count=$(find $folder -type f -name "*-0-*" | wc -l)
    min_count_1=$(find $folder -type f -name "*-1-*" | wc -l)
    min_count_2=$(find $folder -type f -name "*-2-*" | wc -l)

    min_count=$(($min_count < $min_count_1 ? $min_count : $min_count_1))
    min_count=$(($min_count < $min_count_2 ? $min_count : $min_count_2))

    # Downsample files to the minimum count
    for pattern in "-0-" "-1-" "-2-"; do
        files=($(find $folder -type f -name "*$pattern*"))
        total_files=${#files[@]}
        if [ $total_files -gt $min_count ]; then
            files_to_remove=$((total_files - min_count))
            for ((i=0; i<$files_to_remove; i++)); do
                rm "${files[$i]}"
            done
        fi
    done
}

# Downsample files in training and validation folders
downsample_files $training_folder
downsample_files $validation_folder