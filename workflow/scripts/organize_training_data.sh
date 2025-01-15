# this script splits training data based on chromosome.
 
while getopts s: flag
do
    case "${flag}" in
        s) sample=${OPTARG};;
    esac
done

echo "Current directory: $(pwd)"
featvec_dir="$(pwd)/feature_vectors_het"
training_folder=$(pwd)/2L2R
validation_folder=$(pwd)/3L3RX

[ -d $training_folder ] && echo "Output directory exists" || mkdir $training_folder
[ -d $validation_folder ] && echo "Output directory exists" || mkdir $validation_folder

for dir in feature_vectors_het/* #[a-z][0-9][a-z][0-9]
    do echo $dir
    # Copying files to training folder, excluding *readcount* files
    find $dir -type f -name "*2L*" -not -name "*readcount*" | xargs -I {} cp {} $training_folder/
    find $dir -type f -name "*2R*" -not -name "*readcount*" | xargs -I {} cp {} $training_folder/
    touch $training_folder/$sample.txt
    # Copying files to validation folder, excluding *readcount* files
    find $dir -type f -name "*3L*" -not -name "*readcount*" | xargs -I {} cp {} $validation_folder/
    find $dir -type f -name "*3R*" -not -name "*readcount*" | xargs -I {} cp {} $validation_folder/
    find $dir -type f -name "*X*" -not -name "*readcount*" | xargs -I {} cp {} $validation_folder/
    touch $validation_folder/$sample.txt
done
